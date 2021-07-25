use std::{
    collections::{HashMap, HashSet},
    f32::consts::PI,
    ops::Not,
};

use itertools::Itertools;

use crate::{
    types::{BoxIdx, CellIdx, CellVec, EdgeIdx, EdgeVec, VertIdx, VertVec},
    Builder, Shape, Symmetry, V2,
};

/* FUNCTIONS TO EXPORT */

/// Converts this `Builder` into a [`Shape`] and the associated [`Symmetry`]
pub fn build(mut bdr: Builder) -> Result<(Shape, Symmetry), BuildError> {
    // Before we start building the `Shape`, we normalise each box so that the vertices are
    // always specified in clockwise order.  This also detects pathological box shapes (i.e.
    // boxes which are builder-intersecting, non-convex or linear).
    normalise_box_direction(&mut bdr)?;

    // Generate edges between boxes
    let (edges, box_edge_map) = generate_edges(&bdr)?;
    // Generate cell vertices
    let (vert_positions, vert_idx_by_coord) = generate_cell_vertices(&bdr, &edges)?;
    // Generate cells
    let (cells, cell_idx_by_coord) = generate_cells(&bdr, &vert_idx_by_coord);
    // Generate groups
    let groups = generate_groups(&bdr, &edges, &cell_idx_by_coord, &box_edge_map);
    // Package into a `Shape`
    let shape = Shape {
        num_symbols: bdr.box_width * bdr.box_height,
        groups,
        cell_verts: cells,
        verts: vert_positions,
    };

    // Generate the symmetry (currently just assume - incorrectly - that there's no symmetry)
    let symmetry = Symmetry::asymmetric(&shape);
    Ok((shape, symmetry))
}

/// The possible ways that converting a [`Builder`] into a [`Shape`] (using [`Builder::build`])
/// could fail.
#[derive(Debug, Clone)]
pub enum BuildError {
    /// Three vertices in a box form a straight line
    LinearBox { box_idx: BoxIdx },
    /// The built shape contains a box that isn't strictly convex
    NonConvexBox { box_idx: BoxIdx },
    /// The built shape contains a box that self-intersects (i.e. is a bow-tie shape)
    SelfIntersectingBox { box_idx: BoxIdx },

    /// An edge has two boxes attached to the same side of it.  Note that edges often have two
    /// boxes attached, but they are attached in opposite directions (because vertices are always
    /// numbered in a consistent direction):
    /// ```text
    /// 0 ----- (1 2) ----- 3
    /// |   ->    |    ->   |
    /// | ^     | | ^     | |
    /// | |     V | |     V |
    /// |   <-    |    <-   |
    /// 3 ----- (2 1) ----- 0
    /// ```
    /// Overlapping boxes would look like this:
    /// ```text
    ///     0 - 1
    ///    / \ / \
    ///   /   X   \
    ///  /   / \   \
    /// 3 - 2   3 - 2
    /// ```
    OverlappingFacesOnEdge {
        box_idx_1: BoxIdx,
        box_idx_2: BoxIdx,
    },
    /// Two boxes on either side of an edge disagree on the number of cells.
    InconsistentEdge {
        box_idx_left: BoxIdx,
        num_cells_left: usize,

        box_idx_right: BoxIdx,
        num_cells_right: usize,
    },
}

/////////////////////////
// GENERATOR FUNCTIONS //
/////////////////////////

/// Reorder the vertices in each box so that the vertices are always specified in clockwise
/// order.
fn normalise_box_direction(bdr: &mut Builder) -> Result<(), BuildError> {
    // List of indices of boxes which need to be flipped (we can't flip them in the loop
    // without tripping the borrow checker)
    let mut anti_clockwise_box_idxs = Vec::<BoxIdx>::new();
    for (box_idx, vert_idxs) in bdr.boxes.indexed_iter() {
        // Iterate over the positions of all consecutive triple of vertices, and classify that
        // corner according to its direction
        let directions = vert_idxs
            .iter()
            .map(|&idx| bdr.verts[idx].position)
            .circular_tuple_windows()
            .map(|(v1, v2, v3)| {
                // The edges go v1 -> v2 and v2 -> v3
                let d1 = v2 - v1;
                let d2 = v3 - v2;
                // Classify the middle corner (v2) according to its angle
                CornerType::classify(angle_between(d1, d2))
            })
            .collect_vec();

        // Count how many of corners have each classification
        let num_anticlockwise = directions
            .iter()
            .filter(|c| **c == CornerType::AntiClockwise)
            .count();
        let num_linear = directions
            .iter()
            .filter(|c| **c == CornerType::Linear)
            .count();
        let num_clockwise = directions
            .iter()
            .filter(|c| **c == CornerType::Clockwise)
            .count();

        match (num_anticlockwise, num_linear, num_clockwise) {
            // If all corners go anticlockwise, then this box is strictly convex but
            // anticlockwise (and will be correct after being flipped)
            (4, 0, 0) => anti_clockwise_box_idxs.push(box_idx),
            // If all corners go clockwise, then this box doesn't need modification
            (0, 0, 4) => (),

            // If two corners go each way, then the box must b-intersect
            (2, 0, 2) => return Err(BuildError::SelfIntersectingBox { box_idx }),
            // If there are no linear corners, then the only remaining case is that the box
            // isn't strictly convex
            (_, 0, _) => return Err(BuildError::NonConvexBox { box_idx }),
            // If there are any linear corners, then we class the whole box as linear
            (_, _, _) => {
                assert_ne!(num_linear, 0);
                return Err(BuildError::LinearBox { box_idx });
            }
        }
    }
    // Reverse all the anti-clockwise boxes (by swapping the top-left and bottom-right
    // corners)
    for i in anti_clockwise_box_idxs {
        bdr.boxes[i].swap(1, 3);
    }
    Ok(())
}

/// Assuming that the boxes all specify their vertices in clockwise order, generate a list of
/// [`Edge`]s to connect the boxes/vertices into a full graph.
fn generate_edges(bdr: &Builder) -> Result<(EdgeVec<Edge>, BoxEdgeMap), BuildError> {
    let mut box_side_map: BoxEdgeMap = HashMap::new();
    let mut edges = EdgeVec::<Edge>::new();
    // This map tracks pair of (bottom, top) vertices and the edge that joins them (if found so
    // far).  All `Edge`s are stored here, and reversing the vertex pair will test whether or
    // not (the reverse of) the given edge has already been found, in which case we have found
    // the 2nd box of that `Edge`.
    let mut edge_map = HashMap::<(VertIdx, VertIdx), EdgeIdx>::new();
    for (box_idx, vert_idxs) in bdr.boxes.indexed_iter() {
        // These indices will make the current box on the right of
        // `edges[(vert_idx_bottom, vert_idx_top)]`.
        for (side_idx, (&vert_idx_bottom, &vert_idx_top)) in
            vert_idxs.iter().circular_tuple_windows().enumerate()
        {
            // Which side of the box this edge is attaching to
            let side = match side_idx {
                0 => Side::Left,
                1 => Side::Top,
                2 => Side::Right,
                3 => Side::Bottom,
                _ => unreachable!("Boxes should only have 4 sides"),
            };
            // If the reverse of this edge has been seen before, then this the box on its left
            // side
            if let Some(existing_edge_idx) = edge_map.get(&(vert_idx_top, vert_idx_bottom)) {
                let existing_edge = &mut edges[*existing_edge_idx];
                // Check that there aren't two boxes attached to this side of the edge
                if let Some((existing_box_idx, _side)) = existing_edge.box_left {
                    return Err(BuildError::OverlappingFacesOnEdge {
                        box_idx_1: existing_box_idx,
                        box_idx_2: box_idx,
                    });
                }
                // Save which edge goes on this side of the box
                box_side_map.insert((box_idx, side), (*existing_edge_idx, EdgeSide::Left));
                // Add this box to the left side of the existing edge
                existing_edge.box_left = Some((box_idx, side));
                continue;
            }
            // If this edge is new, then add it to the map (asserting that it doesn't already
            // exist)
            //
            // Check that this side of the edge isn't already attached to a box
            if let Some(existing_edge_idx) = edge_map.get(&(vert_idx_bottom, vert_idx_top)) {
                let existing_edge = &edges[*existing_edge_idx];
                return Err(BuildError::OverlappingFacesOnEdge {
                    box_idx_1: existing_edge.box_idx_right,
                    box_idx_2: box_idx,
                });
            }

            // Get the index of the newly added edge
            let new_edge_idx = edges.push(Edge {
                vert_idx_bottom,
                vert_idx_top,
                box_left: None,
                box_idx_right: box_idx,
                box_side_right: side,
            });
            // Add this new edge to the lookup tables
            box_side_map.insert((box_idx, side), (new_edge_idx, EdgeSide::Right));
            edge_map.insert((vert_idx_bottom, vert_idx_top), new_edge_idx);
        }
    }

    // Now that all boxes have been covered, then `edge_map` must be a complete set so extract
    // the `Edge`s and return
    Ok((edges, box_side_map))
}

/// Generate all the cell vertices for the new [`Shape`], along with a [`HashMap`] dictating
/// which vertices appear in which faces
fn generate_cell_vertices(
    bdr: &Builder,
    edges: &EdgeVec<Edge>,
) -> Result<(VertVec<V2>, HashMap<VertCoord, VertIdx>), BuildError> {
    // List of all the vertices of **cells** not boxes.  This starts out as the same as
    // `b.verts`, which means that the indices specified by box corners are valid for both
    // vertex lists
    let mut cell_vert_positions = bdr.verts.map(|v| v.position);
    // This maps triples of `(box_idx, box_coord_right, box_coord_up)` to the index of the
    // vertex at that position.  In this set `(i, 0, 0)` is the bottom-left vertex of box `i`,
    // and `(i, b.box_width, b.box_height)` is the top-right vertex of box `i`.
    let mut cell_vert_indices = HashMap::<(BoxIdx, usize, usize), VertIdx>::new();

    // Add the box's vertices as cell vertices
    for (box_idx, [bl_idx, tl_idx, tr_idx, br_idx]) in bdr.boxes.indexed_iter() {
        // Single-letter variable names make the `insert` statements line up nicely, which I
        // think is a net increase in readability
        let w = bdr.box_width;
        let h = bdr.box_height;
        cell_vert_indices.insert((box_idx, 0, 0), *bl_idx); // Bottom left
        cell_vert_indices.insert((box_idx, 0, h), *tl_idx); // Top left
        cell_vert_indices.insert((box_idx, w, h), *tr_idx); // Top right
        cell_vert_indices.insert((box_idx, w, 0), *br_idx); // Bottom right
    }

    // Add vertices down each edge, checking that boxes on either side of an edge agree on the
    // number of cells
    for e in edges.iter() {
        // Add internal vertices for this edges
        let vert_idxs = get_verts_in_edge(bdr, e, &mut cell_vert_positions)?;
        for (i, vert_idx) in vert_idxs.into_iter().enumerate() {
            // `vert_idxs` doesn't include the corners, so the indices within the box start at
            // 1
            let idx_along_edge = i + 1;
            // This vertex is part of the right box (which always exists)
            let (x, y) = get_sub_box_coord(bdr, idx_along_edge, e.box_side_right, true);
            cell_vert_indices.insert((e.box_idx_right, x, y), vert_idx);
            // This vertex is part of the left box (if it exists)
            if let Some((box_idx_left, box_side_left)) = e.box_left {
                let (x, y) = get_sub_box_coord(bdr, idx_along_edge, box_side_left, false);
                cell_vert_indices.insert((box_idx_left, x, y), vert_idx);
            }
        }
    }

    // Add central vertices of each box
    for (box_idx, [bl_idx, tl_idx, tr_idx, br_idx]) in bdr.boxes.indexed_iter() {
        // Unpack the locations of the corner vertices
        let bl_pos = bdr.verts[*bl_idx].position;
        let tl_pos = bdr.verts[*tl_idx].position;
        let tr_pos = bdr.verts[*tr_idx].position;
        let br_pos = bdr.verts[*br_idx].position;
        // For each internal vertex coordinate
        for x in 1..=bdr.box_width - 1 {
            for y in 1..=bdr.box_height - 1 {
                // Compute the position of the new vertex by lerping in the y and x directions.
                // The results are the same if we were to lerp on the top/bottom edges first.
                let left_pos_ = V2::lerp(bl_pos, tl_pos, y as f32 / bdr.box_height as f32);
                let right_pos = V2::lerp(br_pos, tr_pos, y as f32 / bdr.box_height as f32);
                let pos = V2::lerp(left_pos_, right_pos, x as f32 / bdr.box_width as f32);
                // Add the new vertex to the list ...
                let new_vert_idx = cell_vert_positions.push(pos);
                // ... and add its index to the mapping
                cell_vert_indices.insert((box_idx, x, y), new_vert_idx);
            }
        }
    }

    // Convert `cell_vert_indices` to use `VertCoord` structs
    let cell_vert_indices = cell_vert_indices
        .into_iter()
        .map(|((box_idx, x, y), vert_idx)| (VertCoord { box_idx, x, y }, vert_idx))
        .collect();

    Ok((cell_vert_positions, cell_vert_indices))
}

/// Generates all the cells in each box, linking them to the corresponding vertices
fn generate_cells(
    bdr: &Builder,
    vert_idx_by_coord: &HashMap<VertCoord, VertIdx>,
) -> (CellVec<Vec<VertIdx>>, HashMap<CellCoord, CellIdx>) {
    let mut cells = CellVec::<Vec<VertIdx>>::new();
    let mut cell_vert_idxs = HashMap::<CellCoord, CellIdx>::new();
    for (box_idx, _box) in bdr.boxes.indexed_iter() {
        for x in 0..bdr.box_width {
            for y in 0..bdr.box_height {
                // Compute the coord and vertex indices
                let cell_coord = CellCoord { box_idx, x, y };
                let vert_idxs = [(0, 0), (0, 1), (1, 1), (1, 0)]
                    .iter()
                    .map(|(dx, dy)| {
                        *vert_idx_by_coord
                            .get(&VertCoord {
                                box_idx,
                                x: x + dx,
                                y: y + dy,
                            })
                            .unwrap()
                    })
                    .collect_vec();
                // Add the cell to the lookup tables
                let cell_idx = cells.push(vert_idxs);
                cell_vert_idxs.insert(cell_coord, cell_idx);
            }
        }
    }

    (cells, cell_vert_idxs)
}

/// Generates which groups of cells must not contain the same number
fn generate_groups<'e>(
    bdr: &Builder,
    edges: &'e EdgeVec<Edge>,
    cell_idx_by_coord: &HashMap<CellCoord, CellIdx>,
    box_edge_map: &BoxEdgeMap,
) -> Vec<Vec<CellIdx>> {
    let mut groups = Vec::<Vec<CellIdx>>::new();

    // Group cells into boxes (no graph traversal required)
    for (box_idx, _box) in bdr.boxes.indexed_iter() {
        let box_group = (0..bdr.box_width)
            .cartesian_product(0..bdr.box_height)
            .map(|(x, y)| *cell_idx_by_coord.get(&CellCoord { box_idx, x, y }).unwrap())
            .collect_vec();
        groups.push(box_group);
    }

    /* Group cells into rows/columns (known generally as 'lanes') */

    // In order to make sure that every lane is covered, we need to traverse every possible quad
    // path through the boxes (i.e. paths where we always exit from the opposite edge to the one we
    // came).  If all possible quad paths are traversed, then every edge must be traversed exactly
    // once.  Therefore, we keep a HashSet of untraversed edges.  Until this set is empty, we
    // repeatedly pop an edge and then explore the quad-path which contains that edge (removing the
    // edges from the set as we find them).  Edges can't exist in multiple quad paths, so this
    // algorithm guarantees that every possible quad path must be covered.
    let mut untraversed_edges: HashSet<&'e Edge> = edges.iter().collect();

    while !untraversed_edges.is_empty() {
        // Pick external edges until only internal edges are left (almost all paths will start and
        // end at an external edge, but cyclic paths are possible and have no external edges).
        let mut next_external_edge: Option<&'e Edge> = None;
        for e in &untraversed_edges {
            if e.box_left.is_none() {
                next_external_edge = Some(*e);
            }
        }
        // Traverse the box graph and get which boxes are covered by the path
        let box_path = match next_external_edge {
            Some(e) => traverse_non_cyclic_path(edges, e, &mut untraversed_edges, box_edge_map),
            None => todo!("Cyclic paths aren't implemented yet"),
        };

        // Split this box path into individual lanes of cells, and add those lanes as groups
        let lanes = get_lanes_down_path(bdr, &box_path, cell_idx_by_coord);
        groups.extend(lanes);
    }

    groups
}

/// Traverse a path starting by entering a box at a given [`Edge`].  This returns a list of
/// `(box_idx, side_entered)`.
fn traverse_non_cyclic_path(
    edges: &EdgeVec<Edge>,
    start_edge: &Edge,
    untraversed_edges: &mut HashSet<&Edge>,
    box_edge_map: &BoxEdgeMap,
) -> Vec<(BoxIdx, Side)> {
    // Assert that the start edge only has one side (otherwise the path is non-cyclic)
    assert!(start_edge.box_left.is_none());

    let mut box_path = Vec::<(BoxIdx, Side)>::new();

    let mut cur_edge = start_edge; // Which edge we're on
    let mut cur_box_idx = start_edge.box_idx_right; // The box we're entering
    let mut cur_box_side = start_edge.box_side_right; // The side of the box we're entering

    // This path must be non-cyclic, so it must stop when we reach an edge with only one box (which
    // must be on the right side of that edge).
    loop {
        // By stepping over the box, we have traversed the current edge (checking that it hasn't
        // been traversed before)
        assert!(untraversed_edges.remove(cur_edge));
        box_path.push((cur_box_idx, cur_box_side));
        // Step to the opposite side of the current box
        let exit_side = cur_box_side.opposite();
        let (new_edge_idx, new_edge_approach_side) =
            *box_edge_map.get(&(cur_box_idx, exit_side)).unwrap();
        // Step to the opposite side of the new edge
        let new_edge = &edges[new_edge_idx];
        let new_edge_exit_side = new_edge_approach_side.opposite();
        // Determine which box we've stepped into
        let (new_box_idx, new_box_side) = match new_edge_exit_side {
            // Edges always have a box on their right side
            EdgeSide::Right => (new_edge.box_idx_right, new_edge.box_side_right),
            EdgeSide::Left => {
                if let Some(box_left) = new_edge.box_left {
                    box_left
                } else {
                    // Make sure to mark this last edge as traversed, since otherwise this entire
                    // path would be traversed again (but going in the opposite direction).
                    assert!(untraversed_edges.remove(new_edge));
                    // If the edge doesn't have a left side then we've finished the path, so should
                    // break the loop
                    break;
                }
            }
        };
        cur_edge = new_edge;
        cur_box_idx = new_box_idx;
        cur_box_side = new_box_side;
    }

    box_path
}

//////////////////////
// HELPER FUNCTIONS //
//////////////////////

type BoxEdgeMap = HashMap<(BoxIdx, Side), (EdgeIdx, EdgeSide)>;

/// Convert a path through boxes into a set of lanes of cells which form the groups down that path
fn get_lanes_down_path(
    bdr: &Builder,
    box_path: &[(BoxIdx, Side)],
    cell_idx_by_coord: &HashMap<CellCoord, CellIdx>,
) -> Vec<Vec<CellIdx>> {
    let num_lanes = box_size_in_direction(bdr, box_path[0].1.direction());
    let lane_depth = box_size_in_direction(bdr, !box_path[0].1.direction());
    // The lanes are numbered going clockwise around the perimeter of each box
    let mut lanes = vec![Vec::<CellIdx>::new(); num_lanes];

    for &(box_idx, entry_side) in box_path {
        // Sanity check that the edges being stepped over don't change the number or depth of
        // the lanes.  This should have been caught in `generate_edges` (and an error
        // returned), but it doesn't hurt to double check.
        assert_eq!(
            box_size_in_direction(bdr, entry_side.direction()),
            num_lanes
        );
        assert_eq!(
            box_size_in_direction(bdr, !entry_side.direction()),
            lane_depth
        );
        // We know that the box size is valid, so add this box's cells to the lanes
        for lane_idx in 0..num_lanes {
            for idx_down_lane in 0..lane_depth {
                // Rotate the coordinates according to the side we enter from:
                //
                //   TOP        LEFT        BOTTOM      RIGHT
                // +-----+     +-----+     +-----+     +-----+
                // |0 1 2|     |2 2 2|     |2 1 0|     |0 0 0|
                // |0 1 2|     |1 1 1|     |2 1 0|     |1 1 1|
                // |0 1 2|     |0 0 0|     |2 1 0|     |2 2 2|
                // +-----+     +-----+     +-----+     +-----+
                // (coordinates go from bottom-left corner)
                let (x, y) = match entry_side {
                    Side::Top => (lane_idx, lane_depth - 1 - idx_down_lane),
                    Side::Left => (idx_down_lane, lane_idx),
                    Side::Bottom => (num_lanes - 1 - lane_idx, idx_down_lane),
                    Side::Right => (lane_depth - 1 - idx_down_lane, num_lanes - 1 - lane_idx),
                };
                // Extend the corresponding lane
                lanes[lane_idx].push(*cell_idx_by_coord.get(&CellCoord { box_idx, x, y }).unwrap());
            }
        }
    }

    lanes
}

/// Computes the internal vertices of each [`Edge`], adding them to `vert_vec` and returning
/// the indices.  These vertices go from the 'bottom' to the 'top' of this edge
fn get_verts_in_edge(
    bdr: &Builder,
    edge: &Edge,
    vert_vec: &mut VertVec<V2>,
) -> Result<Vec<VertIdx>, BuildError> {
    // Decide how many cells this edge will be split into (specified by the right-hand box)
    let num_cells = get_num_cells_in_side(bdr, edge.box_side_right);
    // Check that the left-hand box agrees with this (if it exists)
    if let Some((box_idx_left, box_side_left)) = edge.box_left {
        let num_cells_left = get_num_cells_in_side(bdr, box_side_left);
        if num_cells_left != num_cells {
            return Err(BuildError::InconsistentEdge {
                box_idx_left,
                num_cells_left,
                box_idx_right: edge.box_idx_right,
                num_cells_right: num_cells,
            });
        }
    }

    // Split this cell into `num_cells`, computing the locations of the 'internal' vertices -
    // i.e. those which aren't on the corner of a box
    let vert_indices = (1..=num_cells - 1)
        .map(|cells_from_bot| {
            let position = V2::lerp(
                bdr.verts[edge.vert_idx_bottom].position,
                bdr.verts[edge.vert_idx_top].position,
                cells_from_bot as f32 / num_cells as f32,
            );
            vert_vec.push(position)
        })
        .collect_vec();
    Ok(vert_indices)
}

/// Returns the number of cells in a given [`Side`] of any box
fn get_num_cells_in_side(bdr: &Builder, side: Side) -> usize {
    match side.direction() {
        Direction::Vertical => bdr.box_height,
        Direction::Horizontal => bdr.box_width,
    }
}

/// Gets the sub-box (integer) coordinate of a vertex, given its index along a given [`Side`]
/// of a box.  If `reverse` is true, the indices will count anti-clockwise.
fn get_sub_box_coord(
    bdr: &Builder,
    idx_along_edge: usize,
    side: Side,
    clockwise_idx: bool,
) -> (usize, usize) {
    // Reverse the index if we're going anticlockwise
    let idx_along_edge = if clockwise_idx {
        idx_along_edge // The index is already in the clockwise direction
    } else if side.direction() == Direction::Vertical {
        bdr.box_height - idx_along_edge
    } else {
        bdr.box_width - idx_along_edge
    };
    match side {
        Side::Left => (0, idx_along_edge),
        Side::Top => (idx_along_edge, bdr.box_height),
        Side::Right => (bdr.box_width, bdr.box_height - idx_along_edge),
        Side::Bottom => (bdr.box_width - idx_along_edge, 0),
    }
}

/// Gets the number of cells along any [`Edge`] of a box in a given [`Direction`].
fn box_size_in_direction(bdr: &Builder, direction: Direction) -> usize {
    match direction {
        Direction::Horizontal => bdr.box_width,
        Direction::Vertical => bdr.box_height,
    }
}

/* UTILITY TYPES FOR SHAPE GENERATION */

/// The different states of the corner of a box
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
enum CornerType {
    /// The 3 vertices round this corner turn clockwise
    Clockwise,
    /// The 3 vertices lie on a straight line (and therefore don't turn either direction)
    Linear,
    /// The 3 vertices round this corner turn anti-clockwise
    AntiClockwise,
}

impl CornerType {
    fn classify(angle: f32) -> CornerType {
        if angle.abs() < 0.001 || angle.abs() > 2.0 * PI - 0.001 {
            // If the angle is close to 0 or 180 degrees, we call the corner 'linear'
            CornerType::Linear
        } else if angle > 0.0 {
            // If the corner increases the angle, then it turns clockwise
            CornerType::Clockwise
        } else {
            // If the corner decreases the angle, then it turns anti-clockwise
            CornerType::AntiClockwise
        }
    }
}

/// The coordinates of a (cell) vertex within a box of the built shape.  Note that the same vertex
/// can be referred to by many `VertCoordinate` if it is shared between many faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct VertCoord {
    box_idx: BoxIdx,
    x: usize,
    y: usize,
}

/// The coordinates of a cell within a box of the [`Shape`] being built.  This has an identical
/// shape to [`VertCoord`], but the extra type safety of the two types is useful.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct CellCoord {
    box_idx: BoxIdx,
    x: usize,
    y: usize,
}

/// An edge between two **boxes** (unlike [`kuudos::shape::Edge`], which sits between two
/// **cells**.  If an [`Edge`] only has one adjacent box, it will be put on the 'right' side of
/// the [`Edge`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Edge {
    vert_idx_top: VertIdx,
    vert_idx_bottom: VertIdx,

    box_left: Option<(BoxIdx, Side)>,
    box_idx_right: BoxIdx,
    box_side_right: Side,
}

/// The possible directions through a box
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
enum Direction {
    Horizontal,
    Vertical,
}

impl Not for Direction {
    type Output = Direction;

    fn not(self) -> Self::Output {
        match self {
            Direction::Vertical => Direction::Horizontal,
            Direction::Horizontal => Direction::Vertical,
        }
    }
}

/// The possible each edges of each box
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
enum Side {
    Left,
    Top,
    Right,
    Bottom,
}

impl Side {
    /// The direction of the `Edge`'s **line**
    fn direction(self) -> Direction {
        match self {
            Side::Left | Side::Right => Direction::Vertical,
            Side::Top | Side::Bottom => Direction::Horizontal,
        }
    }

    fn opposite(self) -> Side {
        match self {
            Side::Left => Side::Right,
            Side::Top => Side::Bottom,
            Side::Right => Side::Left,
            Side::Bottom => Side::Top,
        }
    }
}

/// The two sides of an edge (assuming that the edge is pointing 'up')
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
enum EdgeSide {
    Right,
    Left,
}

impl EdgeSide {
    fn opposite(self) -> EdgeSide {
        match self {
            EdgeSide::Right => EdgeSide::Left,
            EdgeSide::Left => EdgeSide::Right,
        }
    }
}

/// The angle (going **clockwise**) between the directions of two vectors.  This is always in the
/// region `-PI < x <= PI`.
fn angle_between(v1: V2, v2: V2) -> f32 {
    // The clockwise angle from the x-axis to the direction of `v1` and `v2`
    let angle1 = f32::atan2(v1.y, v1.x);
    let angle2 = f32::atan2(v2.y, v2.x);

    let mut angle_diff = angle2 - angle1;
    // Rotate angle_diff into the range `-2 * PI <= x < 2 * PI`
    angle_diff %= PI * 2.0;
    // If the absolute difference is bigger than PI (i.e. >180 degrees) then rotate it back into
    // the range `-PI < x <= PI`
    if angle_diff > PI {
        angle_diff -= PI * 2.0;
    }
    if angle_diff <= -PI {
        angle_diff += PI * 2.0;
    }

    angle_diff
}
