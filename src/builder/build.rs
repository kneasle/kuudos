//! Code for converting a `Builder` into a [`Shape`]/[`Symmetry`] pair

use std::collections::{HashMap, HashSet};

use itertools::Itertools;

use crate::{
    image::{Elem, StrokeStyle, Style},
    indexed_vec::{BoxIdx, CellIdx, CellVec, EdgeIdx, EdgeVec, IdxType, LinkIdx, VertIdx, VertVec},
    shape::Group,
    Shape, Symmetry, V2,
};

use super::{Builder, Direction, EdgeLinkStyle, LinkSide, RotateDirection, Side};

/// Converts a `Builder` into a [`Shape`] and the associated [`Symmetry`]
pub fn build(mut bdr: Builder) -> Result<(Shape, Symmetry), BuildError> {
    // Before we start building the `Shape`, we normalise each box so that the vertices are
    // always specified in clockwise order.
    normalise_box_direction(&mut bdr);

    // Generate edges between boxes
    let (edges, box_edge_map, link_edge_map) = generate_edges(&bdr)?;
    // Generate cell vertices
    let (verts, vert_idx_by_coord) = generate_cell_vertices(&bdr, &edges)?;
    // Generate cells
    let (cells, cell_idx_by_coord) = generate_cells(&bdr, &vert_idx_by_coord);
    // Generate extra elements to render the edge links
    let extra_elements = generate_extra_elements(&bdr, &link_edge_map, &edges);
    // Generate groups
    let groups = generate_groups(
        &bdr,
        &edges,
        &cell_idx_by_coord,
        &box_edge_map,
        &link_edge_map,
    );
    // Package into a `Shape`
    let num_cells = cells.len();
    let shape = Shape {
        num_symbols: bdr.box_width * bdr.box_height,
        groups,
        cells,
        verts,
        extra_elements,
    };
    // Generate the symmetry (currently just assume - incorrectly - that there's no symmetry)
    let symmetry = generate_symmetry(&bdr, num_cells, &cell_idx_by_coord);
    Ok((shape, symmetry))
}

/// The possible ways that converting a [`Builder`] into a [`Shape`] (using [`Builder::build`])
/// could fail.
#[derive(Debug, Clone)]
pub enum BuildError {
    /// An edge has two or more connections _to the same side_ of it.  Note that edges often have
    /// two attachments, but they are attached in opposite directions (because vertices are always
    /// numbered in a consistent direction):
    /// ```text
    /// 0 ------ (1 2) ------ 3
    /// |   -->    |    -->   |
    /// | ^      | | ^      | |
    /// | |      V | |      V |
    /// |   <--    |    <--   |
    /// 3 ------ (2 1) ------ 0
    /// ```
    /// Overlapping boxes would look like this:
    /// ```text
    ///     0 - 1
    ///    / \ / \
    ///   /   X   \
    ///  /   / \   \
    /// 3 - 2   3 - 2
    /// ```
    OverlappingEdgeConnections {
        connection_1: EdgeConnection,
        connection_2: EdgeConnection,
    },
    /// Two boxes on either side of an edge disagree on the number of cells down that edge
    InconsistentEdge {
        box_idx_left: BoxIdx,
        num_cells_left: usize,

        box_idx_right: BoxIdx,
        num_cells_right: usize,
    },
    /// The boxes on either side of a link disagree on the width of the link
    InconsistentLink {
        link_idx: LinkIdx,
        num_cells_top: usize,
        num_cells_bottom: usize,
    },
}

/////////////////////////
// GENERATOR FUNCTIONS //
/////////////////////////

/// Reorder the vertices in each box so that the vertices are always specified in clockwise
/// order.
fn normalise_box_direction(bdr: &mut Builder) {
    for (_idx, source_box) in bdr.source_boxes_mut() {
        if source_box.rotate_direction == RotateDirection::AntiClockwise {
            source_box.swap_direction();
        }
    }
}

/// Assuming that the boxes all specify their vertices in clockwise order, generate a list of
/// [`Edge`]s to connect the boxes/vertices into a full graph.  This checks for the case where two
/// items are attached to the same side of an edge (see
/// [`BuildError::OverlappingEdgeConnections`]).  Additionally, this verifies that all the edges
/// and links are well defined - i.e. if an edge sits between two items (boxes or links), then
/// both sides agree on how many cells the edge should cover.
fn generate_edges(bdr: &Builder) -> Result<(EdgeVec<Edge>, BoxEdgeMap, LinkEdgeMap), BuildError> {
    let mut box_edge_map: BoxEdgeMap = HashMap::new(); // Which edge lies on each side of each box
    let mut link_edge_map: LinkEdgeMap = HashMap::new(); // Which edge lies on each side of each link
    let mut edges = EdgeVec::<Edge>::new();

    // This map tracks pair of (bottom, top) vertices and the edge that joins them (if found so
    // far).  All `Edge`s are stored here, and reversing the vertex pair will test whether or
    // not (the reverse of) the given edge has already been found, in which case we have found
    // the 2nd box of that `Edge`.
    let mut edge_map = HashMap::<(VertIdx, VertIdx), EdgeIdx>::new();

    for (box_idx, source_box) in bdr.source_boxes() {
        assert_eq!(source_box.rotate_direction, RotateDirection::Clockwise);

        // These indices will make the current box on the right of
        // `edges[(vert_idx_bottom, vert_idx_top)]`.
        for (side_idx, (&vert_idx_bottom, &vert_idx_top)) in source_box
            .vert_idxs
            .iter()
            .circular_tuple_windows()
            .enumerate()
        {
            // Which side of the box this edge is attaching to
            let side = Side::from_index(side_idx);
            let connection = EdgeConnection::Box_(box_idx, side);
            // Number of cells that go down this side of the box
            let length = bdr.box_size_in_direction(side.direction());

            // If the reverse of this edge has been seen before, then this the box on its left
            // side
            if let Some(existing_edge_idx) = edge_map.get(&(vert_idx_top, vert_idx_bottom)) {
                let existing_edge = &mut edges[*existing_edge_idx];
                // Error if the other box doesn't agree on the length of this edge
                if existing_edge.length != length {
                    // At this point, edges have only been connected to boxes, so it's OK to
                    // extract the box index.
                    let existing_box_idx = match existing_edge.connection_right {
                        EdgeConnection::Box_(idx, _side) => idx,
                        // No links have been connected yet
                        EdgeConnection::Link(_, _) => unreachable!(),
                    };
                    return Err(BuildError::InconsistentEdge {
                        // The current box would be on the 'left' of the edge
                        box_idx_left: box_idx,
                        num_cells_left: length,
                        // The existing box is on the right
                        box_idx_right: existing_box_idx,
                        num_cells_right: existing_edge.length,
                    });
                }
                // Check that there aren't two boxes attached to this side of the edge
                if let Some(existing_connection) = existing_edge.connection_left {
                    return Err(BuildError::OverlappingEdgeConnections {
                        connection_1: existing_connection,
                        connection_2: connection,
                    });
                }
                // Save which edge goes on this side of the box
                box_edge_map.insert((box_idx, side), (*existing_edge_idx, EdgeSide::Left));
                // Add this box to the left side of the existing edge
                existing_edge.connection_left = Some(connection);
                continue;
            }
            // If this edge is new, then add it to the map (asserting that it doesn't already
            // exist)
            //
            // Check that this side of the edge isn't already attached to a box
            if let Some(existing_edge_idx) = edge_map.get(&(vert_idx_bottom, vert_idx_top)) {
                let existing_edge = &edges[*existing_edge_idx];
                return Err(BuildError::OverlappingEdgeConnections {
                    connection_1: connection,
                    connection_2: existing_edge.connection_right,
                });
            }

            // Create a new edge and store it
            let new_edge_idx = edges.push(Edge {
                vert_idx_bottom,
                vert_idx_top,
                connection_left: None, // This will get filled later
                connection_right: EdgeConnection::Box_(box_idx, side),
                // Used to check that the other side of the edge agrees with the length
                length,
            });
            // Add this new edge to the lookup tables
            box_edge_map.insert((box_idx, side), (new_edge_idx, EdgeSide::Right));
            edge_map.insert((vert_idx_bottom, vert_idx_top), new_edge_idx);
        }
    }

    // Once all the boxes have created their edges, we create the edges which are attached to the
    // edge links.  This is more complicated than it was for boxes, since edge links don't
    // directly specify their width.  On its own this would be fine, but additionally edge links
    // can be chained together (i.e. creating an edge bordered by two edge links and no boxes).
    //
    // Also, the edge links could (in theory) be created in any order which could mean that we
    // create an edge in the middle of a chain of links without knowing its length.  This would
    // break the contract of `Edge` (which must have a specific length), so in these cases we
    // create a `LengthlessEdge`.  Once we know the length of the `LengthlessEdge`, it is converted
    // to a full `Edge` and stored.
    for (link_idx, link) in bdr.edge_links.indexed_iter() {
        // These are labelled (bottom_vert, top_vert, side of link), so that this link is always on
        // the right-hand side of the edge
        let sides = [
            (
                link.vert_idx_top_left,
                link.vert_idx_top_right,
                LinkSide::Top,
            ),
            (
                link.vert_idx_bottom_right,
                link.vert_idx_bottom_left,
                LinkSide::Bottom,
            ),
        ];

        let mut link_width: Option<usize> = None;
        for (vert_idx_bottom, vert_idx_top, link_side) in sides {
            let edge_connection = EdgeConnection::Link(link_idx, link_side);

            // Connect to the left side of an edge if one already exists on this side
            if let Some(&existing_edge_idx) = edge_map.get(&(vert_idx_top, vert_idx_bottom)) {
                let existing_edge = &mut edges[existing_edge_idx];
                // Use this edge's length to determine the number of lanes in this link
                match link_width {
                    Some(len) => {
                        if len != existing_edge.length {
                            // If this link already has a width that's different to the one
                            // requested by this edge, then this link joins two edges with
                            // different lengths, which is an error
                            return Err(BuildError::InconsistentLink {
                                link_idx,
                                num_cells_top: len, // Top is always checked before bottom
                                num_cells_bottom: existing_edge.length,
                            });
                        }
                    }
                    None => link_width = Some(existing_edge.length),
                }
                // Add this link as the left-hand connection of the edge
                existing_edge.connection_left = Some(edge_connection);
                link_edge_map.insert((link_idx, link_side), (existing_edge_idx, EdgeSide::Left));
            } else {
                unimplemented!("Chaining edge links isn't implemented yet");
            }
        }
    }

    // Now that all boxes and links have created their edges, `edge_map` must be complete so we
    // can extract the `Edge`s and return
    Ok((edges, box_edge_map, link_edge_map))
}

/// Generate all the cell vertices for the new [`Shape`], along with a [`HashMap`] dictating
/// which vertices appear in which faces.
fn generate_cell_vertices(
    bdr: &Builder,
    edges: &EdgeVec<Edge>,
) -> Result<(VertVec<V2>, HashMap<VertCoord, VertIdx>), BuildError> {
    // List of all the vertices of **cells** not boxes.  This starts out as the same as
    // `b.verts`, which means that the indices specified by box corners are valid for both
    // vertex lists
    let mut cell_vert_positions = bdr.verts.clone();
    // This maps triples of `(box_idx, box_coord_right, box_coord_up)` to the index of the
    // vertex at that position.  In this set `(i, 0, 0)` is the bottom-left vertex of box `i`,
    // and `(i, b.box_width, b.box_height)` is the top-right vertex of box `i`.
    let mut cell_vert_indices = HashMap::<(BoxIdx, usize, usize), VertIdx>::new();

    // Add the box's vertices as cell vertices
    for (box_idx, source_box) in bdr.source_boxes() {
        // Single-letter variable/direction names make the `insert` statements line up nicely,
        // which I think is a net increase in readability
        let [bl_idx, tl_idx, tr_idx, br_idx] = source_box.vert_idxs;
        let w = bdr.box_width;
        let h = bdr.box_height;
        cell_vert_indices.insert((box_idx, 0, 0), bl_idx); // Bottom left
        cell_vert_indices.insert((box_idx, 0, h), tl_idx); // Top left
        cell_vert_indices.insert((box_idx, w, h), tr_idx); // Top right
        cell_vert_indices.insert((box_idx, w, 0), br_idx); // Bottom right
    }

    // Add vertices down each edge, checking that boxes on either side of an edge agree on the
    // number of cells
    for edge in edges.iter() {
        // Compute the positions of internal vertices for this edge, and label them with their
        // index from the bottom of this edge
        let vert_indices = (1..=edge.length - 1).map(|num_cells_from_bottom| {
            let position = V2::lerp(
                bdr.verts[edge.vert_idx_bottom],
                bdr.verts[edge.vert_idx_top],
                num_cells_from_bottom as f32 / edge.length as f32,
            );
            (num_cells_from_bottom, cell_vert_positions.push(position))
        });

        for (idx_along_edge, vert_idx) in vert_indices {
            // Make this vertex part of the right box
            if let EdgeConnection::Box_(box_idx, side) = edge.connection_right {
                let (x, y) = get_sub_box_coord(bdr, idx_along_edge, side, true);
                cell_vert_indices.insert((box_idx, x, y), vert_idx);
            }
            // Make this vertex part of the left box
            if let Some(EdgeConnection::Box_(box_idx, side)) = edge.connection_left {
                let (x, y) = get_sub_box_coord(bdr, idx_along_edge, side, false);
                cell_vert_indices.insert((box_idx, x, y), vert_idx);
            }
        }
    }

    // Add central vertices of each box
    for (box_idx, source_box) in bdr.source_boxes() {
        let [bl_idx, tl_idx, tr_idx, br_idx] = source_box.vert_idxs;
        // Unpack the locations of the corner vertices
        let bl_pos = bdr.verts[bl_idx];
        let tl_pos = bdr.verts[tl_idx];
        let tr_pos = bdr.verts[tr_idx];
        let br_pos = bdr.verts[br_idx];
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
    for (box_idx, _box) in bdr.source_boxes() {
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

/// Generate which extra elements have to be drawn to display the edge links
fn generate_extra_elements(
    bdr: &Builder,
    link_edge_map: &LinkEdgeMap,
    edges: &EdgeVec<Edge>,
) -> Vec<Elem> {
    let mut extra_elems = Vec::new();

    // Closure which returns the extra links required to draw a single `EdgeLink`
    for (link_idx, link) in bdr.edge_links.indexed_iter() {
        // Cheeky closure to get the length of one side of the link
        let get_edge_len = |side: LinkSide| {
            let (edge_idx, _) = link_edge_map[&(link_idx, side)];
            edges[edge_idx].length
        };
        // Get the edges on top & bottom of this link
        let num_lanes = get_edge_len(LinkSide::Top);
        // Sanity check that both edges agree on width - this should be converted to an error in
        // `generate_edges`
        assert_eq!(get_edge_len(LinkSide::Bottom), num_lanes);

        // Get the vertex positions
        let vert_pos_tl = bdr.verts[link.vert_idx_top_left];
        let vert_pos_tr = bdr.verts[link.vert_idx_top_right];
        let vert_pos_bl = bdr.verts[link.vert_idx_bottom_left];
        let vert_pos_br = bdr.verts[link.vert_idx_bottom_right];
        // Compute edge normals (i.e. normals facing into the link)
        let top_normal = (vert_pos_tr - vert_pos_tl).normal().normalise();
        let _bottom_normal = (vert_pos_bl - vert_pos_br).normal().normalise();
        // Generate elements
        for lane_idx in 0..num_lanes {
            let lerp_factor = (lane_idx as f32 + 0.5) / num_lanes as f32;
            let top_pos = V2::lerp(vert_pos_tl, vert_pos_tr, lerp_factor);
            let bottom_pos = V2::lerp(vert_pos_bl, vert_pos_br, lerp_factor);
            // Generate the shape of this element
            let shape = match link.style {
                EdgeLinkStyle::Arc => Elem::arc_passing_through(
                    top_pos,
                    top_normal,
                    bottom_pos,
                    Style::JustStroke(StrokeStyle::BoxBorder),
                )
                .unwrap_or(Elem::LineSegment(
                    top_pos,
                    bottom_pos,
                    StrokeStyle::BoxBorder,
                )),
                EdgeLinkStyle::Linear => {
                    Elem::LineSegment(top_pos, bottom_pos, StrokeStyle::BoxBorder)
                }
                EdgeLinkStyle::Hidden => continue, // Don't add hidden links
            };
            extra_elems.push(shape);
        }
    }

    extra_elems
}

/// Generates which groups of cells must not contain the same number
fn generate_groups<'e>(
    bdr: &Builder,
    edges: &'e EdgeVec<Edge>,
    cell_idx_by_coord: &HashMap<CellCoord, CellIdx>,
    box_edge_map: &BoxEdgeMap,
    link_edge_map: &LinkEdgeMap,
) -> Vec<Group> {
    let mut groups = Vec::<Group>::new();

    // Group cells into boxes (no graph traversal required)
    for (box_idx, _box) in bdr.source_boxes() {
        let cells = (0..bdr.box_width)
            .cartesian_product(0..bdr.box_height)
            .map(|(x, y)| cell_idx_by_coord[&CellCoord { box_idx, x, y }])
            .collect_vec();
        groups.push(Group::box_group(cells));
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

    loop {
        // `traverse_path_forward` will only traverse paths starting from one side of an edge.
        // Therefore, calling `traverse_path_forward` on an edge that is in the middle of a
        // non-cyclic path would miss out on part of the path (`traverse_path_forward` actually
        // checks this condition, and will panic on invalid inputs).
        //
        // We guarantee this giving external edges priority when choosing them.  This way, all the
        // external edges will have their paths explored before any internal edge is expanded.  If
        // all the external edges have been traversed, then the only remaining paths must be
        // cyclic.
        let next_external_edge = untraversed_edges
            .iter()
            .find(|e| e.connection_left.is_none());
        let next_edge = match next_external_edge {
            Some(external_edge) => external_edge,
            None => match untraversed_edges.iter().next() {
                Some(internal_edge_in_cyclic_path) => internal_edge_in_cyclic_path,
                // If there are no more edges to traverse, then all paths have been covered
                None => break,
            },
        };

        let box_path = traverse_path_forward(
            edges,
            next_edge,
            &mut untraversed_edges,
            box_edge_map,
            link_edge_map,
        );

        // Split this box path into individual lanes of cells, and add those lanes as groups
        let lanes = get_lanes_down_path(bdr, &box_path, cell_idx_by_coord);
        groups.extend(lanes.into_iter().map(Group::non_box_group));
    }

    groups
}

/// Generate the [`Symmetry`] of the [`Shape`] represented by the [`Builder`].
fn generate_symmetry(
    bdr: &Builder,
    num_cells: usize,
    cell_idx_by_coord: &HashMap<CellCoord, CellIdx>,
) -> Symmetry {
    /* In order to deduce the symmetry, we must first assign each cell into a (not necessarily
     * unique) equivalence class.  We can (and will) use `Builder`'s equivalence classes of faces,
     * but this is complicated by the fact that multiple `Box_`es can exist on top of one another -
     * thus assigning some cells to multiple equivalence classes.  It doesn't make sense to assign
     * many classes to the same cell, so instead the correct thing to do here is to combine these
     * multiple equivalence classes into one.  We implement this by making high-indexed classes
     * dominate lower-indexed ones, and having all classes dominate `None` (this is how `Ord` is
     * implemented for `Option<T>`). */
    let mut equiv_class_by_cell: Vec<Option<usize>> = vec![None; num_cells];
    let mut next_equiv_class_idx = 0;
    for box_idxs in bdr.box_equiv_classes.iter() {
        // Iterate over the cell coords **first**, because cells at the same location in different
        // boxes should be put in the same equivalence class.
        for x in 0..bdr.box_width {
            for y in 0..bdr.box_height {
                // Get the equivalence class index for all the cells in this location
                let equiv_class_idx = next_equiv_class_idx;
                next_equiv_class_idx += 1;
                // For each box in this equiv class of boxes...
                for &box_idx in box_idxs {
                    // Get the source box and the transform from the this box
                    let (source_idx, _source_box, transform) = bdr.get_source_box(box_idx).unwrap();
                    let (transformed_x, transformed_y) =
                        transform.transform_coord(bdr.box_width, bdr.box_height, x, y);
                    // Determine which cell this coordinate belongs to
                    let cell_coord = CellCoord {
                        box_idx: source_idx,
                        x: transformed_x,
                        y: transformed_y,
                    };
                    let cell_idx = cell_idx_by_coord.get(&cell_coord).unwrap();
                    // Note that this cell belongs to this equiv_class
                    let equiv_class = &mut equiv_class_by_cell[cell_idx.to_idx()];
                    *equiv_class = (*equiv_class).max(Some(equiv_class_idx));
                }
            }
        }
    }
    // Ever cell must have come from a box, so no cell can be left without an equivalence class
    // once all boxes have been iterated over
    assert!(equiv_class_by_cell.iter().all(Option::is_some));
    // Build a `Symmetry` from this cell assignment
    Symmetry::from_cell_assignments(equiv_class_by_cell.iter())
}

//////////////////////
// HELPER FUNCTIONS //
//////////////////////

type BoxEdgeMap = HashMap<(BoxIdx, Side), (EdgeIdx, EdgeSide)>;
type LinkEdgeMap = HashMap<(LinkIdx, LinkSide), (EdgeIdx, EdgeSide)>;

/// Traverse a path starting by entering a box at a given [`Edge`].  This returns a list of
/// `(box_idx, side_entered)`.  **NOTE**: This will not traverse backwards from `start_edge`, so
/// it's up to the caller to check that this edge is either 'external' (i.e. is only connected to
/// one box or link) or is part of a cyclic path.  This will panic if a non-conforming edge is
/// given.
fn traverse_path_forward(
    edges: &EdgeVec<Edge>,
    start_edge: &Edge,
    untraversed_edges: &mut HashSet<&Edge>,
    box_edge_map: &BoxEdgeMap,
    link_edge_map: &LinkEdgeMap,
) -> Vec<(BoxIdx, Side)> {
    let mut box_path = Vec::<(BoxIdx, Side)>::new();

    let mut cur_edge = start_edge; // Which edge we're on
    let mut cur_connection = start_edge.connection_right; // Start the path on the right of the edge

    // This path must be non-cyclic, so it must stop when we reach an edge with only one box (which
    // must be on the right side of that edge).
    loop {
        // By stepping over the box/link, we have traversed the current edge.  So mark the current
        // edge as traversed by removing it from the `untraversed_edges` set (checking that it
        // hasn't been traversed/removed before)
        assert!(untraversed_edges.remove(cur_edge));

        // Step across the current `EdgeConnection`
        let (new_edge_idx, new_edge_approach_side) = match cur_connection {
            // If we're looking at a box, then add it to the path and traverse
            EdgeConnection::Box_(cur_box_idx, cur_box_side) => {
                // Add the box to the path
                box_path.push((cur_box_idx, cur_box_side));
                // Step to the opposite side of the box
                let exit_side = cur_box_side.opposite();
                box_edge_map[&(cur_box_idx, exit_side)]
            }
            // Edge links don't form part of the path, so just traverse them
            EdgeConnection::Link(link_idx, link_side) => {
                // Step to the opposite side of the link
                let exit_side = !link_side;
                link_edge_map[&(link_idx, exit_side)]
            }
        };

        // Step across the newly approached edge
        let new_edge = &edges[new_edge_idx];
        let new_edge_exit_side = new_edge_approach_side.opposite();

        // Determine what connection we've stepped into
        let new_connection = match new_edge_exit_side {
            // Edges always have a box on their right side
            EdgeSide::Right => new_edge.connection_right,
            EdgeSide::Left => {
                match new_edge.connection_left {
                    // Step to the left connection if it exists
                    Some(conn_left) => conn_left,
                    // If there's nothing connected to the other side of the edge, then we've
                    // reached the end of the path
                    None => {
                        // Mark this last edge as traversed. If this isn't done, the entire path
                        // would be traversed again (but going in the opposite direction).
                        assert!(untraversed_edges.remove(new_edge));
                        // We now know that this path is non-cyclic, so check that the starting
                        // edge only bordered one side (otherwise we only explored part of the path)
                        assert!(start_edge.connection_left.is_none());
                        // If the edge doesn't have a left side then we've finished the path, so
                        // should break the loop
                        break;
                    }
                }
            }
        };
        cur_edge = new_edge; // Move to the new edge
        cur_connection = new_connection; // Move to the new connection

        // If we reach the start edge again, then we have fully traversed the (cyclic) path and
        // should terminate the loop
        if cur_edge == start_edge && cur_connection == start_edge.connection_right {
            break;
        }
    }

    box_path
}

/// Convert a path through boxes into a set of lanes of cells which form the groups down that path
fn get_lanes_down_path(
    bdr: &Builder,
    box_path: &[(BoxIdx, Side)],
    cell_idx_by_coord: &HashMap<CellCoord, CellIdx>,
) -> Vec<Vec<CellIdx>> {
    let num_lanes = bdr.box_size_in_direction(box_path[0].1.direction());
    let lane_depth = bdr.box_size_in_direction(!box_path[0].1.direction());
    // The lanes are numbered going clockwise around the perimeter of each box
    let mut lanes = vec![Vec::<CellIdx>::new(); num_lanes];

    for &(box_idx, entry_side) in box_path {
        // Sanity check that the edges being stepped over don't change the number or depth of
        // the lanes.  This should have been caught in `generate_edges` (and an error
        // returned), but it doesn't hurt to double check.
        assert_eq!(bdr.box_size_in_direction(entry_side.direction()), num_lanes);
        assert_eq!(
            bdr.box_size_in_direction(!entry_side.direction()),
            lane_depth
        );
        // We know that the box size is valid, so add this box's cells to the lanes
        for (lane_idx, lane) in lanes.iter_mut().enumerate() {
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
                lane.push(cell_idx_by_coord[&CellCoord { box_idx, x, y }]);
            }
        }
    }

    lanes
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

/* UTILITY TYPES FOR SHAPE GENERATION */

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

/// An edge between two **boxes** (not to be confused with [`kuudos::shape::Edge`], which sits
/// between two **cells**).  If an [`Edge`] is only adjacent to one box/link, it will be put on the
/// 'right' side of the [`Edge`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Edge {
    vert_idx_top: VertIdx,
    vert_idx_bottom: VertIdx,

    connection_left: Option<EdgeConnection>, // External edges have nothing on their left side
    connection_right: EdgeConnection,        // Every edge always has something on its right side

    /// The number of cells that go down this `Edge`
    length: usize,
}

/// What can go to the side of an [`Edge`]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EdgeConnection {
    Box_(BoxIdx, Side),
    Link(LinkIdx, LinkSide),
}

/// The two sides of an [`Edge`] (assuming that the edge is pointing 'up')
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
