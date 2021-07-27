use std::collections::HashMap;

use angle::Angle;
use itertools::Itertools;

use crate::types::{CellIdx, CellVec, IdxType, VertIdx, VertVec};
use crate::{Builder, V2};

/// The shape of a sudoku as accepted by Kuudos
#[derive(Debug, Clone)]
pub struct Shape {
    /* SOLUTION CONSTRAINTS */
    /// How many different symbols are allowed in this `Shape`
    pub(crate) num_symbols: usize,
    /// Each group specifies a set of cells which can't contain two of the same symbol
    pub(crate) groups: Vec<Vec<CellIdx>>,

    /* DISPLAY DATA */
    /// The 2D coordinates of the vertices of this shape
    pub(crate) verts: VertVec<V2>,
    /// Each cell is a list of indices of vertices **in clockwise order**.  These will almost
    /// certainly have 4 vertices (and therefore sides), but leaving this generic allows Kuudos to
    /// handle other Sudoku shapes, such as those with hexagonal cells
    pub(crate) cell_verts: CellVec<Vec<VertIdx>>,
}

impl Shape {
    /// Returns the number of cells in this `Shape`
    pub fn num_cells(&self) -> usize {
        self.cell_verts.len()
    }

    /// Returns the number of groups which are shared between two cells
    pub(crate) fn num_groups_shared_between(
        &self,
        cell_idx1: CellIdx,
        cell_idx2: CellIdx,
    ) -> usize {
        self.groups
            .iter()
            .filter(|cells| cells.contains(&cell_idx1) && cells.contains(&cell_idx2))
            .count()
    }

    /// Generate the [`Edge`]s which appear in this `Shape`.
    ///
    /// # Panics
    ///
    /// Panics if more than two cells are adjacent to any [`Edge`]
    pub(crate) fn edges(&self) -> Vec<Edge> {
        // This maps (vert_idx_bottom, vert_idx_top) pairs to the Edges that they correspond to.
        // This allows the 2nd adjacent cell to figure out which edges it joins together
        let mut edges = HashMap::<(VertIdx, VertIdx), Edge>::new();

        for (cell_idx, verts) in self.cell_verts.indexed_iter() {
            // Iterate over the edges of this cell in clockwise order (i.e. with this cell on the
            // **right** side of each edge).
            for (&vert_idx_bottom, &vert_idx_top) in verts.iter().circular_tuple_windows() {
                // Check if the reverse of this edge has already been added
                if let Some(existing_edge) = edges.get_mut(&(vert_idx_top, vert_idx_bottom)) {
                    // Sanity check that this isn't the 3rd cell being attached to this edge
                    assert!(existing_edge.cell_idx_left.is_none());
                    existing_edge.cell_idx_left = Some(cell_idx);
                    continue;
                }
                // If the reversed version of this edge hasn't been found, check that we haven't
                // seen this forward edge before (which would result in an unprintable sudoku)
                assert!(!edges.contains_key(&(vert_idx_bottom, vert_idx_top)));
                // If this edge is new, then add it to the set with this cell as its right cell
                edges.insert(
                    (vert_idx_bottom, vert_idx_top),
                    Edge {
                        vert_idx_bottom,
                        vert_idx_top,
                        cell_idx_left: None,
                        cell_idx_right: cell_idx,
                    },
                );
            }
        }

        edges.into_iter().map(|(_k, v)| v).collect_vec()
    }

    /// Returns the bounding box of this `Shape` as a (min, max) pair of vectors
    pub(crate) fn bbox(&self) -> Option<(V2, V2)> {
        // If the shape has no vertices, then the bounding box isn't defined
        if self.verts.is_empty() {
            return None;
        }

        let mut min_x = f32::MAX;
        let mut min_y = f32::MAX;
        let mut max_x = f32::MIN;
        let mut max_y = f32::MIN;
        for v in self.verts.iter() {
            if v.x < min_x {
                min_x = v.x;
            }
            if v.y < min_y {
                min_y = v.y;
            }
            if v.x > max_x {
                max_x = v.x;
            }
            if v.y > max_y {
                max_y = v.y;
            }
        }
        Some((V2::new(min_x, min_y), V2::new(max_x, max_y)))
    }
}

/// Constructors
impl Shape {
    /// Creates the classic 9x9 sudoku grid that we know and love.
    pub fn classic() -> Shape {
        Self::square(3, 3)
    }

    /// Creates a `Shape` representing a square grid where each box has a given width and height.
    /// The total grid is therefore (box_width * box_height) by (box_width * box_height).
    pub fn square(box_width: usize, box_height: usize) -> Shape {
        assert_ne!(box_width, 0, "`box_width` can't be 0");
        assert_ne!(box_height, 0, "`box_height` can't be 0");

        /* We number the cells as they would be read in a book; so a 4x4 grid would be numbered
         * (numbers are written in hex for convenience):
         *
         * +------+------+
         * | 0  1 | 2  3 |
         * |      |      |
         * | 4  5 | 6  7 |
         * +------+------+
         * | 8  9 | a  b |
         * |      |      |
         * | c  d | e  f |
         * +------+------+
         *
         * Vertices are also labelled that way
         */
        let grid_width = box_height * box_width;
        let cell_idx = |row: usize, col: usize| CellIdx::from_idx((row * grid_width) + col);

        /* GENERATE GROUPS */

        let mut groups = Vec::<Vec<CellIdx>>::with_capacity(grid_width * 3);
        // Rows
        for row in 0..grid_width {
            groups.push((0..grid_width).map(|col| cell_idx(row, col)).collect_vec());
        }
        // Columns
        for col in 0..grid_width {
            groups.push((0..grid_width).map(|row| cell_idx(row, col)).collect_vec());
        }
        // If `box_width` or `box_height` are 1, then the boxes will become either all rows or all
        // columns.  In either case, the boxes aren't needed
        if box_height != 1 && box_width != 1 {
            // Boxes (the grid of boxes will be `box_height` boxes wide, and `box_width` boxes
            // tall)
            for box_col in 0..box_height {
                for box_row in 0..box_width {
                    let mut new_group = Vec::with_capacity(grid_width);
                    for sub_box_row in 0..box_height {
                        for sub_box_col in 0..box_width {
                            new_group.push(cell_idx(
                                box_row * box_height + sub_box_row,
                                box_col * box_width + sub_box_col,
                            ));
                        }
                    }
                    groups.push(new_group);
                }
            }
        }

        /* SHAPE */

        // Generate the vertices, labelled as they would be read in a book
        let mut verts = VertVec::<V2>::with_capacity((grid_width + 1) * (grid_width + 1));
        for y in 0..=grid_width {
            for x in 0..=grid_width {
                verts.push(V2::new(x as f32, y as f32));
            }
        }

        // Generate which vertices go round each cell
        let vert_idx = |x: usize, y: usize| VertIdx::from_idx(y * (grid_width + 1) + x);
        let mut cell_verts = CellVec::<Vec<VertIdx>>::with_capacity(grid_width * grid_width);
        for row in 0..grid_width {
            for col in 0..grid_width {
                cell_verts.push(
                    [(0, 0), (1, 0), (1, 1), (0, 1)]
                        .iter()
                        .map(|&(dx, dy)| vert_idx(col + dx, row + dy))
                        .collect_vec(),
                );
            }
        }

        Self {
            num_symbols: grid_width,
            groups,
            verts,
            cell_verts,
        }
    }

    /// Create a star with a 2x2 diamond-shaped box under each point.  `Shape::star2x2(4, 0)` is
    /// equivalent to `Shape::square(2, 2)`.
    pub fn star2x2(num_points: usize, base_angle: impl Angle<f32> + Copy) -> Shape {
        // Create a new builder with 5-way rotational symmetry
        let mut builder = Builder::new(2, 2, num_points);
        // Create a new parallelogram box for the star's point
        let down = V2::new(0.0, -1.0);
        let rotated_down = builder.rotate_point_by_steps(down, 1);
        builder.add_box_parallelogram(V2::new(0.0, 0.0), down, rotated_down);
        // Rotate the built shape by the requested amount
        builder.rotate(base_angle);
        // Build the shape
        let (shape, _symmetry) = builder.build().unwrap();
        shape // For the time being, throw away the symmetry
    }
}

/// An `Edge` is a line drawn on the grid between two vertices, and adjacent to either one or two
/// cells.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Edge {
    pub(crate) vert_idx_top: VertIdx,
    pub(crate) vert_idx_bottom: VertIdx,
    pub(crate) cell_idx_left: Option<CellIdx>,
    pub(crate) cell_idx_right: CellIdx,
}

/// A grouping of cells of a [`Shape`] into equivalence classes which determine the symmetry of the
/// clues in the resulting puzzle (i.e. in a puzzle grid with this `Symmetry`, all the cells in the
/// same equivalence class must either be all clues or all blanks).
#[derive(Debug, Clone)]
pub struct Symmetry {
    /// For each cell, this maps to an equivalence class index from `0..self.num_equiv_classes`
    cell_equiv_classes: Vec<usize>,
    /// How many unique equivalence classes this `Symmetry` contains
    num_equiv_classes: usize,
}

impl Symmetry {
    /// A `Symmetry` which enforces no symmetry on a given [`Shape`] (i.e. each cell is in its own
    /// equivalence class).
    pub fn asymmetric(shape: &Shape) -> Self {
        Self {
            cell_equiv_classes: (0..shape.num_cells()).collect_vec(),
            num_equiv_classes: shape.num_cells(),
        }
    }
}