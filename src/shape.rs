use itertools::Itertools;
use vector2d::Vector2D;

/// Type alias for 2D floating point vectors (in the geometric sense, unlike [`Vec`])
type V2 = Vector2D<f32>;

/// The shape of a sudoku as accepted by Kuudos
#[derive(Debug, Clone)]
pub struct Shape {
    /* SOLUTION CONSTRAINTS */
    /// How many different symbols are allowed in this `Shape`
    pub(crate) num_symbols: usize,
    /// Each group specifies a set of cells which can't contain two of the same symbol
    pub(crate) groups: Vec<Vec<usize>>,

    /* DISPLAY DATA */
    /// The 2D coordinates of the vertices of this shape
    pub(crate) verts: Vec<V2>,
    /// Each face is a list of indices of vertices **in clockwise order**.  These will almost
    /// certainly have 4 vertices (and therefore sides), but leaving this generic allows Kuudos to
    /// handle other Sudoku shapes, such as those with hexagonal cells
    pub(crate) cell_verts: Vec<Vec<usize>>,
}

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
         *  0  1 | 2  3
         *       |
         *  4  5 | 6  7
         * ------+------
         *  8  9 | a  b
         *       |
         *  c  d | e  f
         *
         * Vertices are also labelled that way
         */
        let grid_width = box_height * box_width;

        let cell_idx = |row: usize, col: usize| (row * grid_width) + col;

        /* GENERATE GROUPS */

        let mut groups = Vec::<Vec<usize>>::with_capacity(grid_width * 3);
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
        let mut verts = Vec::<V2>::with_capacity((grid_width + 1) * (grid_width + 1));
        for y in 0..=grid_width {
            for x in 0..=grid_width {
                verts.push(V2::new(x as f32, y as f32));
            }
        }

        // Generate which vertices go round each cell
        let vert_idx = |x: usize, y: usize| y * (grid_width + 1) + x;
        let mut cell_verts = Vec::<Vec<usize>>::with_capacity(grid_width * grid_width);
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
}
