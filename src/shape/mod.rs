use std::collections::HashMap;
use std::hash::Hash;

use angle::Angle;
use itertools::Itertools;

use crate::image::{Elem, Image, RenderingOpts};
use crate::indexed_vec::{CellIdx, CellVec, IdxType, VertIdx, VertVec};
use crate::utils::Rect2;
use crate::{builder::Builder, V2Ext, V2};

pub mod examples;
mod to_image;

/// The shape of a sudoku as accepted by Kuudos
#[derive(Debug, Clone)]
pub struct Shape {
    /* SOLUTION CONSTRAINTS */
    /// How many different symbols are allowed in this `Shape`
    pub(crate) num_symbols: usize,
    /// Each group specifies a set of cells which can't contain two of the same symbol
    pub(crate) groups: Vec<Group>,

    /* DISPLAY DATA */
    /// The 2D coordinates of the vertices of this shape
    pub(crate) verts: VertVec<V2>,
    /// Each cell is a list of indices of vertices **in clockwise order**.  These will almost
    /// certainly have 4 vertices (and therefore sides), but leaving this generic allows Kuudos to
    /// handle other Sudoku shapes, such as those with hexagonal cells
    pub(crate) cells: CellVec<Vec<VertIdx>>,
    /// Extra elements drawn before the cells.  These can be used for things like rendering extra
    /// links between boxes.
    pub(crate) extra_elements: Vec<Elem>,
}

impl Shape {
    /// Returns the number of cells in this `Shape`
    pub fn num_cells(&self) -> usize {
        self.cells.len()
    }

    /// Write an empty sudoku grid to an SVG string
    pub fn svg_string_empty(&self, opts: &RenderingOpts, scaling: f32) -> String {
        self.svg_string(opts, scaling, &vec![None; self.num_cells()])
    }

    /// Generate an SVG string of a puzzle which has this `Shape`
    pub fn svg_string(&self, opts: &RenderingOpts, scale: f32, clues: &[Option<char>]) -> String {
        self.image(clues).svg_string(scale, opts)
    }

    pub fn image(&self, clues: &[Option<char>]) -> Image {
        to_image::gen_image_with_clues(self, clues)
    }

    /// Returns the number of groups which are shared between two cells
    pub(crate) fn num_groups_shared_between(
        &self,
        cell_idx1: CellIdx,
        cell_idx2: CellIdx,
    ) -> usize {
        self.groups
            .iter()
            .filter(|group| group.cells.contains(&cell_idx1) && group.cells.contains(&cell_idx2))
            .count()
    }

    /// Returns `true` if there is at least one box which contains both cells
    pub(crate) fn do_cells_share_box(&self, cell_idx1: CellIdx, cell_idx2: CellIdx) -> bool {
        self.groups.iter().any(|group| {
            group.is_box && group.cells.contains(&cell_idx1) && group.cells.contains(&cell_idx2)
        })
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

        for (cell_idx, verts) in self.cells.indexed_iter() {
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
    pub fn bbox(&self) -> Option<Rect2> {
        let bbox_from_verts = Rect2::bbox(self.verts.iter().copied());
        let bbox_from_extra_elems = Rect2::union_iter(self.extra_elements.iter().map(Elem::bbox));
        // Take the union of the `Option<Rect2>`s
        Rect2::union_option(bbox_from_verts, bbox_from_extra_elems)
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

        /* We number the cells as they would be read in a book; so a 4x4 grid would be numbered as
         * follows (numbers are written in hex for convenience):
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
         * Vertices are also labelled using this system.  However, there is a `box_width + 1` by
         * `box_height + 1` grid of vertices, so the resulting indices are different. */
        let grid_width = box_height * box_width;
        let cell_idx = |row: usize, col: usize| CellIdx::from_idx((row * grid_width) + col);

        /* GENERATE GROUPS */

        let mut groups = Vec::<Group>::with_capacity(grid_width * 3);
        // Rows
        for row in 0..grid_width {
            let cells = (0..grid_width).map(|col| cell_idx(row, col)).collect_vec();
            groups.push(Group::non_box_group(cells));
        }
        // Columns
        for col in 0..grid_width {
            let cells = (0..grid_width).map(|row| cell_idx(row, col)).collect_vec();
            groups.push(Group::non_box_group(cells));
        }
        // If `box_width` or `box_height` are 1, then the boxes will become either all rows or all
        // columns.  In either case, the boxes aren't needed
        if box_height != 1 && box_width != 1 {
            // Boxes (the grid of boxes will be `box_height` boxes wide, and `box_width` boxes
            // tall)
            for box_col in 0..box_height {
                for box_row in 0..box_width {
                    let mut cells = Vec::with_capacity(grid_width);
                    for sub_box_row in 0..box_height {
                        for sub_box_col in 0..box_width {
                            cells.push(cell_idx(
                                box_row * box_height + sub_box_row,
                                box_col * box_width + sub_box_col,
                            ));
                        }
                    }
                    groups.push(Group::box_group(cells));
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
            cells: cell_verts,
            extra_elements: vec![], // Square grids don't need any extra elements
        }
    }

    /// Create a star with a 2x2 diamond-shaped box under each point.  `Shape::star2x2(4, 0)` is
    /// equivalent to `Shape::square(2, 2)`.
    pub fn star2x2(num_points: usize, base_angle: impl Angle<f32> + Copy) -> Shape {
        // Create a new builder with 5-way rotational symmetry
        let mut builder = Builder::new(2, 2, num_points);
        // Create a new parallelogram box for the star's point
        let down = V2::DOWN;
        let rotated_down = builder.rotate_point_by_steps(down, 1);
        builder
            .add_box_parallelogram(V2::ZERO, down, rotated_down)
            // This unwrap is fine, because unless `num_points` is infinite the parallelogram boxes
            // can't be linear (and are therefore well-defined)
            .unwrap();
        // Rotate the built shape by the requested amount
        builder.rotate(base_angle);
        let (shape, _symmetry) = builder.build().unwrap();
        shape // For the time being, throw away the symmetry
    }
}

#[derive(Debug, Clone)]
pub struct Group {
    is_box: bool,
    pub(crate) cells: Vec<CellIdx>,
}

impl Group {
    pub(crate) fn box_group(cells: Vec<CellIdx>) -> Self {
        Self {
            is_box: true,
            cells,
        }
    }

    pub(crate) fn non_box_group(cells: Vec<CellIdx>) -> Self {
        Self {
            is_box: false,
            cells,
        }
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
    /// List of groups of cells which are symmetrically equivalent.
    ///
    /// **Invariant**: This must represent a partition over the cells - i.e. every [`CellIdx`] must
    /// appear exactly once in these equivalence classes.
    cell_equiv_classes: Vec<Vec<CellIdx>>,
}

impl Symmetry {
    /// A `Symmetry` which enforces no symmetry on a given [`Shape`] (i.e. each cell is in its own
    /// equivalence class).
    pub fn asymmetric(shape: &Shape) -> Self {
        Self {
            cell_equiv_classes: (0..shape.num_cells())
                .map(|idx| vec![CellIdx::from_idx(idx)])
                .collect_vec(),
        }
    }

    /// Creates a new `Symmetry` where cells who's assignments are equal will be in the same
    /// equivalence class.  The equivalence classes are always generated in a deterministic (but
    /// unspecified) order.
    pub fn from_cell_assignments<T: Eq + Hash>(cell_assignments: impl Iterator<Item = T>) -> Self {
        let mut cells_per_assignment = HashMap::<T, Vec<CellIdx>>::new();
        // Invariant: `cells_per_assignment[v]` contains the indices of all the cells yielded so
        // far which were assigned to `v`
        for (idx, v) in cell_assignments.enumerate() {
            cells_per_assignment
                .entry(v)
                .or_insert_with(Vec::new)
                .push(CellIdx::from_idx(idx));
        }
        // Sort the equivalence classes by their minimum cell, thus making sure their order is
        // deterministic
        let mut cell_equiv_classes = cells_per_assignment
            .into_iter()
            .map(|(_k, v)| v)
            .collect_vec();
        for class in &mut cell_equiv_classes {
            class.sort(); // Sort the cell indices in each equiv class
        }
        cell_equiv_classes.sort(); // Sort the equivalence classes themselves

        Self { cell_equiv_classes }
    }

    /// Each sub-list contains indices of cells which must be either all be clues or all be blank.
    pub fn equiv_classes(&self) -> &[Vec<CellIdx>] {
        &self.cell_equiv_classes
    }
}
