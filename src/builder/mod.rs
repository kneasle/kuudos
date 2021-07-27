use std::ops::Not;

use crate::{
    types::{BoxIdx, BoxVec, SymmIdx, SymmVec, VertIdx, VertVec},
    utils::rotate_vec,
    Shape, Symmetry, V2,
};

use angle::{Angle, Deg};
use itertools::Itertools;

mod gen_shape;

/// Re-export [`BuildError`] to the rest of the code.  This is then re-re-exported in `lib.rs`
pub use gen_shape::BuildError;

/// A struct for programmatically building complex `Shapes`
#[derive(Debug, Clone)]
pub struct Builder {
    /// The number of 'columns' in each box of this grid
    box_width: usize,
    /// The number of 'rows' in each box of this grid
    box_height: usize,

    /// The current rotational symmetry factor of the puzzle.  Each object added will be duplicated
    /// and rotated this many times around the origin.  This will also determine the [`Symmetry`]
    /// of the cells added (see [`Symmetry`] for more docs).  This can be changed at any point and
    /// the existing shapes will keep the symmetry they were created with.
    rotational_symmetry_factor: usize,

    /// The coordinates of the vertices added so far
    verts: VertVec<Vert>,
    /// The number given to the next vertex equivalence class
    next_equiv_class: usize,

    /// The boxes added to the shape so far.  These arrays store the vertex indices going clockwise
    /// from the 'bottom left' of the box.  In reality, it's totally fine for people to rotate or
    /// reflect boxes (thus moving the 'bottom left' corner to some other direction), but this way
    /// we always have a consistent way to refer to the sides of boxes.  Additionally, the boxes
    /// may be reflected (possibly in the case of reflectional symmetry) which is fine for a
    /// `Builder` but not for a [`Shape`] (this is handled by [`Builder::build`]).
    boxes: BoxVec<Box_>,
    /// Each `Vec<BoxIdx>` corresponds to a sequence of [`Box_`]es which are 'equivalent' according
    /// to rotational symmetry.
    box_equiv_classes: SymmVec<Vec<BoxIdx>>,
}

impl Builder {
    /* Constructors */

    /// Creates an empty `Builder` with no symmetry pattern
    pub fn new(box_width: usize, box_height: usize, rotational_symmetry_factor: usize) -> Self {
        Self {
            box_width,
            box_height,
            rotational_symmetry_factor,

            verts: VertVec::new(),
            next_equiv_class: 0,

            boxes: BoxVec::new(),
            box_equiv_classes: SymmVec::new(),
        }
    }

    /* Add Boxes */

    /// Adds a new square box with a given position, size and location
    pub fn add_box_square(
        &mut self,
        centre: V2,
        cell_size: f32,
        angle: impl Angle<f32> + Copy,
    ) -> BoxIdx {
        // Vector from the centre of the unrotated box to its top-right corner
        let half_diagonal =
            V2::new(self.box_width as f32, -(self.box_height as f32)) * cell_size / 2.0;
        let bottom_left = centre - rotate_vec(half_diagonal, angle);
        self.add_box_parallelogram(
            bottom_left,
            // Up and right directions are also scaled and rotated
            rotate_vec(V2::new(0.0, -1.0) * cell_size, angle),
            rotate_vec(V2::new(1.0, 0.0) * cell_size, angle),
        )
    }

    /// Add a new box in the shape of a parallelogram
    pub fn add_box_parallelogram(&mut self, bottom_left_corner: V2, up: V2, right: V2) -> BoxIdx {
        // `up` and `right` refer to the sides of the cells, so we scale up to the size of the box
        let box_up = up * self.box_height as f32;
        let box_right = right * self.box_width as f32;
        // Compute the unrotated coordinates of the four corners of the new box(es)
        let unrotated_bl = bottom_left_corner;
        let unrotated_tl = bottom_left_corner + box_up;
        let unrotated_br = bottom_left_corner + box_right;
        let unrotated_tr = bottom_left_corner + box_up + box_right;
        // Create a new box (or rather rotated set of boxes) with these vertices
        self.add_box(unrotated_bl, unrotated_tl, unrotated_tr, unrotated_br)
    }

    /// Connect the edges of two (different) boxes with a new box
    #[must_use]
    pub fn connect_boxes(
        &mut self,
        top_box_idx: BoxIdx,
        top_edge: Side,
        bottom_box_idx: BoxIdx,
        bottom_edge: Side,
    ) -> Option<BoxIdx> {
        // It doesn't make sense to connect two sides of the same box
        assert_ne!(top_box_idx, bottom_box_idx);
        let top_box = self.boxes.get(top_box_idx)?;
        let bottom_box = self.boxes.get(bottom_box_idx)?;
        let (vert_top_right, vert_top_left) = top_box.get_edge_verts(top_edge);
        let (vert_bottom_left, vert_bottom_right) = bottom_box.get_edge_verts(bottom_edge);
        Some(self.add_box(
            self.verts.get(vert_bottom_left).unwrap().position,
            self.verts.get(vert_top_left).unwrap().position,
            self.verts.get(vert_top_right).unwrap().position,
            self.verts.get(vert_bottom_right).unwrap().position,
        ))
    }

    /// Adds a new box to the shape, given 4 arbitrary vertices
    pub fn add_box(
        &mut self,
        bottom_left: V2,
        top_left: V2,
        top_right: V2,
        bottom_right: V2,
    ) -> BoxIdx {
        // Generate each rotation of this box, storing that they are 'equivalent'
        let equiv_class_idx = self.box_equiv_classes.next_idx();
        let mut box_idxs = Vec::<BoxIdx>::with_capacity(self.rotational_symmetry_factor);
        for i in 0..self.rotational_symmetry_factor {
            let new_box = Box_ {
                vert_idxs: [
                    self.get_vert_at(self.rotate_point_by_steps(bottom_left, i as isize)),
                    self.get_vert_at(self.rotate_point_by_steps(top_left, i as isize)),
                    self.get_vert_at(self.rotate_point_by_steps(top_right, i as isize)),
                    self.get_vert_at(self.rotate_point_by_steps(bottom_right, i as isize)),
                ],
                equiv_class: equiv_class_idx,
                rotation_within_equiv_class: i,
            };
            let new_box_idx = self.boxes.push(new_box);
            box_idxs.push(new_box_idx);
        }
        // Add the new equivalence class, checking that all the boxes are using the correct index
        let idx_of_original_box = box_idxs[0];
        assert_eq!(self.box_equiv_classes.push(box_idxs), equiv_class_idx);
        idx_of_original_box
    }

    /* Modifiers (i.e. functions which mutate the existing shape) */

    /// Rotates the whole shape around the origin (i.e. the centre of symmetry) by some angle in
    /// radians
    pub fn rotate(&mut self, angle: impl Angle<f32> + Copy) {
        // Rotate all the vertices in-place
        for vert in self.verts.iter_mut() {
            vert.position = rotate_vec(vert.position, angle);
        }
    }

    /* Getters */

    /// Gets the [`BoxIdx`] of the copy of a box after being rotated by some number of rotation
    /// steps.
    pub fn rotational_copy_of(&self, box_idx: BoxIdx, rotation_steps: isize) -> Option<BoxIdx> {
        let box_ = self.boxes.get(box_idx)?;
        let boxes_in_equiv_class = self.box_equiv_classes.get(box_.equiv_class).unwrap();

        let num_rotation_steps = boxes_in_equiv_class.len();
        let unwrapped_idx = box_.rotation_within_equiv_class as isize + rotation_steps;
        let wrapped_idx = (unwrapped_idx % num_rotation_steps as isize
            + num_rotation_steps as isize)
            % num_rotation_steps as isize;

        assert!(wrapped_idx >= 0);
        Some(boxes_in_equiv_class[wrapped_idx as usize])
    }

    /* Setters */

    /*
    /// Sets the rotational symmetry factor (around the origin) of this `Builder`.  Until this is
    /// next changed, all the newly added boxes will
    /// 1. Be duplicated (up to) this many times
    /// 2. Will generate with this symmetry factor in the resulting grid
    pub fn set_rotational_symmetry(&mut self, rotational_symmetry_factor: usize) {
        assert_eq!(self.rotational_symmetry_factor, 1);
        self.rotational_symmetry_factor = rotational_symmetry_factor;
    }
    */

    /* Helpers */

    /// Rotate a point by some number of rotation steps
    pub fn rotate_point_by_steps(&self, v: V2, num_factors: isize) -> V2 {
        let angle = Deg(360.0) * num_factors as f32 / self.rotational_symmetry_factor as f32;
        rotate_vec(v, angle)
    }

    /* Conversion to `Shape` */

    /// Converts this `Builder` into a [`Shape`] and the associated [`Symmetry`]
    pub fn build(self) -> Result<(Shape, Symmetry), BuildError> {
        gen_shape::build(self)
    }

    /* Internal helpers */

    /// Gets the vertex at a given position, creating a new vertex if necessary.
    /// This also creates symmetric versions of the new vertex and labels them as equivalent.
    fn get_vert_at(&mut self, new_pos: V2) -> VertIdx {
        let existing_vert_idx = self
            .verts
            .idx_of_first(|vert| (vert.position - new_pos).length() < 0.0001);

        // If this vertex already exists, then return its index without creating a new set of
        // vertices
        if let Some(idx) = existing_vert_idx {
            return idx;
        }

        // If the new vertex is at the origin, then all its rotations will be the same and so we
        // add it once
        if new_pos.length_squared() < 0.00001 * 0.00001 {
            return self.verts.push(Vert {
                position: V2::new(0.0, 0.0),
                equiv_class: VertEquivClass::Centre,
            });
        }

        // Because (currently) the symmetry factor must stay constant, if this vertex doesn't
        // already exist, then all of its rotations must also not already exist.  Therefore, we
        // can create all the rotational copies without checking if they already exist.
        let idx_of_new_vert = self.verts.next_idx();
        let equiv_class = self.fresh_equiv_class();
        for i in 0..self.rotational_symmetry_factor {
            self.verts.push(Vert {
                position: self.rotate_point_by_steps(new_pos, i as isize),
                equiv_class: VertEquivClass::NonCentre {
                    equiv_class,
                    equiv_rotation: i,
                },
            });
        }
        idx_of_new_vert
    }

    /// Generates a fresh equivalence class name - i.e. one which hasn't been used before
    fn fresh_equiv_class(&mut self) -> usize {
        let class = self.next_equiv_class;
        self.next_equiv_class += 1;
        class
    }
}

#[derive(Debug, Clone)]
struct Box_ {
    /// The indices of the vertices of this corner.  These are defined to be listed in clockwise
    /// order with the 0th element being the 'bottom-left' corner.
    vert_idxs: [VertIdx; 4],
    equiv_class: SymmIdx,
    rotation_within_equiv_class: usize,
}

impl Box_ {
    /// Swaps the direction of the vertices (so they go from clockwise order to anti-clockwise or
    /// vice versa).
    fn swap_direction(&mut self) {
        self.vert_idxs.swap(1, 3);
    }

    /// Gets the two vertices on either side of an edge of this box.  The pair of vertices are in
    /// clockwise order.
    fn get_edge_verts(&self, side: Side) -> (VertIdx, VertIdx) {
        let edge_idx = match side {
            Side::Left => 0,
            Side::Top => 1,
            Side::Right => 2,
            Side::Bottom => 3,
        };
        self.vert_idxs
            .iter()
            .copied()
            .circular_tuple_windows()
            .nth(edge_idx)
            .unwrap()
    }
}

#[derive(Debug, Clone, Copy)]
struct Vert {
    /// The position of this vertex in 2D space
    position: V2,
    /// This vertex's position within the symmetry of the [`Shape`] being built
    equiv_class: VertEquivClass,
}

#[derive(Debug, Clone, Copy)]
enum VertEquivClass {
    /// The vertex is at the origin, and therefore is simultaneously in every rotation
    Centre,
    /// The vertex isn't at the origin, and therefore belongs to a specific rotation of an
    /// equivalence that it shares with other vertices
    NonCentre {
        /// Which equivalence class this vertex belongs to
        equiv_class: usize,
        /// What rotation this vertex has within its equivalence class.  This allows vertices in the
        /// same equivalence class to be differentiated when computing the [`Symmetry`] of the
        /// resulting [`Shape`]
        equiv_rotation: usize,
    },
}

/// The possible directions an edge can face
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub enum Direction {
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
pub enum Side {
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
