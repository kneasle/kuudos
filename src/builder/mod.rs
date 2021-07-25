use std::f32::consts::PI;

use crate::{
    types::{BoxVec, VertIdx, VertVec},
    Shape, Symmetry, V2,
};

/// Code for converting a `Builder` into a [`Shape`]/[`Symmetry`] pair
mod build;

/// Re-export [`BuildError`] to the rest of the code.  This is then re-re-exported in `lib.rs`
pub use build::BuildError;

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
    /// The boxes added to the shape so far.  These arrays store the vertex indices going clockwise
    /// from the 'bottom left' of the box.  In reality, it's totally fine for people to rotate or
    /// reflect boxes (thus moving the 'bottom left' corner to some other direction), but this way
    /// we always have a consistent way to refer to the sides of boxes.  Additionally, the boxes
    /// may be reflected (possibly in the case of reflectional symmetry) which is fine for a
    /// `Builder` but not for a [`Shape`] (this is handled by [`Builder::build`]).
    boxes: BoxVec<[VertIdx; 4]>,
    /// The number given to the next vertex equivalence class
    next_equiv_class: usize,
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
            boxes: BoxVec::new(),

            next_equiv_class: 0,
        }
    }

    /* Add Boxes */

    /// Add a new box to the grid (in the shape of a parallelogram)
    pub fn add_box_parallelogram(&mut self, bottom_left_corner: V2, up: V2, right: V2) {
        // `up` and `right` refer to the sides of the cells, so we scale up to the size of the box
        let box_up = up * self.box_height as f32;
        let box_right = right * self.box_width as f32;
        // Compute the unrotated coordinates of the four corners of the new box(es)
        let unrotated_bl = bottom_left_corner;
        let unrotated_tl = bottom_left_corner + box_up;
        let unrotated_br = bottom_left_corner + box_right;
        let unrotated_tr = bottom_left_corner + box_up + box_right;
        // For each symmetry element
        for i in 0..self.rotational_symmetry_factor as isize {
            let new_box = [
                self.get_vert_at(self.rotate_point_by_steps(unrotated_bl, i)),
                self.get_vert_at(self.rotate_point_by_steps(unrotated_tl, i)),
                self.get_vert_at(self.rotate_point_by_steps(unrotated_tr, i)),
                self.get_vert_at(self.rotate_point_by_steps(unrotated_br, i)),
            ];
            self.boxes.push(new_box);
        }
    }

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

    /* Modifiers (i.e. functions which mutate the existing shape) */

    /// Rotates the whole shape around the origin (i.e. the centre of symmetry) by some angle in
    /// radians
    pub fn rotate(&mut self, angle: f32) {
        // Rotate all the vertices in-place
        for vert in self.verts.iter_mut() {
            vert.position = rotate_vec(vert.position, angle);
        }
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
        let angle = PI * 2.0 * num_factors as f32 / self.rotational_symmetry_factor as f32;
        rotate_vec(v, angle)
    }

    /* Conversion to `Shape` */

    /// Converts this `Builder` into a [`Shape`] and the associated [`Symmetry`]
    pub fn build(self) -> Result<(Shape, Symmetry), BuildError> {
        build::build(self)
    }
}

/// Rotates a vector **clockwise** by an angle in radians
fn rotate_vec(v: V2, angle: f32) -> V2 {
    let sin = angle.sin();
    let cos = angle.cos();
    // Rotation **clockwise** corresponds to multiplication by the following matrix (which looks
    // like the classic anti-clockwise matrix because our y-axis goes down where the one in maths
    // goes up):
    // | cos(angle)  -sin(angle) |
    // | sin(angle)   cos(angle) |
    V2 {
        x: v.x * cos - v.y * sin,
        y: v.x * sin + v.y * cos,
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

#[cfg(test)]
mod tests {
    use std::f32::consts::PI;

    use crate::V2;

    #[test]
    fn angle_between() {
        fn check(v1: V2, v2: V2, exp_angle: f32) {
            assert_eq!(super::angle_between(v1, v2), exp_angle);
        }

        // The angle from any direction to (any +ve scalar multiple of itself) is 0
        check(V2::new(0.0, -1.0), V2::new(0.0, -1.0), 0.0);
        check(V2::new(0.0, 1.0), V2::new(0.0, 1.0), 0.0);
        check(V2::new(1.0, 1.0), V2::new(2.5, 2.5), 0.0);
        // The angle from any direction to any -ve scalar multiple of itself is +PI
        check(V2::new(0.0, 1.0), V2::new(0.0, -1.0), PI);
        check(V2::new(1.0, -1.0), V2::new(-1.0, 1.0), PI);

        check(V2::new(0.0, -1.0), V2::new(2.5, 0.0), PI / 2.0);
        check(V2::new(0.0, 1.0), V2::new(2.5, 0.0), -PI / 2.0);
        check(V2::new(0.0, 1.0), V2::new(-1.0, -1.0), 0.75 * PI);
        check(V2::new(0.0, 1.0), V2::new(1.0, -1.0), -0.75 * PI);
    }
}
