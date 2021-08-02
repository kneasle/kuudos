use std::ops::Not;

use crate::{
    indexed_vec::{BoxIdx, BoxVec, LinkIdx, LinkVec, SymmIdx, SymmVec, VertIdx, VertVec},
    utils, Shape, Symmetry, V2Ext, V2,
};

use angle::{Angle, Deg};
use itertools::Itertools;

mod gen_shape;
mod svg;

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
    verts: VertVec<V2>,

    /// Each link joins two edges which don't touch each other.  These allow us to space out the
    /// boxes further, and allows for many more shapes
    edge_links: LinkVec<EdgeLink>,

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

            edge_links: LinkVec::new(),

            boxes: BoxVec::new(),
            box_equiv_classes: SymmVec::new(),
        }
    }

    /* Add Boxes */

    /// Adds a new square box with a given centre, size and location
    pub fn add_box_square(
        &mut self,
        centre: V2,
        cell_size: f32,
        angle: impl Angle<f32> + Copy,
    ) -> BoxIdx {
        // Vector from the centre of the unrotated box to its top-right corner
        let half_diagonal =
            V2::new(self.box_width as f32, -(self.box_height as f32)) * cell_size / 2.0;
        let bottom_left = centre - half_diagonal.rotate(angle);
        // Delegate to `add_box_square_by_corner`
        self.add_box_square_by_corner(bottom_left, cell_size, angle)
    }

    /// Adds a new square box with a given centre, size and location
    pub fn add_box_square_by_corner(
        &mut self,
        bottom_left: V2,
        cell_size: f32,
        angle: impl Angle<f32> + Copy,
    ) -> BoxIdx {
        self.add_box_parallelogram(
            bottom_left,
            // Up and right directions are also scaled and rotated
            V2::UP.rotate(angle) * cell_size,
            V2::RIGHT.rotate(angle) * cell_size,
        )
        .unwrap() // Unwrapping here is fine, because square boxes are always well-defined
    }

    /// Add a new box in the shape of a parallelogram
    pub fn add_box_parallelogram(
        &mut self,
        bottom_left_corner: V2,
        up: V2,
        right: V2,
    ) -> Result<BoxIdx, BoxAddError> {
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

    /// Add a new box in the shape of a parallelogram, rotating the box so that the side originally
    /// labelled [`Side::Bottom`] will get the given label
    pub fn add_box_parallelogram_rotated(
        &mut self,
        bottom_left_corner: V2,
        up: V2,
        right: V2,
        map_bottom_to_side: Side,
    ) -> Result<BoxIdx, BoxAddError> {
        let (width, height) = match map_bottom_to_side.direction() {
            Direction::Horizontal => (self.box_width, self.box_height),
            // Swap width and height if the box is being rotated by 90 or 270 degrees
            Direction::Vertical => (self.box_height, self.box_width),
        };
        // `up` and `right` refer to the sides of the cells, so we scale up to the size of the box
        let box_right = right * width as f32;
        let box_up = up * height as f32;
        // Compute the unrotated coordinates of the four corners of the new box(es)
        let unrotated_bl = bottom_left_corner;
        let unrotated_tl = bottom_left_corner + box_up;
        let unrotated_br = bottom_left_corner + box_right;
        let unrotated_tr = bottom_left_corner + box_up + box_right;
        // Rotate the vertex slice
        let mut slice = [unrotated_bl, unrotated_tl, unrotated_tr, unrotated_br];
        slice.rotate_right((map_bottom_to_side.index() + 4 - Side::Bottom.index()) % 4);
        // Create a new box (or rather rotated set of boxes) with these vertices
        self.add_box_from_slice(slice)
    }

    /// Connect the edges of two (different) boxes with a new box
    pub fn connect_edges_with_box(
        &mut self,
        top_box_idx: BoxIdx,
        top_edge: Side,
        bottom_box_idx: BoxIdx,
        bottom_edge: Side,
    ) -> Result<BoxIdx, BoxAddError> {
        assert_ne!(top_box_idx, bottom_box_idx); // You can't connect two sides of the same box

        // Extract vertex locations
        let (vert_top_right, vert_top_left) = self.vert_positions_of_edge(top_box_idx, top_edge)?;
        let (vert_bottom_left, vert_bottom_right) =
            self.vert_positions_of_edge(bottom_box_idx, bottom_edge)?;

        // Add the box
        self.add_box(
            vert_bottom_left,
            vert_top_left,
            vert_top_right,
            vert_bottom_right,
        )
    }

    /// Adds a new box which extends outwards from a given edge.  The new box always has the same
    /// orientation as the box being extruded.
    pub fn extrude_edge(&mut self, box_idx: BoxIdx, side: Side) -> Result<BoxIdx, BoxAddError> {
        self.extrude_edge_with_opts(box_idx, side, 1.0, None)
    }

    /// Adds a new box which extends outwards from a given edge, with extra options.  The new box
    /// always has the same orientation as the box being extruded.
    pub fn extrude_edge_with_opts(
        &mut self,
        box_idx: BoxIdx,
        side: Side,
        extruded_cell_height: f32,
        link_height: Option<f32>,
    ) -> Result<BoxIdx, BoxAddError> {
        // Extract the vertices of the side being extruded
        let box_ = self
            .boxes
            .get(box_idx)
            .ok_or(BoxAddError::InvalidBoxIdx(box_idx))?;
        let (vert_idx1, vert_idx2) = box_.get_edge_verts(side);
        let v1 = self.verts[vert_idx1];
        let v2 = self.verts[vert_idx2];

        // Compute (indirectly) the up and right directions for the new box
        let edge_direction = v2 - v1; // Vector pointing down the full length of the edge
        let normal = -edge_direction.normal().normalise(); // Points outwards from the edge
        let edge_direction_per_cell =
            edge_direction / self.box_size_in_direction(side.direction()) as f32;

        // Create a new box at the right distance from the source box
        let bottom_left_corner = v1 + normal * link_height.unwrap_or(0.0);
        let new_box = self
            .add_box_parallelogram_rotated(
                bottom_left_corner,
                normal * extruded_cell_height,
                edge_direction_per_cell,
                side.opposite(), // Rotate the box so that `Bottom` becomes `side.opposite()`
            )
            // This unwrap is safe because the resulting box will always be well-defined
            .unwrap();
        // Add an edge link, if needed
        if link_height.is_some() {
            self.link_edges(
                box_idx,
                side,
                new_box,
                side.opposite(),
                EdgeLinkStyle::Linear,
            )
            .unwrap();
        }
        // Return the new box
        Ok(new_box)
    }

    /// Adds a new box to the shape, given 4 arbitrary vertices
    pub fn add_box(
        &mut self,
        bottom_left: V2,
        top_left: V2,
        top_right: V2,
        bottom_right: V2,
    ) -> Result<BoxIdx, BoxAddError> {
        self.add_box_from_slice([bottom_left, top_left, top_right, bottom_right])
    }

    /// Adds a new box to the shape, given a slice of 4 arbitrary vertices
    pub fn add_box_from_slice(&mut self, verts: [V2; 4]) -> Result<BoxIdx, BoxAddError> {
        let rotate_direction = classify_box(verts)?;
        // Generate each rotation of this box, storing that they are 'equivalent'.  Rotating a box
        // doesn't change the relative positions of the vertices, so the validity/rotation check
        // doesn't need to be performed each time.
        let equiv_class_idx = self.box_equiv_classes.next_idx();
        let mut box_idxs = Vec::<BoxIdx>::with_capacity(self.rotational_symmetry_factor);
        for i in 0..self.rotational_symmetry_factor {
            let new_box = Box_ {
                vert_idxs: [
                    self.get_vert_at(self.rotate_point_by_steps(verts[0], i as isize)),
                    self.get_vert_at(self.rotate_point_by_steps(verts[1], i as isize)),
                    self.get_vert_at(self.rotate_point_by_steps(verts[2], i as isize)),
                    self.get_vert_at(self.rotate_point_by_steps(verts[3], i as isize)),
                ],
                equiv_class: equiv_class_idx,
                rotation_within_equiv_class: i,
                rotate_direction,
            };
            let new_box_idx = self.boxes.push(new_box);
            box_idxs.push(new_box_idx);
        }
        // Add the new equivalence class, checking that all the boxes are using the correct index
        let idx_of_original_box = box_idxs[0];
        assert_eq!(self.box_equiv_classes.push(box_idxs), equiv_class_idx);
        Ok(idx_of_original_box)
    }

    /* Add edge links */

    /// Add an explicit link between two (non-adjacent) edges
    pub fn link_edges(
        &mut self,
        top_box_idx: BoxIdx,
        top_side: Side,
        bottom_box_idx: BoxIdx,
        bottom_side: Side,
        style: EdgeLinkStyle,
    ) -> Option<LinkIdx> {
        // Compute the vertices on either side of the edge
        let (vert_pos_tr, vert_pos_tl) = self.vert_positions_of_edge(top_box_idx, top_side).ok()?;
        let (vert_pos_bl, vert_pos_br) = self
            .vert_positions_of_edge(bottom_box_idx, bottom_side)
            .ok()?;

        // Add all symmetric versions of this link
        let idx_of_original_link = self.edge_links.next_idx();
        for i in 0..self.rotational_symmetry_factor {
            // Cheeky helper function to get the index of a rotated vertex
            let mut get_rotated_vert_idx =
                |v: V2| self.get_vert_at(self.rotate_point_by_steps(v, i as isize));

            let new_link = EdgeLink {
                vert_idx_top_left: get_rotated_vert_idx(vert_pos_tl),
                vert_idx_top_right: get_rotated_vert_idx(vert_pos_tr),
                vert_idx_bottom_left: get_rotated_vert_idx(vert_pos_bl),
                vert_idx_bottom_right: get_rotated_vert_idx(vert_pos_br),
                style,
            };
            self.edge_links.push(new_link);
        }
        Some(idx_of_original_link)
    }

    /* Modifiers (i.e. functions which mutate the existing shape) */

    /// Rotates the whole shape around the origin (i.e. the centre of symmetry) by some angle in
    /// radians
    pub fn rotate(&mut self, angle: impl Angle<f32> + Copy) {
        // Rotate all the vertices in-place
        for vert in self.verts.iter_mut() {
            *vert = vert.rotate(angle);
        }
    }

    /* Getters */

    /// Gets the [`BoxIdx`] of the copy of a box after being rotated by some number of rotation
    /// steps.
    pub fn rotational_copy_of(&self, box_idx: BoxIdx, rotation_steps: isize) -> Option<BoxIdx> {
        let box_ = self.boxes.get(box_idx)?;
        let boxes_in_equiv_class = &self.box_equiv_classes[box_.equiv_class];

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
        v.rotate(angle)
    }

    /* Conversions */

    /// Generates a debug-able SVG string representing the current state of `self`
    pub fn as_svg(&self, scaling: f32) -> String {
        svg::gen_svg(self, scaling)
    }

    /// Converts this `Builder` into a [`Shape`] and the associated [`Symmetry`]
    pub fn build(self) -> Result<(Shape, Symmetry), BuildError> {
        gen_shape::build(self)
    }

    /* Internal helpers */

    /// Gets the two vertex positions on either end of some [`Side`] of a [`Box_`], in clockwise
    /// order.
    fn vert_positions_of_edge(&self, box_idx: BoxIdx, side: Side) -> Result<(V2, V2), BoxAddError> {
        let (vert_idx1, vert_idx2) = self
            .boxes
            .get(box_idx)
            .ok_or(BoxAddError::InvalidBoxIdx(box_idx))?
            .get_edge_verts(side);
        Ok((self.verts[vert_idx1], self.verts[vert_idx2]))
    }

    /// Gets the vertex at a given position, creating a new vertex if necessary.
    /// This also creates symmetric versions of the new vertex and labels them as equivalent.
    fn get_vert_at(&mut self, new_pos: V2) -> VertIdx {
        let existing_vert_idx = self
            .verts
            .idx_of_first(|vert| (*vert - new_pos).length() < 0.0001);
        // Add a new vertex if there isn't one already at this position
        existing_vert_idx.unwrap_or_else(|| self.verts.push(new_pos))
    }

    /// Gets the number of cells along any [`Edge`] of a box in a given [`Direction`].
    fn box_size_in_direction(&self, direction: Direction) -> usize {
        match direction {
            Direction::Horizontal => self.box_width,
            Direction::Vertical => self.box_height,
        }
    }
}

/// Helper function that takes 4 vertices, and classifies the box that joins those vertices (in the
/// order specified)
fn classify_box(verts: [V2; 4]) -> Result<RotateDirection, BoxAddError> {
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
        fn classify(a: impl Angle<f32>) -> CornerType {
            let angle = a.to_deg().0;
            if angle.abs() < 0.001 || angle.abs() > 360.0 - 0.001 {
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

    // Iterate over the positions of all consecutive triple of vertices (i.e. pairs of edges
    // which join at a vertex), and classify that vertex according to its direction
    let directions = verts
        .iter()
        .circular_tuple_windows()
        .map(|(v1, v2, v3)| {
            // The edges go v1 -> v2 and v2 -> v3
            let d1 = v2 - v1;
            let d2 = v3 - v2;
            // Classify the middle corner (v2) according to its angle
            CornerType::classify(utils::angle_between(d1, d2))
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
        // anticlockwise
        (4, 0, 0) => Ok(RotateDirection::AntiClockwise),
        // If all corners go clockwise, then this box doesn't need modification
        (0, 0, 4) => Ok(RotateDirection::Clockwise),

        // If two corners go each way, then the box must self-intersect
        (2, 0, 2) => Err(BoxAddError::SelfIntersecting),
        // If there are no linear corners, then the only remaining case is that the box
        // isn't strictly convex
        (_, 0, _) => Err(BoxAddError::NonConvex),
        // If there are any linear corners, then we class the whole box as linear. Any case where
        // `num_linear == 0` would be caught by the non-convex case (or higher cases)
        (_, _, _) => Err(BoxAddError::LinearCorner),
    }
}

//////////////////////////
// ELEMENTS OF BUILDERS //
//////////////////////////

/* EDGE LINKS */

/// A link between two non-overlapping edges.  The left/right-ness of vertices doesn't matter - the
/// important thing is that the lefts and rights are consistent at either end (or else the edge
/// link will cause a twist).
#[derive(Debug, Clone)]
struct EdgeLink {
    // Starting edge
    vert_idx_top_left: VertIdx,
    vert_idx_top_right: VertIdx,
    // End edge
    vert_idx_bottom_left: VertIdx,
    vert_idx_bottom_right: VertIdx,
    // What shape this link should take
    style: EdgeLinkStyle,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub enum EdgeLinkStyle {
    /// Draw a circular arc between the edges
    Arc,
    /// Draw a straight line between the edges
    Linear,
    /// Don't draw anything.  This is only useful if (for example) you're making a 2D net of a 3D
    /// shape and want to encode the connections without cluttering the drawing
    Hidden,
}

/// The sides of an [`EdgeLink`] that can connect to boxes
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub enum LinkSide {
    Top,
    Bottom,
}

impl From<LinkSide> for Side {
    fn from(side: LinkSide) -> Self {
        match side {
            LinkSide::Top => Side::Top,
            LinkSide::Bottom => Side::Bottom,
        }
    }
}

impl Not for LinkSide {
    type Output = Self;

    fn not(self) -> Self::Output {
        match self {
            LinkSide::Top => LinkSide::Bottom,
            LinkSide::Bottom => LinkSide::Top,
        }
    }
}

/* BOX */

/// A box of cells.  This is guaranteed to be non-intersecting (and all the box creation methods
/// check this)
#[derive(Debug, Clone)]
struct Box_ {
    /// The indices of the vertices of this corner.  These are defined to be listed in clockwise
    /// order with the 0th element being the 'bottom-left' corner.
    vert_idxs: [VertIdx; 4],
    equiv_class: SymmIdx,
    rotation_within_equiv_class: usize,
    /// In which direction the vertices are specified.  We maintain the invariant that symmetric
    /// versions of boxes should have symmetric edge labellings, so in the case of reflectional
    /// symmetry this requires us to create clockwise/anticlockwise pairs of boxes:
    /// ```text
    ///    (clockwise)     |   (anticlockwise)
    ///                    |
    ///    0 ------- 1     |     1 ------- 0
    ///    |  ---->  |     |     |  <----  |
    ///    | ^     | |     |     | |     ^ |
    ///    | |     V |     |     | V     | |
    ///    |  <----  |     |     |  ---->  |
    ///    3 ------- 2     |     2 ------- 3
    ///                    |
    ///               MIRROR LINE
    /// ```
    rotate_direction: RotateDirection,
}

impl Box_ {
    /// Swaps the direction of the vertices (so they go from clockwise order to anti-clockwise or
    /// vice versa), preserving the lengths of each side.
    fn swap_direction(&mut self) {
        self.vert_idxs.swap(0, 1); // Reverse left edge
        self.vert_idxs.swap(2, 3); // Reverse right edge
        self.rotate_direction = !self.rotate_direction;
    }

    /// Gets the two vertices on either side of an edge of this box.  The pair of vertices are in
    /// clockwise order.
    fn get_edge_verts(&self, side: Side) -> (VertIdx, VertIdx) {
        self.vert_idxs
            .iter()
            .copied()
            .circular_tuple_windows()
            .nth(side.index())
            .unwrap()
    }
}

/// The possible ways that generating a box can fail
#[derive(Debug, Clone, Copy)]
pub enum BoxAddError {
    /// The newly created box would intersect with itself
    SelfIntersecting,
    /// The newly created box has two adjacent edges forming a straight line
    LinearCorner,
    /// The newly created box is not strictly convex
    NonConvex,
    /// A [`BoxIdx`] was given that doesn't correspond to an exiting box
    InvalidBoxIdx(BoxIdx),
}

////////////////////////
// UTILITY DATA TYPES //
////////////////////////

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

/// The possible directions of a rotation (clockwise/anti-clockwise)
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub enum RotateDirection {
    Clockwise,
    AntiClockwise,
}

impl Not for RotateDirection {
    type Output = Self;

    fn not(self) -> Self::Output {
        match self {
            Self::Clockwise => Self::AntiClockwise,
            Self::AntiClockwise => Self::Clockwise,
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

    /// Returns the index from `0..4` which this `Side` will appear in the order specified by a
    /// [`Box_`]
    fn index(self) -> usize {
        match self {
            Side::Left => 0,
            Side::Top => 1,
            Side::Right => 2,
            Side::Bottom => 3,
        }
    }
}
