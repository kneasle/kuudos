use vector2d::Vector2D;

pub mod builder;
mod indexed_vec;
pub mod shape;
pub mod solve;
mod utils;

pub use builder::{Direction, Side};
pub use shape::{Shape, Symmetry};
pub use utils::{regular_polygon_inradius, V2Ext};

/// Type alias for 2D floating point vectors (in the geometric sense, unlike [`Vec`])
pub type V2 = Vector2D<f32>;
