pub mod builder;
pub mod image;
pub mod puzzle_gen;
pub mod shape;
pub mod solve;
mod utils;

pub use builder::{Direction, Side};
pub use shape::{Shape, Symmetry};
pub use utils::{regular_polygon_inradius, V2Ext};

/// Type alias for 2D floating point vectors (these are vectors in the geometric sense, not
/// extensible arrays like [`Vec`])
pub type V2 = vector2d::Vector2D<f32>;
