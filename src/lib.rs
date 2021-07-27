use vector2d::Vector2D;

mod builder;
mod shape;
pub mod solve;
pub mod svg;
mod types;

pub use builder::{Builder, Direction, Side};
pub use shape::{Shape, Symmetry};

/// Type alias for 2D floating point vectors (in the geometric sense, unlike [`Vec`])
pub type V2 = Vector2D<f32>;
