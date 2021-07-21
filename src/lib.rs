use vector2d::Vector2D;

mod shape;
pub mod svg;

pub use shape::Shape;

/// Type alias for 2D floating point vectors (in the geometric sense, unlike [`Vec`])
type V2 = Vector2D<f32>;
