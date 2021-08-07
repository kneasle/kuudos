//! Image specification and rendering utilities.  In essence, this is an intermediate
//! representation (IR) for vector images, much like LLVM's IR: we have multiple 'front ends' (e.g.
//! [`Builder`] and [`Shape`]) which output IR that can be translated to any output format (SVG,
//! HTML canvas, etc.).

mod ir;
mod lowering;
pub mod svg;

pub use ir::{
    ConcreteFillStyle, ConcreteStrokeStyle, ConcreteTextStyle, FillStyle, StrokeStyle, TextAnchor,
    TextStyle,
};
pub use lowering::RenderingOpts;

/// Re-export of [`ir::Image`] with the type params needed by the rest of the code
pub type Image = ir::Image<FillStyle, StrokeStyle, TextStyle>;
/// Version of [`ir::Image`] where all the styles are fully specified
pub type LoweredImage = ir::Image<ConcreteFillStyle, ConcreteStrokeStyle, ConcreteTextStyle>;

/// Re-export of [`ir::Elem`] with the type params needed by the rest of the code
pub type Elem = ir::Elem<FillStyle, StrokeStyle, TextStyle>;
/// Version of [`ir::Elem`] where all the styles are fully specified
pub type LoweredElem = ir::Elem<ConcreteFillStyle, ConcreteStrokeStyle, ConcreteTextStyle>;

/// Re-export of [`ir::Style`] with the type params needed by the rest of the code
pub type Style = ir::Style<FillStyle, StrokeStyle>;
/// Version of [`ir::Style`] where all the styles are fully specified
pub type LoweredStyle = ir::Style<ConcreteFillStyle, ConcreteStrokeStyle>;
