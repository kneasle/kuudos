//! An intermediate representation for vector images

use angle::Rad;
use rgb::RGB8;

use crate::{
    utils::{CircularArc, Rect2},
    V2,
};

use super::{lowering, svg, LoweredImage, RenderingOpts};

/// A full [`Image`], composed of many [`Elem`]ents.  [`Image`]s use the same units as [`Shape`] -
/// i.e. the side length of a cell is roughly 1 unit.
#[derive(Debug, Clone)]
pub struct Image<F, S, T> {
    pub(super) elements: Vec<Elem<F, S, T>>,
}

impl<F, S, T> Image<F, S, T> {
    /// Creates an empty `Image` (i.e. one which contains no [`Elem`]s)
    pub fn empty() -> Self {
        Self { elements: vec![] }
    }

    /// Adds a new [`Elem`] to this `Image`
    pub fn add(&mut self, elem: Elem<F, S, T>) {
        self.elements.push(elem)
    }

    /// Compute the smallest [`Rect2`] which fits around every [`Elem`] in this `Image`.  This
    /// returns `None` if the `Image` contains no [`Elem`]s.
    pub fn bbox(&self) -> Option<Rect2> {
        Rect2::union_iter(self.elements.iter().map(Elem::bbox))
    }

    pub fn elements(&self) -> &[Elem<F, S, T>] {
        self.elements.as_slice()
    }
}

/// Helper methods for easy conversions to various formats
impl Image<FillStyle, StrokeStyle, TextStyle> {
    pub fn svg_string(&self, scale: f32, opts: &RenderingOpts) -> String {
        svg::gen_svg(self, scale, opts).to_string()
    }

    pub fn lower(&self, opts: &RenderingOpts) -> LoweredImage {
        lowering::lower(self, opts)
    }
}

/// The shape of an [`Elem`]
#[derive(Debug, Clone)]
pub enum Elem<F, S, T> {
    LineSegment(V2, V2, S),
    CircularArc(CircularArc, Style<F, S>),
    Text {
        position: V2,
        text: String,
        angle: Rad<f32>,
        style: T,
    },
    Polygon(Vec<V2>, Style<F, S>),
}

impl<F, S, T> Elem<F, S, T> {
    pub fn arc_passing_through(p1: V2, tangent_1: V2, p2: V2, style: Style<F, S>) -> Option<Self> {
        CircularArc::arc_passing_through(p1, tangent_1, p2).map(|arc| Elem::CircularArc(arc, style))
    }

    /// Returns the smallest [`Rect2`] which contains this element
    pub fn bbox(&self) -> Rect2 {
        match self {
            Elem::LineSegment(p1, p2, _) => Rect2::bbox([*p1, *p2]).unwrap(),
            Elem::CircularArc(arc, _) => arc.bbox(),
            Elem::Polygon(verts, _) => Rect2::bbox(verts.iter().copied()).unwrap(),
            // Naively assume that the text elements take up no space
            Elem::Text { position, .. } => Rect2::point(*position),
        }
    }
}

impl<S> Elem<ConcreteFillStyle, S, ConcreteTextStyle> {
    /// Gets the fill style of this `Elem`, if it exists.  It may not exist - for example,
    /// [`Elem::LineSegment`]s can't be filled.
    pub fn fill_style(&self) -> Option<&ConcreteFillStyle> {
        match self {
            Elem::LineSegment(_, _, _) => None, // Line segments can't be filled
            Elem::Text { style, .. } => Some(&style.fill_style),
            Elem::CircularArc(_, style) => style.fill_style(),
            Elem::Polygon(_, style) => style.fill_style(),
        }
    }

    /// Gets the stroke style of this `Elem`, if it exists.
    pub fn stroke_style(&self) -> Option<&S> {
        match self {
            Elem::Text { .. } => None, // Text elements can't be stroked
            Elem::LineSegment(_, _, stroke_style) => Some(stroke_style),
            Elem::CircularArc(_, style) => style.stroke_style(),
            Elem::Polygon(_, style) => style.stroke_style(),
        }
    }
}

/////////////
// STYLING //
/////////////

/// The full styling of an [`Elem`], which is either filled or stroked or both (but invisible
/// elements are not possible).
#[derive(Debug, Clone)]
pub enum Style<F, S> {
    JustFill(F),
    JustStroke(S),
    FillAndStroke(F, S),
}

impl<F, S> Style<F, S> {
    pub fn fill_style(&self) -> Option<&F> {
        match self {
            Self::JustFill(f) => Some(f),
            Self::JustStroke(_) => None,
            Self::FillAndStroke(f, _) => Some(f),
        }
    }

    pub fn stroke_style(&self) -> Option<&S> {
        match self {
            Self::JustFill(_) => None,
            Self::JustStroke(s) => Some(s),
            Self::FillAndStroke(_, s) => Some(s),
        }
    }
}

/// The visual style of the body of an [`Elem`]
#[derive(Debug, Clone)]
pub enum FillStyle {
    /// The fill colour of the cells' background
    Background,
    /// The fill colour of any foreground (i.e. cells' contents)
    Foreground,
    Custom(ConcreteFillStyle),
}

/// A fully specified [`FillStyle`]
#[derive(Debug, Clone)]
pub struct ConcreteFillStyle {
    pub fill_color: RGB8,
}

/// The visual style of the outline of an [`Elem`]
#[derive(Debug, Clone)]
pub enum StrokeStyle {
    /// The border between two cells that are in different boxes, or 'external' edges
    BoxBorder,
    /// The border between two cells which share a 'lane' (row or column), but **not** a box
    CellBorder,
    /// The border between two cells which are adjacent but disconnected (this shouldn't happen,
    /// and is more of an 'error' case than anything else)
    Disconnected,
    Custom(ConcreteStrokeStyle),
}

/// A fully specified [`StrokeStyle`]
#[derive(Debug, Clone)]
pub struct ConcreteStrokeStyle {
    pub line_width: f32,
    pub stroke_color: RGB8,
}

/// The visual style of some [`Elem::Text`]
#[derive(Debug, Clone)]
pub enum TextStyle {
    /// The text is a clue provided by the puzzle
    Clue,
    /// The text is a digit penned by the user
    PennedDigit,
    Custom(ConcreteTextStyle),
}

/// A fully specified [`TextStyle`]
#[derive(Debug, Clone)]
pub struct ConcreteTextStyle {
    pub fill_style: ConcreteFillStyle,
    pub font_size: f32,
    pub font_family: String,
    pub anchor: TextAnchor,
}

/// The horizontal location of the text, relative to its anchor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TextAnchor {
    Left,
    Middle,
    Right,
}
