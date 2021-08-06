use itertools::Itertools;
use rgb::RGB8;

use crate::{V2Ext, V2};

use super::{
    ir::{ConcreteTextStyle, TextStyle},
    ConcreteFillStyle, ConcreteStrokeStyle, Elem, FillStyle, Image, LoweredElem, LoweredImage,
    LoweredStyle, StrokeStyle, Style, TextAnchor,
};

/// 'Lower' an [`Image`] to a [`LoweredImage`] - i.e. use values from [`RenderingOpts`] to make
/// concrete versions for every style.  The [`LoweredImage`] can then be unambiguously converted to
/// an image file, and the [`RenderingOpts`] is no longer needed.
pub(super) fn lower(image: &Image, opts: &RenderingOpts) -> LoweredImage {
    // Lower each element in turn
    LoweredImage {
        elements: image
            .elements
            .iter()
            .map(|e| lower_elem(e, opts))
            .collect_vec(),
    }
}

/// 'Lower' a single [`Elem`] to a [`LoweredElem`], by converting all `XStyle`s to
/// `ConcreteXStyle`s.  This function just lowers the element's styles one by one
fn lower_elem(elem: &Elem, opts: &RenderingOpts) -> LoweredElem {
    match elem {
        Elem::LineSegment(pt1, pt2, stroke_style) => {
            LoweredElem::LineSegment(*pt1, *pt2, lower_stroke_style(stroke_style, opts))
        }
        Elem::CircularArc(arc, style) => LoweredElem::CircularArc(*arc, lower_style(style, opts)),
        Elem::Polygon(pts, style) => LoweredElem::Polygon(pts.to_owned(), lower_style(style, opts)),
        Elem::Text {
            position,
            text,
            angle,
            style,
        } => {
            let lowered_style = lower_text_style(style, opts); // Lower the styling
            let text_nudge =
                (V2::DOWN * opts.text_vertical_nudge * lowered_style.font_size).rotate(*angle);
            LoweredElem::Text {
                // AFAIK there's no way to accurately render text by its centre point.  Therefore,
                // the next best thing is to move all the text 'down' by some proportion of its
                // font size (the exact value of `text_vertical_nudge` should be found per-font by
                // experimentation).
                position: *position + text_nudge,
                text: text.to_owned(),
                angle: *angle,
                style: lowered_style,
            }
        }
    }
}

/// Lowers a [`Style`] by replacing any non-concrete styles (like
/// [`StrokeStyle::CellBorder`]) with a [`ConcreteStrokeStyle`] who's values are populated from
/// `opts`.
fn lower_style(style: &Style, opts: &RenderingOpts) -> LoweredStyle {
    match style {
        Style::JustFill(f) => LoweredStyle::JustFill(lower_fill_style(f, opts)),
        Style::JustStroke(s) => LoweredStyle::JustStroke(lower_stroke_style(s, opts)),
        Style::FillAndStroke(f, s) => {
            LoweredStyle::FillAndStroke(lower_fill_style(f, opts), lower_stroke_style(s, opts))
        }
    }
}

/// Lowers a [`FillStyle`] by replacing any non-concrete styles (like [`FillStyle::Background`])
/// with a [`ConcreteFillStyle`] who's values are populated from `opts`.
fn lower_fill_style(style: &FillStyle, opts: &RenderingOpts) -> ConcreteFillStyle {
    let fill_color = match style {
        FillStyle::Background => opts.cell_fill_color,
        FillStyle::Foreground => opts.clue_color,
        FillStyle::Custom(concrete_style) => return concrete_style.to_owned(), // No lowering needed
    };
    ConcreteFillStyle { fill_color }
}

/// Lowers a [`StrokeStyle`] by replacing any non-concrete styles (like
/// [`StrokeStyle::CellBorder`]) with a [`ConcreteStrokeStyle`] who's values are populated from
/// `opts`.
fn lower_stroke_style(style: &StrokeStyle, opts: &RenderingOpts) -> ConcreteStrokeStyle {
    let (width, stroke_color) = match style {
        StrokeStyle::BoxBorder => (opts.box_border_width, opts.box_border_color),
        StrokeStyle::CellBorder => (opts.cell_border_width, opts.cell_border_color),
        StrokeStyle::Disconnected => (opts.cell_border_width, opts.disconnected_edge_color),
        StrokeStyle::Custom(concrete_style) => return concrete_style.to_owned(), // No lowering needed
    };
    ConcreteStrokeStyle {
        line_width: width,
        stroke_color,
    }
}

/// Lowers a [`TextStyle`] by replacing any non-concrete styles (like [`TextStyle::Clue`]) with a
/// [`ConcreteTextStyle`] who's values are populated from `opts`.
fn lower_text_style(style: &TextStyle, opts: &RenderingOpts) -> ConcreteTextStyle {
    let (fill_color, font_size) = match style {
        TextStyle::Clue => (opts.clue_color, opts.clue_font_size),
        TextStyle::PennedDigit => (opts.pen_color, opts.pen_font_size),
        TextStyle::Custom(concrete_style) => return concrete_style.to_owned(), // No lowering needed
    };
    let text_style = ConcreteTextStyle {
        font_size,
        font_family: opts.font_family.to_owned(),
        anchor: TextAnchor::Middle,
        fill_style: ConcreteFillStyle { fill_color },
    };
    text_style
}

/// Configuration for how a sudoku should be rendered
#[derive(Debug, Clone)]
pub struct RenderingOpts {
    /// What color the cells should be filled
    cell_fill_color: RGB8,

    /// If the cell edges have length ~1 unit, how many units wide should cell borders be
    cell_border_width: f32,
    /// The colour of the cell borders.  Defaults to black
    cell_border_color: RGB8,

    /// If the cell edges have length ~1 unit, how many units wide should the borders of boxes be
    box_border_width: f32,
    /// The colour of the box borders.  Defaults to black
    box_border_color: RGB8,

    /// What font family should be used for populating the sudoku
    font_family: String,
    /// What multiple of `font_size` is added to the y-value of the text in order to centre the
    /// text
    text_vertical_nudge: f32,
    /// What font size should be used for the clues, if the cell edges have a length of ~1 unit
    clue_font_size: f32,
    /// What color the clues should be rendered
    clue_color: RGB8,
    /// What font size should be used for the numbers penned into cells, if the cell edges have a
    /// length of ~1 unit
    pen_font_size: f32,
    /// What color the clues should be rendered
    pen_color: RGB8,

    /// The colour of an edge between two cells which do not share any groups (this edge will be
    /// given the width of a cell border)
    disconnected_edge_color: RGB8,

    /// If the edge lengths are ~1 unit, how many units of space will be reserved round the edge of
    /// the SVG file
    pub(super) padding: f32,
}

impl RenderingOpts {
    pub fn font_size(&self) -> f32 {
        self.clue_font_size
    }

    pub fn text_vertical_nudge(&self) -> f32 {
        self.text_vertical_nudge
    }
}

impl Default for RenderingOpts {
    fn default() -> Self {
        Self {
            cell_fill_color: RGB8::new(255, 255, 255),

            cell_border_width: 0.05, // edge lengths
            cell_border_color: RGB8::new(1, 1, 1) * 200,

            box_border_width: 0.1, // edge lengths
            box_border_color: RGB8::new(0, 0, 0),

            font_family: "monospace".to_owned(),
            text_vertical_nudge: 0.33, // multiples of `font_size`
            clue_font_size: 0.7,       // edge lengths
            clue_color: RGB8::new(0, 0, 0),
            pen_font_size: 0.7, // edge lengths
            pen_color: RGB8::new(0, 200, 0),

            disconnected_edge_color: RGB8::new(255, 0, 0),

            padding: 0.5, // edge lengths
        }
    }
}
