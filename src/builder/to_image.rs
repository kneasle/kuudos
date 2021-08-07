//! Code for converting a `Builder` into a debug-able SVG string

use crate::{
    image::{
        ConcreteFillStyle, ConcreteStrokeStyle, ConcreteTextStyle, LoweredElem, LoweredImage,
        LoweredStyle, TextAnchor,
    },
    utils::CircularArc,
    V2,
};

use super::{Builder, EdgeLinkStyle, RotateDirection};

use angle::Rad;
use itertools::Itertools;
use rgb::RGB8;

const LINE_WIDTH: f32 = 0.2; // (cell) edge lengths
const TEXT_DIST_FROM_EDGE: f32 = 1.0; // multiples of the line width
const FONT_SIZE: f32 = 0.2; // multiple of the length of the box's edge
const FONT_SIZE_MAX: f32 = 0.5; // (cell) edge lengths

const BACKGROUND_COLOR: RGB8 = RGB8::new(255, 255, 255); // Box backgrounds are white
const FOREGROUND_COLOR: RGB8 = RGB8::new(0, 0, 0); // Box borders are solid black
const HIDDEN_COLOR: RGB8 = RGB8::new(170, 170, 170); // Hidden items are rendered as light grey
const ERROR_COLOR: RGB8 = RGB8::new(255, 0, 0); // Errors are rendered as red

/// Generate a simple SVG string representing the current state of a [`Builder`].  Unlike
/// generating an SVG string by converting through a [`Shape`], this conversion will always succeed
/// even if the [`Builder`] cannot be converted to a [`Shape`].
pub fn gen_image(bdr: &Builder) -> LoweredImage {
    let mut image = LoweredImage::empty();

    for (_idx, source_box) in bdr.source_boxes() {
        // Get vertex coordinates
        let vert_coords = source_box
            .vert_idxs
            .iter()
            .map(|idx| bdr.verts[*idx])
            .collect_vec();
        // Generate polygon element for the box background & border
        let poly_elem = LoweredElem::Polygon(
            vert_coords.clone(),
            LoweredStyle::FillAndStroke(
                ConcreteFillStyle {
                    fill_color: BACKGROUND_COLOR,
                },
                ConcreteStrokeStyle {
                    line_width: LINE_WIDTH,
                    stroke_color: FOREGROUND_COLOR,
                },
            ),
        );
        image.add(poly_elem);

        // Add direction labels on all edges
        let edge_names = ["LEFT", "TOP", "RIGHT", "BOTTOM"];
        for ((&bottom_vert, &top_vert), name) in vert_coords
            .iter()
            .circular_tuple_windows()
            .zip_eq(edge_names)
        {
            // Do the vector maths to figure out where the text should go
            let (v1, v2) = match source_box.rotate_direction {
                RotateDirection::Clockwise => (bottom_vert, top_vert),
                RotateDirection::AntiClockwise => (top_vert, bottom_vert), // Flip verts on
                                                                           // anti-clockwise faces
            };
            let d = v2 - v1; // A vector pointing down the line
            let normal = d.normal().normalise(); // A vector pointing into the box
            let centre = v1 + d / 2.0;
            let text_coords = centre + normal * LINE_WIDTH * TEXT_DIST_FROM_EDGE;
            let text_angle = Rad((-d).angle());
            // Compute the text size (made harder because NaNs are stupid and can't be ordered)
            let font_size = std::cmp::min_by(d.length() * FONT_SIZE, FONT_SIZE_MAX, |a, b| {
                a.partial_cmp(b).expect("Text size should not be NaN")
            });

            // Create the text element
            let text_style = ConcreteTextStyle {
                fill_style: ConcreteFillStyle {
                    fill_color: FOREGROUND_COLOR,
                },
                font_size,
                font_family: "monospace".to_owned(),
                anchor: TextAnchor::Middle,
            };
            let text_elem = LoweredElem::Text {
                position: text_coords,
                angle: text_angle,
                text: name.to_owned(),
                style: text_style,
            };
            // Add text element to image
            image.add(text_elem);
        }
    }

    enum LineType {
        Full,   // The link is meant to be a line
        Error,  // The link is being rendered as a line as a fallback for some other shape
        Hidden, // The link is meant to be hidden, but is visible for debugging
    }

    enum LinkShape {
        CircularArc(CircularArc), // Draw a circular arc (angles are taken clockwise from the +X axis)
        Line(V2, V2),             // Draw a straight line segment between two points.
    }

    // Edge links
    for link in bdr.edge_links.iter() {
        // Get the vertex positions
        let start_vert_left = bdr.verts[link.vert_idx_top_left];
        let start_vert_right = bdr.verts[link.vert_idx_top_right];
        let end_vert_left = bdr.verts[link.vert_idx_bottom_left];
        let end_vert_right = bdr.verts[link.vert_idx_bottom_right];
        // Get the midpoints and normals of the two ends of the link (with the normals pointing
        // away from the line)
        let (start_pos, start_norm) = midpoint_and_normal(start_vert_right, start_vert_left);
        let (end_pos, _end_norm) = midpoint_and_normal(end_vert_left, end_vert_right);

        // Decide how to display this link
        let (element, line_type) = match link.style {
            EdgeLinkStyle::Arc => {
                match CircularArc::arc_passing_through(start_pos, start_norm, end_pos) {
                    // Draw the arc if it is well defined
                    Some(arc) => (LinkShape::CircularArc(arc), LineType::Full),
                    // If the arc doesn't exist, fall back on an error line
                    None => (LinkShape::Line(start_pos, end_pos), LineType::Error),
                }
            }
            EdgeLinkStyle::Linear => (LinkShape::Line(start_pos, end_pos), LineType::Full),
            EdgeLinkStyle::Hidden => (LinkShape::Line(start_pos, end_pos), LineType::Hidden),
        };

        // Generate the SVG element for this link
        let stroke_style = ConcreteStrokeStyle {
            line_width: LINE_WIDTH,
            stroke_color: match line_type {
                LineType::Full => FOREGROUND_COLOR,
                LineType::Hidden => HIDDEN_COLOR,
                LineType::Error => ERROR_COLOR,
            },
        };
        let link_elem = match element {
            LinkShape::CircularArc(arc) => {
                LoweredElem::CircularArc(arc, LoweredStyle::JustStroke(stroke_style))
            }
            LinkShape::Line(p1, p2) => LoweredElem::LineSegment(p1, p2, stroke_style),
        };
        // Add link element to the image
        image.add(link_elem);
    }

    image
}

fn midpoint_and_normal(pos1: V2, pos2: V2) -> (V2, V2) {
    // Compute values
    let midpoint = V2::lerp(pos1, pos2, 0.5);
    let normal = (pos1 - pos2).normal().normalise(); // Vector facing away from the box
    (midpoint, normal)
}
