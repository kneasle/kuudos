//! Code for converting a `Builder` into a debug-able SVG string

use crate::{
    shape::LinkShape,
    utils::{CircularArc, Rect2},
    V2Ext, V2,
};

use super::{Builder, EdgeLinkStyle, RotateDirection};

use angle::{Angle, Rad};
use itertools::Itertools;
use simple_xml_builder::XMLElement;

const PADDING: f32 = 0.5; // (cell) edge lengths
const LINE_WIDTH: f32 = 0.2; // (cell) edge lengths
const TEXT_DIST_FROM_EDGE: f32 = 1.0; // multiples of the line width
const FONT_SIZE: f32 = 0.2; // multiple of the length of the box's edge
const FONT_SIZE_MAX: f32 = 0.5; // (cell) edge lengths

/// Generate a simple SVG string representing the current state of a [`Builder`].  Unlike
/// generating an SVG string by converting through a [`Shape`], this conversion will always succeed
/// even if the [`Builder`] cannot be converted to a [`Shape`].
pub fn gen_svg(bdr: &Builder, scaling: f32) -> String {
    let stroke_width = LINE_WIDTH * scaling;

    // Transform the vertices so that their minimum point is (padding, padding) and they're scaled
    // up by the `scaling` factor
    let bbox = Rect2::bbox(bdr.verts.iter().copied()) // This bbox is in **untansformed** space
        .unwrap_or_else(|| Rect2::from_min_size(V2::ZERO, V2::ZERO)); // Default to the empty rect
    let padding_vec = V2::ONE * PADDING;
    let transform = |pt: V2| (pt - bbox.min() + padding_vec) * scaling;
    let transformed_verts = bdr.verts.map(|v| transform(*v));
    // Compute the dimensions of the resulting SVG image
    let img_dimensions = (bbox.max() - bbox.min() + padding_vec * 2.0) * scaling;

    // Root SVG element
    let mut root = XMLElement::new("svg");
    root.add_attribute("width", &img_dimensions.x.to_string());
    root.add_attribute("height", &img_dimensions.y.to_string());

    // Boxes
    let stroke_width_str = stroke_width.to_string();
    for (_idx, source_box) in bdr.source_boxes() {
        let transformed_vert_coords = source_box
            .vert_idxs
            .iter()
            .map(|idx| transformed_verts[*idx])
            .collect_vec();

        // Create the polygon SVG element
        let coord_string = transformed_vert_coords
            .iter()
            .map(|vert| format!("{},{}", vert.x, vert.y))
            .join(" ");
        let mut border_elem = XMLElement::new("polygon");
        border_elem.add_attribute("points", &coord_string);
        border_elem.add_attribute("fill", "white");
        border_elem.add_attribute("stroke", "black");
        border_elem.add_attribute("stroke-linejoin", "round");
        border_elem.add_attribute("stroke-width", &stroke_width_str);
        root.add_child(border_elem);

        // Label each edge of the box with its direction
        let edge_names = ["LEFT", "TOP", "RIGHT", "BOTTOM"];
        for ((&bottom_vert, &top_vert), name) in transformed_vert_coords
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
            let text_coords = centre + normal * stroke_width * TEXT_DIST_FROM_EDGE;
            let text_angle_deg = Rad((-d).angle()).to_deg().0;
            // Compute the text size (made harder because NaNs are stupid and can't be ordered)
            let text_size =
                std::cmp::min_by(d.length() * FONT_SIZE, scaling * FONT_SIZE_MAX, |a, b| {
                    a.partial_cmp(b).expect("Text size should not be NaN")
                });
            // Create the SVG text element
            let transform_str = format!(
                "translate({},{}) rotate({})",
                text_coords.x, text_coords.y, text_angle_deg
            );
            let mut text_elem = XMLElement::new("text");
            text_elem.add_attribute("transform", &transform_str);
            text_elem.add_attribute("font-size", &text_size.to_string());
            text_elem.add_attribute("font-family", "monospace");
            text_elem.add_attribute("text-anchor", "middle");
            text_elem.add_text(name);
            root.add_child(text_elem);
        }
    }

    enum EdgeLinkLineType {
        Full,   // The link is meant to be a line
        Error,  // The link is being rendered as a line as a fallback for some other shape
        Hidden, // The link is meant to be hidden, but is visible for debugging
    }

    // Edge links
    for link in bdr.edge_links.iter() {
        // Get the vertex positions
        let start_vert_left = transformed_verts[link.vert_idx_top_left];
        let start_vert_right = transformed_verts[link.vert_idx_top_right];
        let end_vert_left = transformed_verts[link.vert_idx_bottom_left];
        let end_vert_right = transformed_verts[link.vert_idx_bottom_right];

        // Get the midpoints and normals of the two ends of the link (with the normals pointing
        // away from the line)
        let (start_pos, start_norm) = midpoint_and_normal(start_vert_right, start_vert_left);
        let (end_pos, _end_norm) = midpoint_and_normal(end_vert_left, end_vert_right);

        // Compute how this edge link should be displayed
        let (element, line_type) = match link.style {
            EdgeLinkStyle::Arc => LinkShape::arc_passing_through(start_pos, start_norm, end_pos)
                .map(|arc| (arc, EdgeLinkLineType::Full)) // Draw the arc if it exists
                .unwrap_or_else(|| {
                    // If the arc doesn't exist, fall back on an error line
                    (LinkShape::Line(start_pos, end_pos), EdgeLinkLineType::Error)
                }),
            EdgeLinkStyle::Linear => (LinkShape::Line(start_pos, end_pos), EdgeLinkLineType::Full),
            EdgeLinkStyle::Hidden => (
                LinkShape::Line(start_pos, end_pos),
                EdgeLinkLineType::Hidden,
            ),
        };

        // Generate the SVG element for this link
        let stroke_color = match line_type {
            EdgeLinkLineType::Full => "black",
            EdgeLinkLineType::Hidden => "rgb(170,170,170)",
            EdgeLinkLineType::Error => "red",
        };
        match element {
            LinkShape::CircularArc(CircularArc {
                centre,
                radius,
                start_angle,
                end_angle,
            }) => {
                // Draw a circular arc between the points
                let mut line_elem = XMLElement::new("path");
                line_elem.add_attribute(
                    "d",
                    &CircularArc {
                        centre,
                        radius,
                        start_angle,
                        end_angle,
                    }
                    .svg_path_str(),
                );
                line_elem.add_attribute("fill", "none");
                line_elem.add_attribute("stroke", stroke_color);
                line_elem.add_attribute("stroke-linecap", "round");
                line_elem.add_attribute("stroke-width", &stroke_width_str);
                root.add_child(line_elem);
            }
            LinkShape::Line(p1, p2) => {
                // Add a line between the two edge centres
                let mut line_elem = XMLElement::new("line");
                line_elem.add_attribute("x1", &p1.x.to_string());
                line_elem.add_attribute("y1", &p1.y.to_string());
                line_elem.add_attribute("x2", &p2.x.to_string());
                line_elem.add_attribute("y2", &p2.y.to_string());
                line_elem.add_attribute("stroke", stroke_color);
                line_elem.add_attribute("stroke-linecap", "round");
                line_elem.add_attribute("stroke-width", &stroke_width_str);
                root.add_child(line_elem);
            }
        }
    }

    root.to_string()
}

fn midpoint_and_normal(pos1: V2, pos2: V2) -> (V2, V2) {
    // Compute values
    let midpoint = V2::lerp(pos1, pos2, 0.5);
    let normal = (pos1 - pos2).normal().normalise(); // Vector facing away from the box
    (midpoint, normal)
}
