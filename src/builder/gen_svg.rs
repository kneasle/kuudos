//! Code for converting a `Builder` into a debug-able SVG string

use crate::{
    types::{BoxIdx, VertVec},
    utils::{self, circle_passing_through},
    Side, V2Ext, V2,
};

use super::{Builder, EdgeLinkStyle};

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

    // This bounding box is in **untransformed** space
    let (bbox_min, bbox_max) =
        utils::bbox(bdr.verts.iter().map(|v| v.position)).unwrap_or_else(|| (V2::ZERO, V2::ZERO));
    let padding_vec = V2::ONE * PADDING;

    // Transform the vertices so that their min-point is (padding, padding) and they're scaled up
    // by the `scaling` factor
    let transform = |pt: V2| (pt - bbox_min + padding_vec) * scaling;
    let transformed_verts = bdr.verts.map(|v| transform(v.position));
    // Compute the dimensions of the resulting SVG image
    let img_dimensions = (bbox_max - bbox_min + padding_vec * 2.0) * scaling;

    // Root SVG element
    let mut root = XMLElement::new("svg");
    root.add_attribute("width", &img_dimensions.x.to_string());
    root.add_attribute("height", &img_dimensions.y.to_string());

    // Boxes
    let stroke_width_str = stroke_width.to_string();
    for box_ in bdr.boxes.iter() {
        let vert_coords = box_
            .vert_idxs
            .iter()
            .map(|idx| transformed_verts[*idx])
            .collect_vec();
        let coord_string = vert_coords
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
        for ((&v1, &v2), name) in vert_coords
            .iter()
            .circular_tuple_windows()
            .zip_eq(edge_names)
        {
            let d = v2 - v1; // A vector pointing down the line
            let normal = d.normal().normalise(); // A vector pointing into the box
            let centre = v1 + d / 2.0;
            let text_coords = centre + normal * stroke_width * TEXT_DIST_FROM_EDGE;

            let transform_str = format!(
                "translate({},{}) rotate({})",
                text_coords.x,
                text_coords.y,
                Rad((-d).angle()).to_deg().0
            );
            let text_size =
                std::cmp::min_by(d.length() * FONT_SIZE, scaling * FONT_SIZE_MAX, |a, b| {
                    a.partial_cmp(b).expect("Text size should not be NaN")
                });

            let mut text_elem = XMLElement::new("text");
            text_elem.add_attribute("transform", &transform_str);
            text_elem.add_attribute("font-size", &text_size.to_string());
            text_elem.add_attribute("font-family", "monospace");
            text_elem.add_attribute("text-anchor", "middle");
            text_elem.add_text(name);
            root.add_child(text_elem);
        }
    }

    // Edge links
    for link in bdr.edge_links.iter() {
        // Get the vertices on either side of the link
        let (start_pos, start_norm) =
            edge_midpoint_and_normal(bdr, &transformed_verts, link.start_box_idx, link.start_side);
        let (end_pos, _end_norm) =
            edge_midpoint_and_normal(bdr, &transformed_verts, link.end_box_idx, link.end_side);

        // Compute the shape of the ellipse, or default to lines
        let circle = match link.style {
            EdgeLinkStyle::Arc => circle_passing_through(start_pos, start_norm, end_pos),
            EdgeLinkStyle::Linear => None,
        };
        if let Some((centre, radius, start_angle, end_angle)) = circle {
            // Draw a circular arc between the points
            let mut line_elem = XMLElement::new("path");
            line_elem.add_attribute(
                "d",
                &utils::svg_circle_arc_path_str(centre, radius, start_angle, end_angle),
            );
            line_elem.add_attribute("fill", "none");
            line_elem.add_attribute("stroke", "black");
            line_elem.add_attribute("stroke-linecap", "round");
            line_elem.add_attribute("stroke-width", &stroke_width_str);
            root.add_child(line_elem);
        } else {
            // Add a line between the two edge centres
            let mut line_elem = XMLElement::new("line");
            line_elem.add_attribute("x1", &start_pos.x.to_string());
            line_elem.add_attribute("y1", &start_pos.y.to_string());
            line_elem.add_attribute("x2", &end_pos.x.to_string());
            line_elem.add_attribute("y2", &end_pos.y.to_string());
            line_elem.add_attribute("stroke", "black");
            line_elem.add_attribute("stroke-linecap", "round");
            line_elem.add_attribute("stroke-width", &stroke_width_str);
            root.add_child(line_elem);
        }
    }

    root.to_string()
}

fn edge_midpoint_and_normal(
    bdr: &Builder,
    vert_positions: &VertVec<V2>,
    box_idx: BoxIdx,
    side: Side,
) -> (V2, V2) {
    // Extract vertex positions
    let (vert1, vert2) = bdr.boxes[box_idx].get_edge_verts(side);
    let pos1 = vert_positions[vert1];
    let pos2 = vert_positions[vert2];
    // Compute values
    let midpoint = V2::lerp(pos1, pos2, 0.5);
    let normal = (pos1 - pos2).normal().normalise(); // Vector facing away from the box
    (midpoint, normal)
}
