//! Code to convert an [`Image`] to an SVG string

use angle::Angle;
use itertools::Itertools;
use simple_xml_builder::XMLElement;

use crate::{image::TextAnchor, utils::Rect2, V2Ext, V2};

use super::{
    lowering::lower, ConcreteFillStyle, ConcreteStrokeStyle, ConcreteTextStyle, Image, LoweredElem,
    LoweredImage, RenderingOpts,
};

/// Generate an SVG root element for an unlowered [`Image`]
pub fn gen_svg(image: &Image, scale: f32, opts: &RenderingOpts) -> XMLElement {
    // Lower image then delegate to `gen_svg_from_lowered`
    gen_svg_from_lowered(&lower(image, opts), opts.margin, scale)
}

/// Generate an SVG root element for a lowered [`Image`] (i.e. one where all styles are stated
/// concretely).
pub fn gen_svg_from_lowered(lowered_image: &LoweredImage, margin: f32, scale: f32) -> XMLElement {
    let bbox = lowered_image
        .bbox()
        // An empty image has a 0-size bbox at the origin
        .unwrap_or_else(|| Rect2::point(V2::ZERO));
    let margin_vec = V2::ONE * margin;
    let translation = margin_vec - bbox.min(); // Maps the `min` corner of the bbox to inside the margin

    // Generate the root SVG element
    let mut root = XMLElement::new("svg");
    let image_dimensions = (bbox.max() - bbox.min() + margin_vec * 2.0) * scale;
    root.add_attribute("width", &image_dimensions.x.to_string());
    root.add_attribute("height", &image_dimensions.y.to_string());
    root.add_attribute("xmlns", "http://www.w3.org/2000/svg");
    // Translate all `Elem`s to SVG's `XMLElement`s
    for e in lowered_image.elements() {
        root.add_child(gen_svg_elem(e, translation, scale));
    }

    root
}

/// Creates an [`XMLElement`] for a given [`LoweredElem`]
fn gen_svg_elem(elem: &LoweredElem, translation: V2, scale: f32) -> XMLElement {
    // Helper closure to transform (i.e. translate then scale) a point
    let transform_point = |pt: &V2| (*pt + translation) * scale;

    // Create an unstyled XML element
    let mut xml_elem = match elem {
        // For `LineSegment`, we use the SVG `<line>` element
        LoweredElem::LineSegment(untransformed_pt_1, untransformed_pt_2, _) => {
            // Transform the points
            let p1 = transform_point(untransformed_pt_1);
            let p2 = transform_point(untransformed_pt_2);
            // Create the element
            let mut elem = XMLElement::new("line");
            elem.add_attribute("x1", &p1.x.to_string());
            elem.add_attribute("y1", &p1.y.to_string());
            elem.add_attribute("x2", &p2.x.to_string());
            elem.add_attribute("y2", &p2.y.to_string());
            elem
        }
        // For `CircularArc`, we use an SVG `<path>` element
        LoweredElem::CircularArc(untransformed_arc, _) => {
            // Transform the arc
            let arc = untransformed_arc.transform(translation, scale);
            // Create the SVG element
            let mut elem = XMLElement::new("path");
            elem.add_attribute("d", &arc.svg_path_str());
            elem
        }
        // For `Text`, we use (drumroll, please) an SVG `<text>` element
        LoweredElem::Text {
            position,
            text,
            angle,
            style: ConcreteTextStyle {
                font_size,
                font_family,
                anchor,
                .. // Other styles will be filled later
            },
        } => {
            // Generate SVG values
            let anchor_str = match anchor {
                TextAnchor::Left => "left",
                TextAnchor::Middle => "middle",
                TextAnchor::Right => "right",
            };
            // Compute the transform (i.e. translation then rotation) of the text element
            let transformed_pos = transform_point(position);
            let transform_str = format!(
                "translate({},{}) rotate({})",
                transformed_pos.x,
                transformed_pos.y,
                angle.to_deg().0
            );
            // Create the text element
            let mut elem = XMLElement::new("text");
            elem.add_attribute("transform", &transform_str);
            elem.add_attribute("font-size", &(font_size * scale).to_string());
            elem.add_attribute("font-family", &font_family);
            elem.add_attribute("text-anchor", anchor_str);
            elem.add_text(text);
            elem
        }
        // For `Polygon`, we use an SVG `<polygon>` element
        LoweredElem::Polygon(vertices, _) => {
            // Transform the vertices & generate the SVG point string.  This string is a
            // whitespace-delimited list of vertices, expressed as `x,y` pairs.  So the unit square
            // at the origin would have `coord_string = "0,0 0,1 1,1 1,0"`
            let coord_string = vertices
                .iter()
                .map(|vert| {
                    let vert = transform_point(vert);
                    format!("{},{}", vert.x, vert.y)
                })
                .join(" ");
            // Create the SVG element
            let mut cell_elem = XMLElement::new("polygon");
            cell_elem.add_attribute("points", &coord_string);
            cell_elem
        }
    };
    // Add styles
    add_fill_style_attrs(elem.fill_style(), &mut xml_elem);
    add_stroke_style_attrs(elem.stroke_style(), scale, &mut xml_elem);
    // Return the new element
    xml_elem
}

/// Add SVG attributes to give an [`XMLElement`] a given `FillStyle`
fn add_fill_style_attrs(style: Option<&ConcreteFillStyle>, xml_elem: &mut XMLElement) {
    match style {
        Some(s) => xml_elem.add_attribute("fill", &s.fill_color.to_string()),
        None => xml_elem.add_attribute("fill", "none"), // We need `fill="none"` to disable the fill
    }
}

/// Add SVG attributes to give an [`XMLElement`] a given `StrokeStyle`
fn add_stroke_style_attrs(
    style: Option<&ConcreteStrokeStyle>,
    scale: f32,
    xml_elem: &mut XMLElement,
) {
    match style {
        Some(s) => {
            xml_elem.add_attribute("stroke", &s.stroke_color.to_string());
            xml_elem.add_attribute("stroke-width", &(s.line_width * scale).to_string());
            xml_elem.add_attribute("stroke-linecap", "round"); // Always round off the line ends
            xml_elem.add_attribute("stroke-linejoin", "round"); // Always round off internal corners
        }
        None => xml_elem.add_attribute("stroke", "none"), // Put `stroke="none"` if no stroke
    }
}
