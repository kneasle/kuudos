use itertools::Itertools;
use rgb::RGB8;
use simple_xml_builder::XMLElement;

use super::{DrawElement, LinkShape, Shape};
use crate::{utils, V2Ext, V2};

/// Write a populated sudoku grid to an SVG string
pub fn gen_svg_string(
    shape: &Shape,
    opts: &RenderingOpts,
    scaling: f32,
    assignment: &[Option<char>],
) -> String {
    // This bounding box is in **untransformed** space
    let bbox = shape.bbox().expect("Shape can't have no vertices");
    let padding_vec = V2::ONE * opts.padding;

    let transform = |pt: &V2| (*pt - bbox.min() + padding_vec) * scaling;

    // Transform the vertices so that their min-point is (padding, padding) and they're scaled up
    // by the `scaling` factor
    let transformed_verts = shape.verts.map(transform);
    // Compute the dimensions of the resulting SVG image
    let img_dimensions = (bbox.max() - bbox.min() + padding_vec * 2.0) * scaling;

    // Root SVG element
    let mut root = XMLElement::new("svg");
    root.add_attribute("width", &img_dimensions.x.to_string());
    root.add_attribute("height", &img_dimensions.y.to_string());

    // Add extra elements behind the cells
    for elem in &shape.extra_elements {
        root.add_child(extra_elem_to_svg(elem, padding_vec - bbox.min(), scaling));
    }

    // Add the cell backgrounds
    for vert_idxs in shape.cells.iter() {
        let coord_string = vert_idxs
            .iter()
            .map(|&idx| {
                let vert = transformed_verts[idx];
                format!("{},{}", vert.x, vert.y)
            })
            .join(" ");

        let mut cell_elem = XMLElement::new("polygon");
        cell_elem.add_attribute("points", &coord_string);
        cell_elem.add_attribute("fill", &opts.cell_fill_color.to_string());
        root.add_child(cell_elem);
    }

    // Add the cell contents
    let text_size_str = (opts.font_size * scaling).to_string();
    for (vert_idxs, value) in shape.cells.iter().zip_eq(assignment) {
        if let Some(c) = *value {
            let verts = vert_idxs.iter().map(|idx| shape.verts[*idx]);
            let mut untransformed_centre =
                utils::centroid(verts).expect("Can't have cell with no versions");
            untransformed_centre.y += opts.font_size * opts.text_vertical_nudge;
            let centre = transform(&untransformed_centre);

            let mut content_elem = XMLElement::new("text");
            content_elem.add_attribute("x", &centre.x.to_string());
            content_elem.add_attribute("y", &centre.y.to_string());
            content_elem.add_attribute("font-size", &text_size_str);
            content_elem.add_attribute("font-family", &opts.font_family);
            content_elem.add_attribute("text-anchor", "middle");
            content_elem.add_text(&c.to_string());

            root.add_child(content_elem);
        }
    }

    #[derive(Debug, Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash)]
    enum EdgeType {
        /// The edge is between two disconnected cells
        Disconnected,
        CellBorder,
        BoxBorder,
    }

    // Sort the edges so that all the box borders come after all the cell borders (to prevent
    // overlap issues)
    let mut edges = shape
        .edges()
        .into_iter()
        .map(|edge| {
            let edge_type = match edge.cell_idx_left {
                Some(cell_left) => {
                    if shape.num_groups_shared_between(cell_left, edge.cell_idx_right) == 0 {
                        EdgeType::Disconnected // Adjacent cells which share no groups
                    } else if shape.do_cells_share_box(cell_left, edge.cell_idx_right) {
                        EdgeType::CellBorder // Adjacent cells which are in the same box
                    } else {
                        EdgeType::BoxBorder // Adjacent cells which share a group but aren't in the
                                            // same box
                    }
                }
                None => EdgeType::BoxBorder, // External edges
            };
            (edge_type, edge)
        })
        .collect_vec();
    // Sort the edges by type with `Disconnected` < `CellBorder` < `BoxBorder`
    edges.sort_by_key(|(ty, _edge)| *ty);

    // Add the edges
    let box_border_width_str = (opts.box_border_width * scaling).to_string();
    let cell_border_width_str = (opts.cell_border_width * scaling).to_string();
    for (edge_type, edge) in edges {
        // Extract the coordinates of the top/bottom vertices of this edge
        let vert_top = transformed_verts[edge.vert_idx_top];
        let vert_bottom = transformed_verts[edge.vert_idx_bottom];
        // Determine the edge's style, based on whether or not it's a box border
        let (color, width_str): (_, &str) = match edge_type {
            EdgeType::Disconnected => (opts.disconnected_edge_color, &cell_border_width_str),
            EdgeType::CellBorder => (opts.cell_border_color, &cell_border_width_str),
            EdgeType::BoxBorder => (opts.box_border_color, &box_border_width_str),
        };

        let mut edge_elem = XMLElement::new("line");
        // Coords
        edge_elem.add_attribute("x1", &vert_top.x.to_string());
        edge_elem.add_attribute("y1", &vert_top.y.to_string());
        edge_elem.add_attribute("x2", &vert_bottom.x.to_string());
        edge_elem.add_attribute("y2", &vert_bottom.y.to_string());
        // Stroke style
        edge_elem.add_attribute("stroke", &color.to_string());
        edge_elem.add_attribute("stroke-width", width_str);
        edge_elem.add_attribute("stroke-linecap", "round");
        root.add_child(edge_elem);
    }

    root.to_string()
}

/// Helper function to generate an [`XMLElement`] to represent a given [`DrawElement`]
fn extra_elem_to_svg(elem: &DrawElement, translation: V2, scaling: f32) -> XMLElement {
    match elem {
        DrawElement::EdgeLink(shape) => {
            let mut xml_elem = match shape.transform(translation, scaling) {
                LinkShape::Line(p1, p2) => {
                    let mut elem = XMLElement::new("line");
                    elem.add_attribute("x1", &p1.x.to_string());
                    elem.add_attribute("y1", &p1.y.to_string());
                    elem.add_attribute("x2", &p2.x.to_string());
                    elem.add_attribute("y2", &p2.y.to_string());
                    elem
                }
                LinkShape::CircularArc {
                    centre,
                    radius,
                    start_angle,
                    end_angle,
                } => {
                    let mut elem = XMLElement::new("path");
                    elem.add_attribute(
                        "d",
                        &utils::svg_circle_arc_path_str(centre, radius, start_angle, end_angle),
                    );
                    elem
                }
            };
            xml_elem.add_attribute("fill", "none");
            xml_elem.add_attribute("stroke", "black");
            xml_elem.add_attribute("stroke-width", "4");
            xml_elem.add_attribute("stroke-linecap", "round");
            xml_elem
        }
    }
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
    /// What font size should be used for the numbers written in the cells, if the cell edges have
    /// a length of ~1 unit
    font_size: f32,
    /// What multiple of `font_size` is added to the y-value of the text in order to centre the
    /// text
    text_vertical_nudge: f32,

    /// The colour of an edge between two cells which do not share any groups (this edge will be
    /// given the width of a cell border)
    disconnected_edge_color: RGB8,

    /// If the edge lengths are ~1 unit, how many units of space will be reserved round the edge of
    /// the SVG file
    padding: f32,
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
            font_size: 0.7,            // edge lengths
            text_vertical_nudge: 0.33, // multiples of `font_size`

            disconnected_edge_color: RGB8::new(255, 0, 0),

            padding: 0.5, // edge lengths
        }
    }
}
