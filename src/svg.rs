use itertools::Itertools;
use rgb::RGB8;
use simple_xml_builder::XMLElement;

use crate::{Shape, V2};

/// Code to write sudoku grids to SVG files
pub fn gen_svg_string(shape: &Shape, opts: &RenderingOpts, scaling: f32) -> String {
    // This bounding box is in **untransformed** space
    let (bbox_min, bbox_max) = shape.bbox().expect("Shape should have at least one vertex");
    let padding_vec = V2::new(opts.padding, opts.padding);

    // Transform the vertices so that their min-point is (padding, padding) and they're scaled up
    // by the `scaling` factor
    let transformed_verts = shape
        .verts
        .iter()
        .map(|&v| (v - bbox_min + padding_vec) * scaling)
        .collect_vec();
    // Compute the dimensions of the resulting SVG image
    let img_dimensions = (bbox_max - bbox_min + padding_vec * 2.0) * scaling;

    // Root SVG element
    let mut root = XMLElement::new("svg");
    root.add_attribute("width", &img_dimensions.x.to_string());
    root.add_attribute("height", &img_dimensions.y.to_string());

    // Add the cell backgrounds
    for vert_idxs in &shape.cell_verts {
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

    // Sort the edges so that all the box borders come after all the cell borders (to prevent
    // overlap issues)
    let mut edges = shape
        .edges()
        .into_iter()
        .map(|edge| {
            let is_box_border = match edge.cell_idx_left {
                // If the edge sits between two faces, then it is internal and is therefore bold iff it
                // is on the border between two boxes (i.e. is only connected by a row or column).
                Some(cell_left) => {
                    shape.num_groups_shared_between(cell_left, edge.cell_idx_right) == 1
                }
                // Edges with only one cell must be at the edge of the grid, and are therefore bold
                None => true,
            };
            (is_box_border, edge)
        })
        .collect_vec();
    // Sort the edges by `is_box_border`, using the fact that `true` > `false`
    edges.sort_by_key(|(is_box_border, _edge)| *is_box_border);

    // Add the edges
    let box_border_width_str = (opts.box_border_width * scaling).to_string();
    let cell_border_width_str = (opts.cell_border_width * scaling).to_string();
    for (is_box_border, edge) in edges {
        // Extract the coordinates of the top/bottom vertices of this edge
        let vert_top = transformed_verts[edge.vert_idx_top];
        let vert_bottom = transformed_verts[edge.vert_idx_bottom];
        // Determine the edge's style, based on whether or not it's a box border
        let (color, width_str): (_, &str) = if is_box_border {
            (opts.box_border_color, &box_border_width_str)
        } else {
            (opts.cell_border_color, &cell_border_width_str)
        };

        let mut edge_elem = XMLElement::new("line");
        // Coords
        edge_elem.add_attribute("x1", &vert_top.x.to_string());
        edge_elem.add_attribute("y1", &vert_top.y.to_string());
        edge_elem.add_attribute("x2", &vert_bottom.x.to_string());
        edge_elem.add_attribute("y2", &vert_bottom.y.to_string());
        // Stroke style
        edge_elem.add_attribute("stroke", &color.to_string());
        edge_elem.add_attribute("stroke-linecap", "round");
        edge_elem.add_attribute("stroke-width", width_str);
        root.add_child(edge_elem);
    }

    root.to_string()
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

            padding: 0.5, // edge lengths
        }
    }
}
