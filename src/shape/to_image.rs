use angle::Rad;
use itertools::Itertools;

use super::{CellIdx, CellVec, Shape};
use crate::{
    image::{Elem, FillStyle, Image, StrokeStyle, Style, TextStyle},
    utils, V2,
};

/// Write a populated sudoku grid to an SVG string
pub fn gen_image_with_clues(shape: &Shape, clues: &CellVec<Option<char>>) -> Image {
    gen_image(shape, single_digit_per_cell(TextStyle::Clue, clues))
}

pub fn gen_image(
    shape: &Shape,
    mut get_cell_contents: impl FnMut(CellIdx, &[V2]) -> Vec<Elem>,
) -> Image {
    // Create image
    let mut image = Image::empty();

    // Add extra elements behind the cells
    for extra_elem in &shape.extra_elements {
        image.add(extra_elem.clone());
    }

    // Add the cell backgrounds
    for vert_idxs in shape.cells.iter() {
        let verts = vert_idxs.iter().map(|idx| shape.verts[*idx]).collect_vec();
        image.add(Elem::Polygon(verts, Style::JustFill(FillStyle::Background)));
    }

    // Add the cell contents
    for (cell_idx, vert_idxs) in shape.cells.iter_enumerated() {
        let verts = vert_idxs.iter().map(|idx| shape.verts[*idx]).collect_vec();
        image.add_iter(get_cell_contents(cell_idx, &verts));
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
                        StrokeStyle::Disconnected // Adjacent cells which share no groups
                    } else if shape.do_cells_share_box(cell_left, edge.cell_idx_right) {
                        StrokeStyle::CellBorder // Adjacent cells which are in the same box
                    } else {
                        StrokeStyle::BoxBorder // Adjacent cells which share a group but aren't in
                                               // the same box
                    }
                }
                None => StrokeStyle::BoxBorder, // External edges
            };
            (edge_type, edge)
        })
        .collect_vec();
    // Sort the edges by type with `Disconnected` < `CellBorder` < `BoxBorder` (i.e. add bolder
    // borders higher up the image stack).
    edges.sort_by_key(|(ty, _edge)| match ty {
        StrokeStyle::Disconnected => 0,
        StrokeStyle::CellBorder => 1,
        StrokeStyle::BoxBorder => 2,
        &StrokeStyle::Custom(_) => unreachable!("`Shape::gen_image` won't generate custom styles"),
    });

    // Add the edges
    for (style, edge) in edges {
        // Extract the coordinates of the top/bottom vertices of this edge
        let vert_top = shape.verts[edge.vert_idx_top];
        let vert_bottom = shape.verts[edge.vert_idx_bottom];
        image.add(Elem::LineSegment(vert_bottom, vert_top, style));
    }
    image
}

///////////////////////////////////////////
// HELPER FUNCTIONS FOR POPULATING CELLS //
///////////////////////////////////////////

/// Returns a closure that generates a places a single 'penned' digit in the middle of each cell.
pub fn single_digit_per_cell(
    text_style: TextStyle,
    digits: &CellVec<Option<char>>,
) -> impl FnMut(CellIdx, &[V2]) -> Vec<Elem> + '_ {
    move |idx, verts| match digits[idx] {
        Some(c) => char_in_cell_center(text_style.clone(), verts, c),
        None => vec![],
    }
}

/// Generates the image [`Elem`]ents to place a single [`char`] in the center of a cell bounded by
/// some vertices.
pub fn char_in_cell_center(text_style: TextStyle, verts: &[V2], digit: char) -> Vec<Elem> {
    let cell_centre =
        utils::centroid(verts.iter().copied()).expect("Can't have cell with no versions");
    vec![Elem::Text {
        position: cell_centre,
        text: digit.to_string(),
        angle: Rad(0.0),
        style: text_style,
    }]
}
