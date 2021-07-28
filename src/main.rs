#![allow(dead_code)]

use angle::Deg;
use itertools::Itertools;
use kuudos::{
    builder::{Builder, EdgeLinkStyle},
    solve::{clues_from_str, solve},
    svg::{gen_svg_string, RenderingOpts},
    Shape, Side, V2Ext, V2,
};

const VALUE_NAMES: &str = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

fn main() {
    let (shape, clues) = if true {
        // let s = Shape::star2x2(5, PI);
        let builder = nine_mens_morris();
        std::fs::write("bdr.svg", builder.as_svg(40.0)).unwrap();

        let (s, _symm) = builder.build().unwrap();
        let clues = vec![None; s.num_cells()];
        (s, clues)
    } else {
        let s = Shape::classic();
        let clues = clues_from_str(
            // "004300209 005009001 070060043 006002087 190007400 050083000 600000105 003508690 042910300",
            // "..24..... .......98 ..8..16.. .5....826 ....4.7.. 7..8...5. ..167.... 3...92... ..9......",
            "4.....8.5 .3....... ...7..... .2.....6. ....8.4.. ....1.... ...6.3.7. 5..2..... 1.4......",
        );
        (s, clues)
    };

    let soln = solve(&shape, &clues, false).unwrap();

    let display_type = DisplayType::Solution;

    // Write the shape to an SVG file
    let svg = gen_svg_string(
        &shape,
        &RenderingOpts::default(),
        40.0,
        // Label the cells with alpha-numeric characters
        &match display_type {
            DisplayType::Clues => clues
                .iter()
                .map(|x| x.map(|v| VALUE_NAMES.chars().nth(v).unwrap()))
                .collect_vec(),
            DisplayType::Solution => soln
                .iter()
                .map(|digit| VALUE_NAMES.chars().nth(*digit).unwrap())
                .map(Some)
                .collect_vec(),
            DisplayType::CellNames => VALUE_NAMES
                .chars()
                .cycle()
                .take(shape.num_cells())
                .map(Some)
                .collect_vec(),
        },
    );
    std::fs::write("classic.svg", svg).unwrap();
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum DisplayType {
    Clues,
    Solution,
    CellNames,
}

fn nine_mens_morris() -> Builder {
    let mut bdr = Builder::new(3, 3, 4);

    // Add line(s) of connected cells
    let line_box_1 = bdr.add_box_square(V2::UP * 3.0, 1.0, Deg(0.0));
    let line_box_2 = bdr.extrude_edge(line_box_1, Side::Top).unwrap();
    let line_box_3 = bdr.extrude_edge(line_box_2, Side::Top).unwrap();
    // Add diagonal line(s) of 'floating' cells
    let up_right = V2::UP + V2::RIGHT;
    bdr.add_box_square(up_right * 3.0, 1.0, Deg(0.0));
    let corner_box_2 = bdr.add_box_square(up_right * 6.0, 1.0, Deg(0.0));
    let corner_box_3 = bdr.add_box_square(up_right * 9.0, 1.0, Deg(0.0));
    // Link the lines of boxes
    bdr.link_edges(
        line_box_2,
        Side::Right,
        corner_box_2,
        Side::Left,
        EdgeLinkStyle::Linear,
    );
    bdr.link_edges(
        line_box_3,
        Side::Right,
        corner_box_3,
        Side::Left,
        EdgeLinkStyle::Linear,
    );
    bdr.link_edges(
        corner_box_2,
        Side::Bottom,
        bdr.rotational_copy_of(line_box_2, 1).unwrap(),
        Side::Left,
        EdgeLinkStyle::Linear,
    );
    bdr.link_edges(
        corner_box_3,
        Side::Bottom,
        bdr.rotational_copy_of(line_box_3, 1).unwrap(),
        Side::Left,
        EdgeLinkStyle::Linear,
    );

    bdr
}

fn race_track() -> Builder {
    // Create a Builder with 3-fold symmetry
    let mut bdr = Builder::new(3, 3, 3);

    // Create a new parallelogram box to make a hexagon
    let hex_box = bdr.add_box_parallelogram(V2::ZERO, V2::UP, bdr.rotate_point_by_steps(V2::UP, 1));
    // Extrude a square face off one side
    let outer_box = bdr.extrude_edge(hex_box, Side::Top).unwrap();
    // Link the outer boxes together
    bdr.link_edges(
        outer_box,
        Side::Right,
        bdr.rotational_copy_of(outer_box, 1).unwrap(),
        Side::Left,
        EdgeLinkStyle::Arc,
    );

    // Rotate the hexagon so that the side we extruded faces upwards
    bdr.rotate(Deg(-30.0));
    bdr
}

fn triangle() -> Builder {
    let mut builder = Builder::new(3, 3, 3);
    let b = builder.add_box_square(V2::UP * 3.7, 1.0, Deg(-45.0));
    builder
        .connect_edges_with_box(
            b,
            Side::Bottom,
            builder.rotational_copy_of(b, 1).unwrap(),
            Side::Left,
        )
        .unwrap();
    builder
}
