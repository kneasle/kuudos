use itertools::Itertools;
use kuudos::{
    solve::{clues_from_str, solve},
    svg::{gen_svg_string, RenderingOpts},
    Shape,
};

const VALUE_NAMES: &str = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

fn main() {
    let (shape, clues) = if true {
        let s = Shape::star2x2(5, std::f32::consts::PI);
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
#[allow(dead_code)]
enum DisplayType {
    Clues,
    Solution,
    CellNames,
}
