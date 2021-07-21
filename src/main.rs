use itertools::Itertools;
use kuudos::{
    svg::{gen_svg_string, RenderingOpts},
    Shape,
};

fn main() {
    let s = Shape::star2x2(5, std::f32::consts::PI);
    // let s = Shape::classic();

    // Write the shape to an SVG file
    let svg = gen_svg_string(
        &s,
        &RenderingOpts::default(),
        40.0,
        // Label the cells with alpha-numeric characters
        &"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
            .chars()
            .cycle()
            .take(s.num_cells())
            .map(Some)
            .collect_vec(),
    );
    std::fs::write("classic.svg", svg).unwrap();
}
