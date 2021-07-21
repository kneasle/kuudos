use kuudos::{gen_svg_string, RenderingOpts, Shape};

fn main() {
    let s = Shape::square(3, 3);

    // Write the shape to an SVG file
    let svg = gen_svg_string(&s, &RenderingOpts::default(), 40.0);
    std::fs::write("classic.svg", svg).unwrap();
}
