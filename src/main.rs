use kuudos::{gen_svg_string, RenderingOpts, Shape};

fn main() {
    let s = Shape::square(3, 2);

    // Write the shape to an SVG file
    let svg = gen_svg_string(&s, &RenderingOpts::default(), 20.0);
    println!("{}", svg);
    std::fs::write("classic.svg", svg).unwrap();
}
