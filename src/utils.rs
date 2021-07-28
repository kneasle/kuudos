//! Miscellaneous utility functions, usually geometric operations.

use crate::V2;

use angle::{Angle, Rad};
use num_complex::Complex32;

/// Returns the bounding box of a set of points as a (min, max) pair of vectors.  Returns `None` if
/// the iterator didn't yield any points.
pub fn bbox(points: impl IntoIterator<Item = V2>) -> Option<(V2, V2)> {
    let mut is_iter_empty = true;
    // If the shape has no vertices, then the bounding box isn't defined
    let mut min_x = f32::MAX;
    let mut min_y = f32::MAX;
    let mut max_x = f32::MIN;
    let mut max_y = f32::MIN;
    for v in points {
        is_iter_empty = false;
        if v.x < min_x {
            min_x = v.x;
        }
        if v.y < min_y {
            min_y = v.y;
        }
        if v.x > max_x {
            max_x = v.x;
        }
        if v.y > max_y {
            max_y = v.y;
        }
    }
    if is_iter_empty {
        // If the iterator yielded no elements, then the bbox is not defined
        None
    } else {
        Some((V2::new(min_x, min_y), V2::new(max_x, max_y)))
    }
}

/// Compute the average of a set of vectors
pub fn centroid(vs: impl IntoIterator<Item = V2>) -> Option<V2> {
    let mut sum = V2::ZERO;
    let mut num_elems = 0usize;
    for v in vs {
        sum += v;
        num_elems += 1;
    }
    if num_elems == 0 {
        // The centroid of no points is undefined
        None
    } else {
        Some(sum / num_elems as f32)
    }
}

/// Computes the centre and radius of the circle which passes through `p1` with tangent `tang_1`,
/// and also passes through `p2`.  Returns `None` if `p2` lies on the tangent.
pub fn circle_passing_through(p1: V2, tang_1: V2, p2: V2) -> Option<(V2, f32, Rad<f32>, Rad<f32>)> {
    // Given two points on a circle, the centre must lie on the locus of points equidistant from
    // each point.  Given that we have a tangent and therefore its normal, the centre of the circle
    // will be the intersection between these lines.
    let eq_line_direction = (p2 - p1).normal(); // The normal of p1 -> p2
    let eq_line_pt = V2::lerp(p1, p2, 0.5); // The midpoint of (p1, p2)
    let normal_to_tangent = tang_1.normal();
    // Compute the centre and radius
    let centre = intersect_lines(p1, normal_to_tangent, eq_line_pt, eq_line_direction)?;
    let radius = (p1 - centre).length();
    // Compute the angles `p1` and `p2`
    let start_angle = Rad((p1 - centre).angle());
    let end_angle = Rad((p2 - centre).angle());
    Some((centre, radius, start_angle, end_angle))
}

/// Computes the centre and the maximum and minimum radii (in that order) of an ellipse which
/// passes through two points in given directions.  Returns `None` if the directions are parallel
/// and therefore no such (unique) ellipse exists.
#[allow(dead_code)]
pub fn ellipse_passing_through(p1: V2, tang_1: V2, p2: V2, tang_2: V2) -> Option<(V2, V2, V2)> {
    let tangent_intersection = intersect_lines(p1, tang_1, p2, tang_2)?;
    let symmetry_over_p1 = V2::lerp(tangent_intersection, p1, 2.0);
    let symmetry_over_p2 = V2::lerp(tangent_intersection, p2, 2.0);
    // We are trying to find the Steiner in-ellipse (see
    // https://en.wikipedia.org/wiki/Steiner_inellipse) of the triangle with vertices
    // `(tangent_intersection, symmetry_over_p1, symmetry_over_p2)`.  This fact is noted in
    // [this stackexchange answer](https://math.stackexchange.com/a/2668084).
    //
    // The centre of a Steiner in-ellipse is the centroid of the triangle:
    let centre = centroid(vec![
        tangent_intersection,
        symmetry_over_p1,
        symmetry_over_p2,
    ])
    .unwrap();
    // Computing the major and minor axis is much more difficult, and we do it in two parts:
    // 1. Compute the direction of the major axis by computing one of the ellipse's foci
    // 2. Compute the lengths of the axes directly (using the equations from the wiki page)
    //
    // The foci are computed using Marden's Theorem
    // (https://en.wikipedia.org/wiki/Marden%27s_theorem), which states that:
    //
    // > If a polynomial has complex roots at the vertices of our triangle,
    // > then the foci of the Steiner in-ellipse are the complex roots of the derivative of that
    //   polynomial.
    //
    // If the vertices of our triangle correspond to complex numbers `a`, `b` and `c` then we can
    // easily construct a polynomial `p(z)` with these roots:
    //
    //     p(z) = (z - a)(z - b)(z - c)
    //
    // Any non-zero scalar multiplication of this function also works, but this is simpler.  Now we
    // need the derivative for this function, which we could do by hand - or we could simply give
    // the problem to Wolfram Alpha with
    // https://www.wolframalpha.com/input/?i=derivative+of+%28x+-+a%29%28x+-+b%29%28x+-+c%29.  Not
    // only does Wolfram Alpha give the derivative, it also computes the formula for the roots:
    //
    //       p'(z) = 0
    //
    //     when 3z = a + b + c ± sqrt(a^2 + b^2 + c^2 - (ab + bc + ca))
    //
    // We only need one of the roots, so we arbitrarily choose whichever one is returned by
    // `num_complex`.
    let a = Complex32::new(tangent_intersection.x, tangent_intersection.y);
    let b = Complex32::new(symmetry_over_p1.x, symmetry_over_p1.y);
    let c = Complex32::new(symmetry_over_p2.x, symmetry_over_p2.y);
    let val_in_sqrt = a.powu(2) + b.powu(2) + c.powu(2) - (a * b + b * c + c * a);
    let three_z = a + b + c + val_in_sqrt.sqrt();
    let z = three_z / 3.0;
    let focus = V2::new(z.re, z.im); // Convert back out of complex numbers

    // If the focus is on top of the centre, then the ellipse is actually a circle
    if (focus - centre).length_squared() < 0.00001 * 0.00001 {
        let radius = (p1 - centre).length(); // Both p1 and p2 lie on the circle
        return Some((centre, V2::UP * radius, V2::DOWN));
    }

    let _direction_of_major_axis = (focus - centre).normalise();
    todo!()
}

/// Compute the intersection between two lines, returning `None` if the lines are parallel
pub fn intersect_lines(p1: V2, d1: V2, p2: V2, d2: V2) -> Option<V2> {
    /// Convert both lines into the form `ax + by = c`
    fn convert_line(p: V2, d: V2) -> (f32, f32, f32) {
        let normal = d.normal();
        // Every point `q` on the line must satisfy `(p - q).dot(normal) = 0`
        //  => p.dot(normal) - q.dot(normal) = 0
        //  =>             0 - q.dot(normal) = 0 - p.dot(normal)
        //  =>                 q.dot(normal) = p.dot(normal)
        // `q.dot(normal)` expands to `q.x * normal.x + q.y * normal.x`, so:
        //   q.x * normal.x + q.y * normal.x = p.dot(normal)
        //   q.x * a        + q.y * b        = c
        //
        //  => a = normal.x, b = normal.y, c = p.dot(normal)
        (normal.x, normal.y, V2::dot(p, normal))
    }

    let (a, b, c) = convert_line(p1, d1); // ax + by = c
    let (d, e, f) = convert_line(p2, d2); // dx + ey = f

    // At this point, we are solving the linear simultaneous equations
    //      `ax + by = c`
    //  and `dx + ey = f`.
    // These reduces to the 2x2 matrix equation
    //
    //   |a  b|   |x|   |c|
    //   |    | * | | = | |
    //   |d  e|   |y|   |f|
    //
    // and re-arranging gives us:
    //
    //            |x|   |a  b|-1   |c|
    //            | | = |    |   * | |
    //            |y|   |d  e|     |f|
    //
    // expanding the matrix inverse and re-arranging:
    //
    //            |x|      1    | e  -b|   |c|
    //            | | = ------- |      | * | |
    //            |y|   ae - bd |-d   a|   |f|
    //
    //                     1    |  ec - bf |
    //                = ------- |          |
    //                  ae - bd | -dc + af |
    //
    //                     1      | e|    |-b|
    //                = ------- (c|  | + f|  |)
    //                  ae - bd   |-d|    | a|
    let determinant = a * e - b * d;
    if determinant == 0.0 {
        None
    } else {
        let part_in_brackets = V2::new(e, -d) * c + V2::new(-b, a) * f;
        Some(part_in_brackets / determinant)
    }
}

/// Generate the SVG path string for a circular arc
pub fn svg_circle_arc_path_str(
    centre: V2,
    radius: f32,
    start_angle: impl Angle<f32> + Copy,
    end_angle: impl Angle<f32> + Copy,
) -> String {
    let right_radius = V2::RIGHT * radius;
    let start_pt = centre + right_radius.rotate(start_angle);
    let end_pt = centre + right_radius.rotate(end_angle);
    format!(
        "M {} {} A {} {} 0 0 1 {} {}",
        start_pt.x, start_pt.y, radius, radius, end_pt.x, end_pt.y
    )
}

/// Extension trait to add more methods to [`V2`]
pub trait V2Ext {
    const ZERO: Self;
    const ONE: Self;

    const UP: Self;
    const DOWN: Self;
    const RIGHT: Self;
    const LEFT: Self;

    fn rotate(self, angle: impl Angle<f32> + Copy) -> Self;
}

impl V2Ext for V2 {
    const ZERO: Self = V2 { x: 0.0, y: 0.0 };
    const ONE: Self = V2 { x: 1.0, y: 1.0 };

    const UP: Self = V2 { x: 0.0, y: -1.0 };
    const DOWN: Self = V2 { x: 0.0, y: 1.0 };
    const RIGHT: Self = V2 { x: 1.0, y: 0.0 };
    const LEFT: Self = V2 { x: -1.0, y: 0.0 };

    /// Rotates a vector **clockwise** by an angle in radians
    fn rotate(self, angle: impl Angle<f32> + Copy) -> Self {
        let sin = angle.sin();
        let cos = angle.cos();
        // Rotation **clockwise** corresponds to multiplication by the following matrix (which looks
        // like the classic anti-clockwise matrix because our y-axis goes down where the one in maths
        // goes up):
        // | cos(angle)  -sin(angle) |
        // | sin(angle)   cos(angle) |
        V2 {
            x: self.x * cos - self.y * sin,
            y: self.x * sin + self.y * cos,
        }
    }
}
