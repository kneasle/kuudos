//! Miscellaneous utility functions, usually geometric operations.

use std::{cmp::Ordering, f32::consts::PI};

use crate::V2;

use angle::{Angle, Rad};
use num_complex::Complex32;

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

/// The angle (going **clockwise**) between the directions of two vectors.  This is always in the
/// region `-PI < x <= PI`.
pub fn angle_between(v1: V2, v2: V2) -> Rad<f32> {
    // The clockwise angle from the x-axis to the direction of `v1` and `v2`
    let angle1 = f32::atan2(v1.y, v1.x);
    let angle2 = f32::atan2(v2.y, v2.x);

    let mut angle_diff = angle2 - angle1;
    // Rotate angle_diff into the range `-2 * PI <= x < 2 * PI`
    angle_diff %= PI * 2.0;
    // If the absolute difference is bigger than PI (i.e. >180 degrees) then rotate it back into
    // the range `-PI < x <= PI`
    if angle_diff > PI {
        angle_diff -= PI * 2.0;
    }
    if angle_diff <= -PI {
        angle_diff += PI * 2.0;
    }

    Rad(angle_diff)
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

/// Returns the radius of the largest circle which fits within a regular polygon with a given side
/// length.
pub fn regular_polygon_inradius(sides: usize, side_length: f32) -> f32 {
    /* We do trig on the following triangle
     *
     *          half-edge
     *          +--------+
     *          |_|     /
     *          |      /
     *          |     /
     * inradius |    /
     *          |   / outradius
     *          |  /
     *          | /
     *          |/  <---- angle = 360deg / sides * 2 = (PI / sides) radians
     *          + centre
     */
    let angle = PI / sides as f32;
    side_length * 0.5 / angle.tan()
}

////////////////////////
// UTILITY DATA TYPES //
////////////////////////

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

/// An axis-aligned 2D rectangle
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rect2 {
    min: V2,
    max: V2,
}

#[allow(dead_code)]
impl Rect2 {
    /* CONSTRUCTORS */

    pub fn from_min_max(min: V2, max: V2) -> Self {
        Self { min, max }
    }

    pub fn from_min_size(min: V2, size: V2) -> Self {
        Self {
            min,
            max: min + size,
        }
    }

    pub fn from_centre_size(centre: V2, size: V2) -> Self {
        Self {
            min: centre - size * 0.5,
            max: centre + size * 0.5,
        }
    }

    /// Creates a zero-size [`Rect2`] at a given point
    pub fn point(pt: V2) -> Self {
        Self::from_min_size(pt, V2::ZERO)
    }

    /// Returns the bounding box of a set of points.  Returns `None` if the iterator yields no
    /// points.
    pub fn bbox(points: impl IntoIterator<Item = V2>) -> Option<Self> {
        let mut point_iter = points.into_iter();

        let mut bbox = Self::point(point_iter.next()?);
        for pt in point_iter {
            bbox.expand_to_include(pt);
        }
        Some(bbox)
    }

    /* GETTERS */

    pub fn min(self) -> V2 {
        self.min
    }

    pub fn max(self) -> V2 {
        self.max
    }

    pub fn size(self) -> V2 {
        self.max - self.min
    }

    pub fn centre(self) -> V2 {
        V2::lerp(self.min, self.max, 0.5)
    }

    /* OPERATIONS */

    /// Creates a `Rect2` that contains both `self` and `other`
    pub fn union(self, other: Rect2) -> Self {
        Self {
            min: min_v2(self.min, other.min),
            max: max_v2(self.max, other.max),
        }
    }

    /// Creates a `Rect2` that contains both `self` and `other`
    pub fn union_iter(rects: impl IntoIterator<Item = Rect2>) -> Option<Self> {
        let mut rect_iter = rects.into_iter();

        let mut combined_union = rect_iter.next()?;
        for rect in rect_iter {
            combined_union = combined_union.union(rect);
        }
        Some(combined_union)
    }

    /// Expands `self` to so that it contains some point
    pub fn expand_to_include(&mut self, pt: V2) {
        self.min = min_v2(self.min, pt);
        self.max = max_v2(self.max, pt);
    }
}

/// Component-wise min of two vectors
pub fn min_v2(a: V2, b: V2) -> V2 {
    V2::new(min_f32(a.x, b.x), min_f32(a.y, b.y))
}

/// Component-wise max of two vectors
pub fn max_v2(a: V2, b: V2) -> V2 {
    V2::new(max_f32(a.x, b.x), max_f32(a.y, b.y))
}

/// Returns the smaller of `a` and `b`, panicking if either are NaN
pub fn min_f32(a: f32, b: f32) -> f32 {
    match a.partial_cmp(&b) {
        Some(Ordering::Less) => a,
        Some(Ordering::Equal) => a,
        Some(Ordering::Greater) => b,
        None => panic!(),
    }
}

/// Returns the smaller of `a` and `b`, panicking if either are NaN
pub fn max_f32(a: f32, b: f32) -> f32 {
    match a.partial_cmp(&b) {
        Some(Ordering::Less) => b,
        Some(Ordering::Equal) => a,
        Some(Ordering::Greater) => a,
        None => panic!(),
    }
}

#[cfg(test)]
mod tests {
    use std::fmt::Debug;

    use angle::{Angle, Deg};

    use crate::{V2Ext, V2};

    #[test]
    fn angle_between() {
        fn check(v1: V2, v2: V2, exp_angle: impl Angle<f32> + Debug) {
            assert_eq!(super::angle_between(v1, v2), exp_angle.to_rad());
        }

        // The angle from any direction to (any +ve scalar multiple of itself) is 0
        check(V2::LEFT, V2::LEFT, Deg(0.0));
        check(V2::DOWN, V2::DOWN, Deg(0.0));
        check(V2::ONE, V2::ONE * 2.5, Deg(0.0));
        // The angle from any direction to any -ve scalar multiple of itself is +180*
        check(V2::DOWN, V2::UP, Deg(180.0));
        check(V2::RIGHT, V2::LEFT, Deg(180.0));
        check(V2::ONE, -V2::ONE, Deg(180.0));
        check(V2::new(1.0, -1.0), V2::new(-1.0, 1.0), Deg(180.0));

        check(V2::UP, V2::RIGHT * 2.5, Deg(90.0));
        check(V2::DOWN, V2::RIGHT * 2.5, Deg(-90.0));
        check(V2::DOWN, -V2::ONE, Deg(135.0));
        check(V2::DOWN, V2::new(1.0, -1.0), Deg(-135.0));
    }
}
