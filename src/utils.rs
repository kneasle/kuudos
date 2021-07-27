//! Miscellaneous utility functions, usually related to vectors.

use crate::V2;

use angle::Angle;

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
    let mut sum = V2::new(0.0, 0.0);
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

/// Rotates a vector **clockwise** by an angle in radians
pub fn rotate_vec(v: V2, angle: impl Angle<f32> + Copy) -> V2 {
    let sin = angle.sin();
    let cos = angle.cos();
    // Rotation **clockwise** corresponds to multiplication by the following matrix (which looks
    // like the classic anti-clockwise matrix because our y-axis goes down where the one in maths
    // goes up):
    // | cos(angle)  -sin(angle) |
    // | sin(angle)   cos(angle) |
    V2 {
        x: v.x * cos - v.y * sin,
        y: v.x * sin + v.y * cos,
    }
}
