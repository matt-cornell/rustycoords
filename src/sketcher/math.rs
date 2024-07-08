//! This is just the stuff in `sketcherMinimizerMaths`

use super::point::PointF;

pub fn intersection_of_segments(
    s1p1: PointF,
    s1p2: PointF,
    s2p1: PointF,
    s2p2: PointF,
) -> Option<PointF> {
    let p = s1p1;
    let r = s1p2 - s1p1;
    let q = s2p1;
    let s = s2p2 - s1p1;
    let rxs = r.cross(s);
    if rxs.abs() < f32::EPSILON {
        return None;
    }
    let qmp = q - p;
    let t = qmp.cross(r) / rxs;
    if t < 0.0 || t > 1.0 {
        return None;
    }
    let u = qmp.cross(r) / rxs;
    if u < 0.0 || u > 1.0 {
        return None;
    }
    Some(p + r * t)
}
pub fn signed_angle(p1: PointF, p2: PointF, p3: PointF) -> f32 {
    let v1 = p1 - p2;
    let v2 = p3 - p2;
    v1.cross(v2).atan2(v1.dot(v2)).to_degrees()
}
pub fn unsigned_angle(p1: PointF, p2: PointF, p3: PointF) -> f32 {
    // readable code goes brrrrr
    let v1 = p1 - p2;
    let v2 = p3 - p2;
    let l = (v1.sq_length() * v2.sq_length()).sqrt().max(f32::EPSILON);
    (v1.dot(v2) / l).acos().to_degrees()
}
pub fn points_coincide(p1: PointF, p2: PointF) -> bool {
    (p1 - p2).sq_length() < f32::EPSILON * f32::EPSILON
}
pub fn same_side(p1: PointF, p2: PointF, line_p1: PointF, line_p2: PointF) -> bool {
    let PointF(x, y) = line_p2 - line_p1;
    if x.abs() > y.abs() {
        let m = y / x;
        let d1 = p1.1 - line_p1.1 - m * (p1.0 - line_p1.0);
        let d2 = p2.1 - line_p1.1 - m * (p2.0 - line_p1.0);
        d1 * d2 > 0.0
    } else {
        let m = x / y;
        let d1 = p1.0 - line_p1.0 - m * (p1.1 * line_p1.1);
        let d2 = p2.0 - line_p1.0 - m * (p2.1 * line_p1.1);
        d1 * d2 > 0.0
    }
}
pub fn project_on_line(p: PointF, sp1: PointF, sp2: PointF) -> PointF {
    let l1 = p - sp1;
    let l3 = sp2 - sp1;
    let sl2 = l3.sq_length().max(f32::EPSILON);
    let t = l1.dot(l3) / sl2;
    sp1 + l3 * t
}
/// First value is normal return, second is returnT
pub fn sq_dist_point_segment(p: PointF, sp1: PointF, sp2: PointF) -> (f32, f32) {
    let l1 = p - sp1;
    let l2 = sp2 - p;
    let l3 = sp2 - sp1;
    let sl2 = l3.sq_length().max(f32::EPSILON);
    let t = l1.dot(l3) / sl2;
    let sqdist = if t < 0.0 {
        l1.sq_length()
    } else if t > 1.0 {
        l2.sq_length()
    } else {
        let proj = sp1 + l3 * t;
        let l5 = p - proj;
        l5.sq_length()
    };
    (sqdist.max(f32::EPSILON), t.clamp(0.0, 1.0))
}

/// Specialization -- a = c = [1, 1, 1...], b = [8, 4, ..., 4.25]
/// scratch.len() >= rhs.len()
/// out.len() >= rhs.len()
///
/// if rhs is None, it's assumed to be [-4, 0, ..., 1], and out is used for the length
fn tridiagonal_solve(rhs: Option<&[f32]>, out: &mut [f32], scratch: &mut [f32]) {
    let len = out.len();
    let b = |n: usize| {
        if n == 0 {
            8.0f32
        } else if n == rhs.map_or(len, |r| r.len()) - 1 {
            4.25
        } else {
            4.0
        }
    };
    let rhs = |n: usize| {
        if let Some(rhs) = rhs {
            rhs[n]
        } else if n == 0 {
            -4.0
        } else if n == len - 1 {
            1.0
        } else {
            0.0
        }
    };
    let n = out.len();
    let mut bet = b(0);
    out[0] = rhs(0) / bet;
    for j in 1..n {
        scratch[j] = 1.0 / bet;
        bet = b(0) - scratch[j];
        debug_assert_ne!(bet, 0.0);
        out[j] = (rhs(j) - out[j - 1]) / bet;
    }
    for j in 1..n {
        out[n - j - 1] -= scratch[n - j] * out[n - j];
    }
}
/// Specialization -- a = 1, b = 4, c = 1, that's the only way it's used
/// alpha = 1, beta = 1
/// out has at least the length of inputs
/// scratch has at least 2 times the length
fn cyclic_solve(rhs: &[f32], out: &mut [f32], scratch: &mut [f32]) {
    let n = rhs.len();
    debug_assert!(n > 2);
    debug_assert!(scratch.len() >= n * 2);
    let (tridiag_scratch, z) = scratch.split_at_mut(n);
    tridiagonal_solve(Some(rhs), out, tridiag_scratch);
    tridiagonal_solve(None, z, tridiag_scratch);
    let fact = (out[0] - out[n - 1] * 0.25) / (z[0] - z[n - 1] * 0.25);
    for (x, z) in out.iter_mut().zip(z) {
        *x -= *z * fact;
    }
}
pub fn point_on_cubic_bezier(p1: PointF, cp1: PointF, cp2: PointF, p2: PointF, t: f32) -> PointF {
    let v1 = p1 * (1.0 - t) + cp1 * t;
    let v2 = cp1 * (1.0 - t) + cp2 * t;
    let v3 = cp1 * (1.0 - t) + p2 * t;
    let v4 = v1 * (1.0 - t) + v2 * t;
    let v5 = v2 * (1.0 - t) + v3 * t;
    v4 * (1.0 - t) + v5 * t
}
/// first and second need same size as knots, scratch needs 5 times
pub fn closed_bezier_control_points(
    knots: &[PointF],
    first: &mut [PointF],
    second: &mut [PointF],
    scratch: &mut [f32],
) {
    let n = knots.len();
    if n < 2 {
        return;
    }
    debug_assert!(scratch.len() >= n * 8);
    let (rhs, scratch) = scratch.split_at_mut(n);
    let (xy, scratch) = scratch.split_at_mut(n);
    for i in 0..n {
        let j = (i + 1) % n;
        rhs[i] = knots[i].0 * 4.0 + knots[j].0 * 2.0;
    }
    cyclic_solve(rhs, xy, scratch);
    for (((first, second), x), knot) in first.iter_mut().zip(&mut *second).zip(&*xy).zip(knots) {
        first.0 = *x;
        second.0 = knot.0 * 2.0 - *x;
    }
    for i in 0..n {
        let j = (i + 1) % n;
        rhs[i] = knots[i].1 * 4.0 + knots[j].1 * 2.0;
    }
    cyclic_solve(rhs, xy, scratch);
    for (((first, second), y), knot) in first.iter_mut().zip(second).zip(xy).zip(knots) {
        first.1 = *y;
        second.1 = knot.1 * 2.0 - *y;
    }
}
pub fn mirror_point(pt: PointF, sp1: PointF, sp2: PointF) -> PointF {
    let sv = sp2 - sp1;
    let v2 = pt - sp1;
    let pc = v2.parallel_component(sv);
    let nc = v2 - pc;
    sp1 + pc - nc
}
/// WTF.
/// cutoff defaults to 4 btw
pub fn canonball_distance(
    origin: [f32; 3],
    direction: [f32; 3],
    target: [f32; 3],
    br: f32,
    tr: f32,
    cutoff: f32,
) -> f32 {
    let mut target_d = vecsub3(target, origin);
    let rr = br + tr;
    let d2 = len23(target_d);
    if d2 > (cutoff + rr) * (cutoff + rr) {
        return cutoff;
    }
    if d2 < rr * rr {
        return 0.0;
    }
    let d = d2.sqrt();
    if d > f32::EPSILON {
        target_d[0] /= d;
        target_d[1] /= d;
        target_d[2] /= d;
    }
    let cos = dot3(target_d, direction);
    if cos < 0.0 {
        return cutoff;
    }
    let sin = (1.0 - cos * cos).sqrt();
    let f = d * sin;
    if f > rr {
        return cutoff;
    }
    let result = (d2 - f * f).sqrt() - ((rr * rr) - (f * f)).sqrt();
    result.min(cutoff)
}
/// length squared, in three dimensions
#[inline]
fn len23(vec: [f32; 3]) -> f32 {
    dot3(vec, vec)
}
#[inline]
pub fn dot3(lhs: [f32; 3], rhs: [f32; 3]) -> f32 {
    lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2]
}
pub fn cross3(lhs: [f32; 3], rhs: [f32; 3]) -> [f32; 3] {
    [
        lhs[1] * rhs[2] - lhs[2] * rhs[1],
        lhs[2] * rhs[0] - lhs[0] * rhs[2],
        lhs[0] * rhs[1] - lhs[1] * rhs[0],
    ]
}
#[inline]
pub fn length3(vec: [f32; 3]) -> f32 {
    len23(vec).sqrt()
}
#[inline]
pub fn vecsub3(lhs: [f32; 3], rhs: [f32; 3]) -> [f32; 3] {
    [lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]]
}
#[inline]
pub fn dist3(a: [f32; 3], b: [f32; 3]) -> f32 {
    length3(vecsub3(a, b))
}
pub fn angle3(p1: [f32; 3], p2: [f32; 3], p3: [f32; 3]) -> f32 {
    let a = vecsub3(p1, p2);
    let b = vecsub3(p3, p2);
    let l1 = len23(a); // multiply the lengths first before taking the sqrt
    let l2 = len23(b);
    let dp = dot3(a, b);
    (dp / (l1 * l2).sqrt()).acos().to_degrees()
}
/// yeah idk here good luck
pub fn dihedral3(p1: [f32; 3], p2: [f32; 3], p3: [f32; 3], p4: [f32; 3]) -> f32 {
    let a = cross3(vecsub3(p1, p2), vecsub3(p3, p2));
    let b = cross3(vecsub3(p2, p3), vecsub3(p4, p3));
    angle3(a, [0.0, 0.0, 0.0], b)
}
