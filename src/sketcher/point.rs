use std::{iter::Sum, ops::*};

#[inline]
fn round_to_precision(f: f32, prec: u8) -> f32 {
    f.mul_add(10.0f32.powi(prec as _), 0.5).floor() * 0.1f32.powi(prec as _)
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct PointF(pub f32, pub f32);
impl PointF {
    /// Return the square of the length of the vector.
    pub fn sq_length(self) -> f32 {
        self.0 * self.0 + self.1 * self.1
    }
    /// Return the length of the vector.
    pub fn length(self) -> f32 {
        let slen = self.sq_length();
        if slen > f32::EPSILON {
            slen.sqrt()
        } else {
            0.0
        }
    }
    /// Normalize the vector to have a length of 1 if it has a nonzero length.
    pub fn normalize(&mut self) {
        let len = self.length();
        if len > f32::EPSILON {
            self.0 /= len;
            self.1 /= len;
        }
    }
    /// Rotate a vector given the sine and cosine of an angle.
    pub fn rotate(&mut self, sin: f32, cos: f32) {
        let Self(x, y) = *self;
        self.0 = x * cos + y * sin;
        self.1 = -x * sin + y * cos;
    }
    /// Find the parallel component of this vector along a given axis.
    pub fn parallel_component(self, axis: Self) -> Self {
        let dot = self.dot(axis);
        axis * dot / axis.sq_length()
    }
    /// Round to a given precision.
    pub fn round(&mut self, prec: u8) {
        self.0 = round_to_precision(self.0, prec);
        self.1 = round_to_precision(self.1, prec);
    }

    /// Find the dot product of this vector with another.
    pub fn dot(self, other: Self) -> f32 {
        self.0 * other.0 + self.1 * other.1
    }
    /// Find the cross product of this vector with another.
    pub fn cross(self, other: Self) -> f32 {
        self.0 * other.1 - self.1 * other.0
    }
}
impl Add for PointF {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0, self.1 + rhs.1)
    }
}
impl Sub for PointF {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0, self.1 - rhs.1)
    }
}
impl Neg for PointF {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0, -self.1)
    }
}
impl Mul<f32> for PointF {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self::Output {
        Self(self.0 * rhs, self.1 * rhs)
    }
}
impl Div<f32> for PointF {
    type Output = Self;

    fn div(self, rhs: f32) -> Self::Output {
        Self(self.0 / rhs, self.1 / rhs)
    }
}
impl AddAssign for PointF {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
        self.1 += rhs.1;
    }
}
impl SubAssign for PointF {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0;
        self.1 -= rhs.1;
    }
}
impl MulAssign<f32> for PointF {
    fn mul_assign(&mut self, rhs: f32) {
        self.0 *= rhs;
        self.1 *= rhs;
    }
}
impl DivAssign<f32> for PointF {
    fn div_assign(&mut self, rhs: f32) {
        self.0 /= rhs;
        self.1 /= rhs;
    }
}
impl Sum for PointF {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self(0.0, 0.0), Add::add)
    }
}
