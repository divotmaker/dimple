use std::ops::{Add, Div, Mul, Neg, Sub, SubAssign};

/// A 3D vector of `f64`, in the crate frame: +X downrange, +Y up, +Z right
/// of the target line.
///
/// Convert to and from your own math library through `[f64; 3]`:
/// `Vec3::from([x, y, z])` and `<[f64; 3]>::from(v)`.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Vec3 {
    /// Downrange component.
    pub x: f64,
    /// Vertical component, positive up.
    pub y: f64,
    /// Lateral component, positive right of the target line.
    pub z: f64,
}

impl Vec3 {
    /// All zeros.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);
    /// Unit vector downrange.
    pub const X: Self = Self::new(1.0, 0.0, 0.0);
    /// Unit vector up.
    pub const Y: Self = Self::new(0.0, 1.0, 0.0);
    /// Unit vector right of the target line.
    pub const Z: Self = Self::new(0.0, 0.0, 1.0);

    /// A vector from its components.
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Dot product.
    #[must_use]
    pub fn dot(self, rhs: Self) -> f64 {
        (self.x * rhs.x) + (self.y * rhs.y) + (self.z * rhs.z)
    }

    /// Cross product, right-handed.
    #[must_use]
    pub fn cross(self, rhs: Self) -> Self {
        Self {
            x: self.y * rhs.z - rhs.y * self.z,
            y: self.z * rhs.x - rhs.z * self.x,
            z: self.x * rhs.y - rhs.x * self.y,
        }
    }

    /// Length.
    #[must_use]
    pub fn length(self) -> f64 {
        self.dot(self).sqrt()
    }

    /// Squared length.
    #[must_use]
    pub fn length_squared(self) -> f64 {
        self.dot(self)
    }

    /// The unit vector in the same direction. Non-finite for a zero vector.
    #[must_use]
    pub fn normalize(self) -> Self {
        self * self.length().recip()
    }

    /// The unit vector in the same direction, or [`ZERO`](Self::ZERO) when
    /// the length is zero or not finite.
    #[must_use]
    pub fn normalize_or_zero(self) -> Self {
        let rcp = self.length().recip();
        if rcp.is_finite() && rcp > 0.0 {
            self * rcp
        } else {
            Self::ZERO
        }
    }

    /// Linear interpolation: `self` at `s = 0`, `rhs` at `s = 1`.
    #[must_use]
    pub fn lerp(self, rhs: Self, s: f64) -> Self {
        self * (1.0 - s) + rhs * s
    }
}

impl From<[f64; 3]> for Vec3 {
    fn from([x, y, z]: [f64; 3]) -> Self {
        Self { x, y, z }
    }
}

impl From<Vec3> for [f64; 3] {
    fn from(v: Vec3) -> Self {
        [v.x, v.y, v.z]
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl SubAssign for Vec3 {
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;
    fn mul(self, rhs: Vec3) -> Vec3 {
        Vec3::new(self * rhs.x, self * rhs.y, self * rhs.z)
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}
