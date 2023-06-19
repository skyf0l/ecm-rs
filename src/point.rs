use rug::Integer;

/// Montgomery form of Points in an elliptic curve.
///
/// In this form, the addition and doubling of points
/// does not need any y-coordinate information thus
/// decreasing the number of operations.
/// Using Montgomery form we try to perform point addition
/// and doubling in least amount of multiplications.
///
/// The elliptic curve used here is of the form
/// `(E : b*y**2*z = x**3 + a*x**2*z + x*z**2)`.
/// The `a_24` parameter is equal to `(a + 2)/4`.
///
/// References
/// ----------
/// - http://www.hyperelliptic.org/tanja/SHARCS/talks06/Gaj.pdf
#[derive(Debug, Clone, Default)]
pub struct Point {
    /// X coordinate of the Point
    pub x_cord: Integer,
    /// Z coordinate of the Point
    pub z_cord: Integer,
    /// Parameter of the elliptic curve in Montgomery form
    pub a_24: Integer,
    /// modulus
    pub modulus: Integer,
}

impl Point {
    /// Initial parameters for the Point struct.
    ///
    /// # Parameters
    ///
    /// - `x_cord`: X coordinate of the Point
    /// - `z_cord`: Z coordinate of the Point
    /// - `a_24`: Parameter of the elliptic curve in Montgomery form
    /// - `mod`: modulus
    pub fn new(x_cord: Integer, z_cord: Integer, a_24: Integer, modulus: Integer) -> Point {
        Point {
            x_cord,
            z_cord,
            a_24,
            modulus,
        }
    }

    /// Adds two points `self` and `Q` where `diff = self - Q`.
    ///
    /// This algorithm requires 6 multiplications. The assumption is that `self.x_cord * Q.x_cord * (self.x_cord - Q.x_cord) != 0`.
    /// Using this algorithm speeds up the addition by reducing the number of multiplications required.
    ///
    /// The `mont_ladder` algorithm is constructed in a way that the difference between intermediate points is always equal to the initial point.
    /// So, we always know what the difference between the point is.
    ///
    /// # Parameters
    ///
    /// - `Q`: Point on the curve in Montgomery form.
    /// - `diff`: `self - Q`
    pub fn add(&self, q: &Point, diff: &Point) -> Point {
        let u = (self.x_cord.clone() - self.z_cord.clone()) * (q.x_cord.clone() + q.z_cord.clone());
        let v = (self.x_cord.clone() + self.z_cord.clone()) * (q.x_cord.clone() - q.z_cord.clone());
        let add = u.clone() + v.clone();
        let subt = u - v;
        let x_cord = (diff.z_cord.clone() * add.clone() * add) % self.modulus.clone();
        let z_cord = (diff.x_cord.clone() * subt.clone() * subt) % self.modulus.clone();

        Point::new(x_cord, z_cord, self.a_24.clone(), self.modulus.clone())
    }

    /// Doubles a point in an elliptic curve in Montgomery form.
    pub fn double(&self) -> Point {
        let u = (self.x_cord.clone() + self.z_cord.clone()).square();
        let v = (self.x_cord.clone() - self.z_cord.clone()).square();
        let diff = u.clone() - v.clone();
        let x_cord = (u * v.clone()) % self.modulus.clone();
        let z_cord = (diff.clone() * (v + self.a_24.clone() * diff)) % &self.modulus;

        Point::new(x_cord, z_cord, self.a_24.clone(), self.modulus.clone())
    }

    /// Scalar multiplication of a point in Montgomery form
    /// using Montgomery Ladder Algorithm.
    /// A total of 11 multiplications are required in each step of this
    /// algorithm.
    ///
    /// # Parameters
    ///
    /// - `k`: The positive integer multiplier
    pub fn mont_ladder(&self, k: Integer) -> Point {
        let mut q = self.clone();
        let mut r = self.double();

        for i in format!("{:b}", k)[1..].chars() {
            if i == '1' {
                q = r.add(&q, self);
                r = r.double();
            } else {
                r = q.add(&r, self);
                q = q.double();
            }
        }
        q
    }
}

impl PartialEq for Point {
    /// Two points are equal if X/Z of both points are equal.
    fn eq(&self, other: &Self) -> bool {
        if self.a_24 != other.a_24 || self.modulus != other.modulus {
            false
        } else {
            self.x_cord.clone() * self.z_cord.clone().invert(&self.modulus).unwrap() % &self.modulus
                == other.x_cord.clone() * other.z_cord.clone().invert(&self.modulus).unwrap()
                    % &self.modulus
        }
    }
}
