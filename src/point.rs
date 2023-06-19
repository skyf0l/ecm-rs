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
        let u = Integer::from(&self.x_cord - &self.z_cord) * Integer::from(&q.x_cord + &q.z_cord);
        let v = Integer::from(&self.x_cord + &self.z_cord) * Integer::from(&q.x_cord - &q.z_cord);
        let add = Integer::from(&u + &v);
        let subt = u - v;
        let x_cord = Integer::from(&diff.z_cord * &add) * &add % &self.modulus;
        let z_cord = Integer::from(&diff.x_cord * &subt) * &subt % &self.modulus;

        Point::new(x_cord, z_cord, self.a_24.clone(), self.modulus.clone())
    }

    /// Doubles a point in an elliptic curve in Montgomery form.
    pub fn double(&self) -> Point {
        let u = Integer::from(&self.x_cord + &self.z_cord).square();
        let v = Integer::from(&self.x_cord - &self.z_cord).square();
        let diff = Integer::from(&u - &v);
        let x_cord = (u * &v) % &self.modulus;
        let z_cord = ((v + &self.a_24 * &diff) * diff) % &self.modulus;

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
    pub fn mont_ladder(&self, k: &Integer) -> Point {
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
            self.z_cord.clone().invert(&self.modulus).unwrap() * &self.x_cord % &self.modulus
                == other.z_cord.clone().invert(&self.modulus).unwrap() * &other.x_cord
                    % &self.modulus
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::Integer;

    #[test]
    fn test_point_add() {
        let p1 = Point::new(11.into(), 16.into(), 7.into(), 29.into());
        let p2 = Point::new(13.into(), 10.into(), 7.into(), 29.into());
        let p3 = p2.add(&p1, &p1);

        assert_eq!(p3.x_cord, Integer::from(23));
        assert_eq!(p3.z_cord, Integer::from(17));
    }

    #[test]
    fn test_point_double() {
        let p1 = Point::new(11.into(), 16.into(), 7.into(), 29.into());
        let p2 = p1.double();

        assert_eq!(p2.x_cord, Integer::from(13));
        assert_eq!(p2.z_cord, Integer::from(10));
    }

    #[test]
    fn test_point_mont_ladder() {
        let p1 = Point::new(11.into(), 16.into(), 7.into(), 29.into());
        let p3 = p1.mont_ladder(&3.into());

        assert_eq!(p3.x_cord, Integer::from(23));
        assert_eq!(p3.z_cord, Integer::from(17));
    }

    #[test]
    fn test_point() {
        let modulus = 101.into();
        let a: Integer = 10.into();
        let a_24: Integer = (a + Integer::from(2)) * Integer::from(4).invert(&modulus).unwrap();

        let p1 = Point::new(10.into(), 17.into(), a_24.clone(), modulus.clone());
        let p2 = p1.double();
        assert_eq!(
            p2,
            Point::new(68.into(), 56.into(), a_24.clone(), modulus.clone())
        );
        let p4 = p2.double();
        assert_eq!(
            p4,
            Point::new(22.into(), 64.into(), a_24.clone(), modulus.clone())
        );
        let p8 = p4.double();
        assert_eq!(
            p8,
            Point::new(71.into(), 95.into(), a_24.clone(), modulus.clone())
        );
        let p16 = p8.double();
        assert_eq!(
            p16,
            Point::new(5.into(), 16.into(), a_24.clone(), modulus.clone())
        );
        let p32 = p16.double();
        assert_eq!(
            p32,
            Point::new(33.into(), 96.into(), a_24.clone(), modulus.clone())
        );

        // p3 = p2 + p1
        let p3 = p2.add(&p1, &p1);
        assert_eq!(
            p3,
            Point::new(1.into(), 61.into(), a_24.clone(), modulus.clone())
        );
        // p5 = p3 + p2 or p4 + p1
        let p5 = p3.add(&p2, &p1);
        assert_eq!(
            p5,
            Point::new(49.into(), 90.into(), a_24.clone(), modulus.clone())
        );
        assert_eq!(p5, p4.add(&p1, &p3));
        // # p6 = 2*p3
        let p6 = p3.double();
        assert_eq!(
            p6,
            Point::new(87.into(), 43.into(), a_24.clone(), modulus.clone())
        );
        assert_eq!(p6, p4.add(&p2, &p2));
        // # p7 = p5 + p2
        let p7 = p5.add(&p2, &p3);
        assert_eq!(
            p7,
            Point::new(69.into(), 23.into(), a_24.clone(), modulus.clone())
        );
        assert_eq!(p7, p4.add(&p3, &p1));
        assert_eq!(p7, p6.add(&p1, &p5));
        // # p9 = p5 + p4
        let p9 = p5.add(&p4, &p1);
        assert_eq!(p9, Point::new(56.into(), 99.into(), a_24, modulus));
        assert_eq!(p9, p6.add(&p3, &p3));
        assert_eq!(p9, p7.add(&p2, &p5));
        assert_eq!(p9, p8.add(&p1, &p7));

        assert_eq!(p5, p1.mont_ladder(&5.into()));
        assert_eq!(p9, p1.mont_ladder(&9.into()));
        assert_eq!(p16, p1.mont_ladder(&16.into()));
        assert_eq!(p9, p3.mont_ladder(&3.into()));
    }
}
