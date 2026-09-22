#[cfg(any(test, feature = "bench"))]
use rug::Assign;
use rug::Integer;
use std::rc::Rc;

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
/// `a_24` and the modulus are shared by all the points of a curve (cloning a point doesn't copy
/// them).
///
/// This is the interface type of the stages, and its operations are a simple reference
/// implementation (for tests and benchmarks only): the stages compute in Montgomery
/// representation, see `crate::curve`.
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
    pub a_24: Rc<Integer>,
    /// modulus
    pub modulus: Rc<Integer>,
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
            a_24: Rc::new(a_24),
            modulus: Rc::new(modulus),
        }
    }

    /// Point with the given coordinates, on the same curve as `self`.
    pub(crate) fn on_same_curve(&self, x_cord: Integer, z_cord: Integer) -> Point {
        Point {
            x_cord,
            z_cord,
            a_24: Rc::clone(&self.a_24),
            modulus: Rc::clone(&self.modulus),
        }
    }

    /// Empty integer large enough for the intermediate results of `add` and `double` (up to 5
    /// times the size of the modulus, before reduction): computing in place in it never needs
    /// to reallocate.
    #[cfg(any(test, feature = "bench"))]
    fn scratch(&self) -> Integer {
        Integer::with_capacity(5 * self.modulus.significant_bits() as usize + 64)
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
    #[cfg(any(test, feature = "bench"))]
    pub fn add(&self, q: &Point, diff: &Point) -> Point {
        let n: &Integer = &self.modulus;
        let (mut u, mut v, mut t) = (self.scratch(), self.scratch(), self.scratch());

        // u = (x1 - z1) * (x2 + z2)
        u.assign(&self.x_cord - &self.z_cord);
        t.assign(&q.x_cord + &q.z_cord);
        u *= &t;
        // v = (x1 + z1) * (x2 - z2)
        v.assign(&self.x_cord + &self.z_cord);
        t.assign(&q.x_cord - &q.z_cord);
        v *= &t;

        // x = diff.z * (u + v)^2
        t.assign(&u + &v);
        t.square_mut();
        t *= &diff.z_cord;
        t %= n;
        // z = diff.x * (u - v)^2
        u -= &v;
        u.square_mut();
        u *= &diff.x_cord;
        u %= n;

        self.on_same_curve(t, u)
    }

    /// Doubles a point in an elliptic curve in Montgomery form.
    #[cfg(any(test, feature = "bench"))]
    pub fn double(&self) -> Point {
        let n: &Integer = &self.modulus;
        let (mut u, mut v, mut diff) = (self.scratch(), self.scratch(), self.scratch());

        // u = (x + z)^2, v = (x - z)^2
        u.assign(&self.x_cord + &self.z_cord);
        u.square_mut();
        u %= n;
        v.assign(&self.x_cord - &self.z_cord);
        v.square_mut();
        v %= n;
        diff.assign(&u - &v);

        // x = u * v
        u *= &v;
        u %= n;
        // z = (v + a_24 * diff) * diff
        v += &*self.a_24 * &diff;
        v *= &diff;
        v %= n;

        self.on_same_curve(u, v)
    }

    /// Scalar multiplication of a point in Montgomery form
    /// using Montgomery Ladder Algorithm.
    /// A total of 11 multiplications are required in each step of this
    /// algorithm.
    ///
    /// # Parameters
    ///
    /// - `k`: The positive integer multiplier
    #[cfg(any(test, feature = "bench"))]
    pub fn mont_ladder(&self, k: &Integer) -> Point {
        let mut q = self.clone();
        let mut r = self.double();

        // Bits of `k` from the most significant one, which is skipped.
        for bit in (0..k.significant_bits().saturating_sub(1)).rev() {
            if k.get_bit(bit) {
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
        // X1/Z1 = X2/Z2 without any modular inverse: X1*Z2 = X2*Z1.
        self.a_24 == other.a_24
            && self.modulus == other.modulus
            && Integer::from(&self.x_cord * &other.z_cord) % &*self.modulus
                == Integer::from(&other.x_cord * &self.z_cord) % &*self.modulus
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
