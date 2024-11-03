//! Package rhilbert provides functions to convert back and forth between
//! Hilbert Curve values and points in n-space.
//!
//! This uses an extension to the algorithm described in

//! "Alternative Algorithm for Hillbert's Space-Filling Curve," A R
//! Butz, IEEE Trans. on Computers, April 1971, pp 424-426. (NOTE: Butz
//! was a rabid Holocaust denier but this does not mean we can't use this work.)
//!
//! The extension added here allows for conversion back and forth and not just
//! in one direction.
//!
//! For dimensionality n and order m:
//! - The Hilbert number varies from 0 to H = (2 ** (n * m)) - 1
//! - The integer grid varies from 0 to G = (2 ** m) - 1 along each axis
//! TODO - Examples:
//! - n = 2, m = 4 -> H = 255, G = 15
//! - n = 3, m = 5 -> H = 32767, G = 31
//!
//! For n = 2, m = 2 there are H=15 unique Hilbert numbers and the X
//! and Y coordinates vary between 0 and G=3.
//!
//! The correspondence between Hilbert numbers and coordinates is shown
//! here for n = 2 and m = 2:
//! ```text
//!     3| 5---6   9--10
//!      | |   |   |   |
//!     2| 4   7---8  11
//!    Y | |           |
//!     1| 3---2  13--12
//!      |     |   |
//!     0| 0---1  14--15
//!      +--------------
//!        0   1   2   3
//!            X
//! ```
//! ```text
//!    Hilbert    Coordinate
//!    Number      [X Y]
//!    -------    ----------
//!      0   <->   [0 0]
//!      1   <->   [1 0]
//!      2   <->   [1 1]
//!      3   <->   [0 1]
//!      4   <->   [0 2]
//!      5   <->   [0 3]
//!      6   <->   [1 3]
//!      7   <->   [1 2]
//!      8   <->   [2 2]
//!      9   <->   [2 3]
//!     10   <->   [3 3]
//!     11   <->   [3 2]
//!     12   <->   [3 1]
//!     13   <->   [2 1]
//!     14   <->   [2 0]
//!     15   <->   [3 0]
//! ```

pub struct HilbertCurve {
    n: u8,
    m: u8,
    mask: u32,
    rho: Vec<u32>,
    j: Vec<u8>,
    sigma: Vec<u32>,
    tao: Vec<u32>,
    delta: Vec<u32>,
    gamma: Vec<u32>,
    omega: Vec<u32>,
    alpha: Vec<u32>,
}

impl HilbertCurve {
    /// Create a new HilbertCurve with dimensionality n and order m.
    pub fn new(n: u8, m: u8) -> Self {
        HilbertCurve {
            // Dimensionality of the space
            n,
            // Order of the Hilbert curve
            m,
            // Make an mask with n bits that will be used to split/recombine the
            // hilbert value and to perform rotations.
            mask: (1u32 << n) - 1,
            // Scratch arrays derived from the Butz paper. No need to reinitialize
            // these between calls since the algorithms initialize exhaustively.
            rho: vec![0u32; m as usize],
            j: vec![0u8; m as usize],
            sigma: vec![0u32; m as usize],
            tao: vec![0u32; m as usize],
            delta: vec![0u32; m as usize],
            gamma: vec![0u32; m as usize],
            omega: vec![0u32; m as usize],
            alpha: vec![0u32; m as usize],
        }
    }

    /// TODO - Fix this documentation to be more rustacean.
    ///  MaxValue returns the maximum value of the Hilbert number for a Hilbert curve of dimensionality n and order m.
    /// This will be:
    /// (2 ** (n * m)) - 1
    pub fn max_value(&self) -> u32 {
        (1u32 << (self.n * self.m)) - 1
    }

    /// MaxCoord returns the maximum coordinate on each axis for a Hilbert curve of order m.
    /// This will be:
    /// (2 ** m) - 1
    pub fn max_coord(&self) -> u32 {
        (1u32 << self.m) - 1
    }

    /// TODO: Fix this documentation to be more rustacean.
    ///
    /// Convert a Hilbert value, h, to a point in n-space
    ///
    /// TODO - Examples:
    /// - n = 2, m = 4 -> H = 255, G = 15
    /// - n = 3, m = 5 -> H = 32767, G = 31
    pub fn value_to_point(&mut self, h: u32) -> Vec<u32> {
        // self.reset();

        // This is a running sum of j[i] used to determine how much to
        // rotate right sigma and tao to produce delta and gamma,
        // respectively.
        let mut rotation = 0u8;

        for i in 0..self.m as usize {
            // Break up the m*n bits of the hilbert value into an
            // array of m n-bit values.
            // rho will contain h split up into m chunks of n bits.
            let tmp = (h >> ((self.m - i as u8 - 1) * self.n)) & self.mask;
            self.rho[i] = tmp;

            // Each sigma is an n-bit value such that
            //
            //	sigma[i](0) = rho[i](0)
            //	sigma[i](1) = rho[i](1) ^ rho[i](0)
            //	sigma[i](2) = rho[i](2) ^ rho[i](1)
            //	...
            //	sigma[i](n-1) = rho[i](n-1) ^ rho[i](n-2)
            //
            // where the notation a[b](c) represents the c'th bit of a[b].
            self.sigma[i] = tmp ^ (tmp >> 1);

            // Each j is an integer between 0 and n-1 equal to thesubscript of
            // the principal position for rho[i].  The principal position is
            // the last position bit position in rho[i] such that bit[j] !=
            // bit[n-1].  If all bits are equal, then j == n - 1.  In the Butz
            // paper these values are between 1 and n but starting at 0
            // simplifies things.  For example, the principal positions are
            // marked with a * below and the resulting value of j is shown on
            // the right...
            //
            //	1	1	1	1	1*	-> 4
            //	1	0*	1	1	1	-> 1
            //	0	0	1	1*	0	-> 3
            //	0	0	0	0	0*	-> 4
            //
            let parity = tmp & 1;
            self.j[i] = self.n - 1;
            for k in 1..self.n {
                if ((tmp >> k) & 1) != parity {
                    self.j[i] = self.n - k - 1;
                    break;
                }
            }

            // Each tao is an n-bit value obtained by complementing the
            // corresponding sigma in the (n-1)th position and then, if and
            // only if the resulting value is of odd parity, complementing in
            // the principal position.  Hence, each tao is always of even
            // parity.  Note that the parity of sigma[i] is given by the bit
            // rho[i](n-1).
            self.tao[i] = self.sigma[i] ^ 1;
            if parity == 0 {
                self.tao[i] ^= 1 << (self.n - self.j[i] - 1);
            }

            // This is sigma~ in the Butz paper.  Each delta[i] is an n-bit
            // value obtained by rotating right sigma[i] by a number of
            // positions equal to j[0] + j[1] + ... + j[i-1].  We keep track of
            // this quantity in rotation.
            self.delta[i] =
                ((self.sigma[i] << (self.n - rotation)) | (self.sigma[i] >> rotation)) & self.mask;

            // This is tao~ in the Butz paper.  Each gamma[i] is an n-bit value
            // obtained by rotating right tao[i] by a number of positions equal
            // to j[0] + j[1] + ... + j[i-1].  We keep track of this quantity
            // in rotation.
            self.gamma[i] =
                ((self.tao[i] << (self.n - rotation)) | (self.tao[i] >> rotation)) & self.mask;

            rotation += self.j[i];
            if rotation > self.n {
                rotation -= self.n;
            }
        }

        // Each omega is an n-bit value where omega[0] = 0 and
        // omega[i] = omega[i-1] ^ gamma[i-1]
        self.omega[0] = 0;
        // Each alpha is an n-bit value where alpha[i] = omega[i] ^ delta[i]
        self.alpha[0] = self.delta[0];
        for i in 1..self.m as usize {
            self.omega[i] = self.omega[i - 1] ^ self.gamma[i - 1];
            self.alpha[i] = self.omega[i] ^ self.delta[i];
        }

        // The result, a, is an array of length n of m-bit value obtained by bit transposing
        // the alpha array of length m of n-bit values.
        let mut a = vec![0u32; self.n as usize];
        for (k, a_k) in a.iter_mut().enumerate() {
            *a_k = (0..self.m as usize).fold(0, |acc, i| {
                acc | (((self.alpha[i] >> (self.n - k as u8 - 1)) & 1) << (self.m - i as u8 - 1))
            });
        }
        a
    }

    /// Converts a point in n-space a on a Hilbert Curve of
    /// dimensionality n and order m hilbert number. The dimensionality of the
    /// input vector, a, is used even if it doesn't match the dimensionality of
    /// the HilbertCurve.
    /// TODO - implement try_point_to_value to enforce the correct dimensionality.
    pub fn point_to_value(&mut self, a: Vec<u32>) -> u32 {
        // Get the dimensionality of the vector
        let n = a.len() as u8;

        // Make mask for relevant bits
        let mask = (1u32 << n) - 1;

        // Compute the alpha values from a[]
        for (i, alpha_i) in self.alpha.iter_mut().enumerate() {
            *alpha_i = (0..n as usize).fold(0, |acc, k| {
                acc | (((a[k] >> (self.m - i as u8 - 1)) & 1) << (n - k as u8 - 1))
            });
        }

        // Calculate values for i == 0 in preperation for loop
        self.omega[0] = 0;
        self.delta[0] = self.alpha[0];
        self.sigma[0] = self.delta[0];

        // This computes rho[0] from sigma[0] as the inverse of
        // sigma[0] = rho[0] ^ (rho[0] >>> 1)
        let mut bit = 1u32 << (n - 1);
        self.rho[0] = self.sigma[0] & bit;
        for _ in 0..n - 1 {
            bit >>= 1;
            self.rho[0] |= (self.sigma[0] ^ (self.rho[0] >> 1)) & bit;
        }

        let mut rotation = 0u8;

        let mut r = self.rho[0];

        for i in 1..self.m as usize {
            let p = i - 1;

            // Calculate principal position of previous rho
            self.j[p] = n - 1;
            let tmp = self.rho[p];
            let parity = tmp & 1;
            for k in 1..n {
                if ((tmp >> k) & 1) != parity {
                    self.j[p] = n - k - 1;
                    break;
                }
            }

            // Calculate previous tao
            self.tao[p] = self.sigma[p] ^ 1;
            if parity == 0 {
                self.tao[p] ^= 1 << (n - self.j[p] - 1);
            }

            // Calculate gamma as a right rotation of tao
            self.gamma[p] = ((self.tao[p] >> rotation) | (self.tao[p] << (n - rotation))) & mask;

            // Calculate omega
            self.omega[i] = self.omega[p] ^ self.gamma[p];

            // Calculate delta
            self.delta[i] = self.alpha[i] ^ self.omega[i];

            rotation += self.j[p];
            if rotation > n {
                rotation -= n
            }

            // Calculate sigma as a left rotation of delta
            self.sigma[i] =
                ((self.delta[i] << rotation) | (self.delta[i] >> (n - rotation))) & mask;

            // This computes rho[i] from sigma[i] as the inverse of
            // sigma[i] = rho[i] ^ (rho[i] >>> 1)
            bit = 1 << (n - 1);
            self.rho[i] = self.sigma[i] & bit;
            for _ in (0..n - 1).rev() {
                bit >>= 1;
                self.rho[i] |= (self.sigma[i] ^ (self.rho[i] >> 1)) & bit;
            }

            r = (r << n) | self.rho[i]
        }
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_max_value() {
        let hc = HilbertCurve::new(2, 4);
        assert_eq!(hc.max_value(), 255);
        let hc = HilbertCurve::new(3, 5);
        assert_eq!(hc.max_value(), 32767);
    }

    #[test]
    fn test_max_coord() {
        let hc = HilbertCurve::new(2, 4);
        assert_eq!(hc.max_coord(), 15);
        let hc = HilbertCurve::new(3, 5);
        assert_eq!(hc.max_coord(), 31);
    }

    #[test]
    fn test_conversions() {
        // Vary the dimensionality
        for n in 1..4 {
            // Vary the order
            for m in 1..7 {
                let mut hc = HilbertCurve::new(n, m);

                let max_value = hc.max_value();
                let max_coord = hc.max_coord();

                // println!(
                //     "dimensionality n = {}, order m = {}, max_coord = {}, max_value = {}",
                //     n, m, max_coord, max_value
                // );

                // Allocate storage for all coordinates in hilbert value order and as a set
                let mut coords = vec![Vec::new(); (max_value + 1) as usize];
                let mut set = std::collections::HashSet::new();

                for h in 0..=max_value {
                    coords[h as usize] = hc.value_to_point(h);
                    for &c in &coords[h as usize] {
                        assert!(
                            c <= max_coord,
                            "Coordinates {:?} > {}",
                            coords[h as usize],
                            max_coord
                        );
                    }
                    set.insert(coords[h as usize].clone());
                    let h_test = hc.point_to_value(coords[h as usize].clone());
                    assert_eq!(
                        h, h_test,
                        "Hilbert value mismatch {} -> {:?} -> {}",
                        h, coords[h as usize], h_test
                    );
                }

                assert_eq!(
                    coords.len(),
                    set.len(),
                    "Not all {} coordinates are unique",
                    max_value + 1
                );

                for i in 1..coords.len() {
                    let mut change_count = 0;
                    let a = &coords[i];
                    let b = &coords[i - 1];
                    for j in 0..n as usize {
                        if a[j] != b[j] {
                            change_count += 1;
                            let diff = if a[j] > b[j] {
                                a[j] - b[j]
                            } else {
                                b[j] - a[j]
                            };
                            assert_eq!(diff, 1, "{:?} and {:?} differ by {}", b, a, diff);
                        }
                    }
                    assert_eq!(
                        change_count, 1,
                        "Too many coordinates change between {:?} and {:?}",
                        a, b
                    );
                }
            }
        }
    }
}
