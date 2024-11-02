/// MaxValue returns the maximum value of the Hilbert number for a Hilbert curve of dimensionality n and order m.
/// This will be:
///	(2 ** (n * m)) - 1
pub fn max_value(n: u8, m: u8) -> u32 {
    (1u32 << (n * m)) - 1
}

/// MaxCoord returns the maximum coordinate on each axis for a Hilbert curve of order m.
/// This will be:
///	(2 ** m) - 1
pub fn max_coord(m: u8) -> u32 {
    return (1u32 << m) - 1;
}

/// TODO: Fix this documentation to be more rustacean.
///
/// Convert a Hilbert value for a point in n-space
///
/// This basically uses the algorithm described in
///
/// "Alternative Algorithm for Hillbert's Space-Filling Curve," A R Butz,
/// IEEE Trans. on Computers, April 1971, pp 424-426.
///
/// For dimensionality n and order m:
/// - The Hilbert number varies from 0 to H = (2 ** (n * m)) - 1
/// - The integer grid varies from 0 to G = (2 ** m) - 1 along each axis
/// - Examples:
/// - n = 2, m = 4 -> H = 255, G = 15
/// - n = 3, m = 5 -> H = 32767, G = 31
///
/// @param n dimensionality of the space.
/// @param m order of the Hilbert curve.
/// @param h Hilbert value to be converted.
/// @return a long array of length n holding the coordinates of the point.
///
pub fn hilbert_value_to_point(n: u8, m: u8, h: u32) -> Vec<u32> {
    // Make an mask with n bits that will be used to chop up the
    // hilbert value and to perform rotations.
    let mask = (1u32 << n) - 1;

    // This will contain h split up into m chunks of n bits.
    let mut rho = vec![0u32; m as usize];

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
    let mut j = vec![0u8; m as usize];

    // Each sigma is an n-bit value such that
    //
    //	sigma[i](0) = rho[i](0)
    //	sigma[i](1) = rho[i](1) ^ rho[i](0)
    //	sigma[i](2) = rho[i](2) ^ rho[i](1)
    //	...
    //	sigma[i](n-1) = rho[i](n-1) ^ rho[i](n-2)
    //
    // where the notation a[b](c) represents the c'th bit of a[b].
    let mut sigma = vec![0u32; m as usize];

    // Each tao is an n-bit value obtained by complementing the
    // corresponding sigma in the (n-1)th position and then, if and
    // only if the resulting value is of odd parity, complementing in
    // the principal position.  Hence, each tao is always of even
    // parity.  Note that the parity of sigma[i] is given by the bit
    // rho[i](n-1).
    let mut tao = vec![0u32; m as usize];

    // This is sigma~ in the Butz paper.  Each delta[i] is an n-bit
    // value obtained by rotating right sigma[i] by a number of
    // positions equal to j[0] + j[1] + ... + j[i-1].  We keep track of
    // this quantity in rotation.
    let mut delta = vec![0u32; m as usize];

    // This is tao~ in the Butz paper.  Each gamma[i] is an n-bit value
    // obtained by rotating right tao[i] by a number of positions equal
    // to j[0] + j[1] + ... + j[i-1].  We keep track of this quantity
    // in rotation.
    let mut gamma = vec![0u32; m as usize];

    // Each omega is an n-bit value where omega[0] = 0 and
    // omega[i] = omega[i-1] ^ gamma[i-1]
    let mut omega = vec![0u32; m as usize];

    // Each alpha is an n-bit value where alpha[i] = omega[i] ^ delta[i]
    let mut alpha = vec![0u32; m as usize];

    // An array of length n of m-bit value obtained by bit transposing
    // the alpha array of length m of n-bit values.
    let mut a = vec![0u32; n as usize];

    // This is a running sum of j[i] used to determine how much to
    // rotate right sigma and tao to produce delta and gamma,
    // respectively.
    let mut rotation = 0u8;

    for i in 0..m as usize {
        // Break up the m*n bits of the hilbert value into an
        // array of m n-bit values.
        let tmp = (h >> ((m - i as u8 - 1) * n)) & mask;
        rho[i] = tmp;

        // Calculate sigma from rho
        sigma[i] = tmp ^ (tmp >> 1);

        // Calculate principal position of rho
        let parity = tmp & 1;
        j[i] = n - 1;
        for k in 1..n {
            if ((tmp >> k) & 1) != parity {
                j[i] = n - k - 1;
                break;
            }
        }

        // Calculate tao
        tao[i] = sigma[i] ^ 1;
        if parity == 0 {
            tao[i] ^= 1 << (n - j[i] - 1);
        }

        // Calculate delta as a right rotation of sigma
        delta[i] = ((sigma[i] << (n - rotation)) | (sigma[i] >> rotation)) & mask;

        // Calculate gamma as a right rotation of tao
        gamma[i] = ((tao[i] << (n - rotation)) | (tao[i] >> rotation)) & mask;

        rotation += j[i];
        if rotation > n {
            rotation -= n;
        }
    }

    // Calculate omegas and alphas
    omega[0] = 0;
    alpha[0] = delta[0];
    for i in 1..m as usize {
        omega[i] = omega[i - 1] ^ gamma[i - 1];
        alpha[i] = omega[i] ^ delta[i];
    }

    // Calculate the result
    for k in 0..n as usize {
        a[k] = 0;
        let shift = n - k as u8 - 1;
        for i in 0..m as usize {
            a[k] |= ((alpha[i] >> shift) & 1) << (m - i as u8 - 1);
        }
    }

    // #    print "h", h
    // #    _print_vec("i", range(m), m)
    // #    _print_vec("rho", rho, m)
    // #    _print_vec("j", j, m)
    // #    _print_vec("sigma", sigma, m)
    // #    _print_vec("tao", tao, m)
    // #    _print_vec("delta", delta, m)
    // #    _print_vec("gamma", gamma, m)
    // #    _print_vec("omega", omega, m)
    // #    _print_vec("alpha", alpha, m)

    a
}

/// PointToValue converts a point in n-space a on a Hilbert Curve of
/// dimensionality n and order m hilbert number.
pub fn point_to_hilbert_value(m: u8, a: Vec<u32>) -> u32 {
    // Get the dimensionality of the vector
    let n = a.len() as u8;

    // Make mask for relevant bits
    let mask = (1u32 << n) - 1;

    let mut rho = vec![0u32; m as usize];
    let mut j = vec![0u8; m as usize];
    let mut sigma = vec![0u32; m as usize];
    let mut tao = vec![0u32; m as usize];
    let mut delta = vec![0u32; m as usize];
    let mut gamma = vec![0u32; m as usize];
    let mut omega = vec![0u32; m as usize];
    let mut alpha = vec![0u32; m as usize];

    // Compute the alpha values from a[]
    for i in 0..m as usize {
        alpha[i] = 0;
        let shift = m - i as u8 - 1;
        for k in 0..n as usize {
            alpha[i] |= ((a[k] >> shift) & 1) << (n - k as u8 - 1);
        }
    }

    // Calculate values for i == 0 in preperation for loop
    omega[0] = 0;
    delta[0] = alpha[0];
    sigma[0] = delta[0];

    // This computes rho[0] from sigma[0] as the inverse of
    // sigma[0] = rho[0] ^ (rho[0] >>> 1)
    let mut bit = 1u32 << (n - 1);
    rho[0] = sigma[0] & bit;
    for k in (0..n - 1).rev() {
        bit >>= 1;
        rho[0] |= (sigma[0] ^ (rho[0] >> 1)) & bit;
    }

    let mut rotation = 0u8;

    let mut r = rho[0];

    for i in 1..m as usize {
        let p = i - 1;

        // Calculate principal position of previous rho
        j[p] = n - 1;
        let tmp = rho[p];
        let parity = tmp & 1;
        for k in 1..n {
            if ((tmp >> k) & 1) != parity {
                j[p] = n - k - 1;
                break;
            }
        }

        // Calculate previous tao
        tao[p] = sigma[p] ^ 1;
        if parity == 0 {
            tao[p] ^= 1 << (n - j[p] - 1);
        }

        // Calculate gamma as a right rotation of tao
        gamma[p] = ((tao[p] >> rotation) | (tao[p] << (n - rotation))) & mask;

        // Calculate omega
        omega[i] = omega[p] ^ gamma[p];

        // Calculate delta
        delta[i] = alpha[i] ^ omega[i];

        rotation += j[p];
        if rotation > n {
            rotation -= n
        }

        // Calculate sigma as a left rotation of delta
        sigma[i] = ((delta[i] << rotation) | (delta[i] >> (n - rotation))) & mask;

        // This computes rho[i] from sigma[i] as the inverse of
        // sigma[i] = rho[i] ^ (rho[i] >>> 1)
        bit = 1 << (n - 1);
        rho[i] = sigma[i] & bit;
        for _ in (0..n - 1).rev() {
            bit >>= 1;
            rho[i] |= (sigma[i] ^ (rho[i] >> 1)) & bit;
        }

        r = (r << n) | rho[i]
    }

    //fmt.Printf("mask %v\n", mask)
    //fmt.Printf("alpha %v\n", alpha)
    //fmt.Printf("omega %v\n", omega)
    //fmt.Printf("delta %v\n", delta)
    //fmt.Printf("sigma %v\n", sigma)
    //fmt.Printf("rho %v\n", rho)
    //fmt.Printf("tao %v\n", tao)
    //fmt.Printf("gamma %v\n", gamma)

    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_max_value() {
        assert_eq!(max_value(2, 4), 255);
        assert_eq!(max_value(3, 5), 32767);
    }

    #[test]
    fn test_max_coord() {
        assert_eq!(max_coord(4), 15);
        assert_eq!(max_coord(5), 31);
    }

    #[test]
    fn test_conversions() {
        // Vary the dimensionality
        for n in 1..4 {
            // Vary the order
            for m in 1..7 {
                let max_value = max_value(n, m);
                let max_coord = max_coord(m);

                println!("dimensionality n = {}, order m = {}, max_coord = {}, max_value = {}", n, m, max_coord, max_value);

			// Allocate storage for all coordinates in hilbert value order and as a set
                let mut coords = vec![Vec::new(); (max_value + 1) as usize];
                let mut set = std::collections::HashSet::new();

                for h in 0..=max_value {
                    coords[h as usize] = hilbert_value_to_point(n, m, h);
                    for &c in &coords[h as usize] {
                        assert!(c <= max_coord, "Coordinates {:?} > {}", coords[h as usize], max_coord);
                    }
                    set.insert(coords[h as usize].clone());
                    let h_test = point_to_hilbert_value(m, coords[h as usize].clone());
                    assert_eq!(h, h_test, "Hilbert value mismatch {} -> {:?} -> {}", h, coords[h as usize], h_test);
                }

                assert_eq!(coords.len(), set.len(), "Not all {} coordinates are unique", max_value + 1);

                for i in 1..coords.len() {
                    let mut change_count = 0;
                    let a = &coords[i];
                    let b = &coords[i - 1];
                    for j in 0..n as usize {
                        if a[j] != b[j] {
                            change_count += 1;
                            let diff = if a[j] > b[j] { a[j] - b[j] } else { b[j] - a[j] };
                            assert_eq!(diff, 1, "{:?} and {:?} differ by {}", b, a, diff);
                        }
                    }
                    assert_eq!(change_count, 1, "Too many coordinates change between {:?} and {:?}", a, b);
                }
            }
        }
    }
}

