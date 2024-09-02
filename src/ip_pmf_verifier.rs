use rand::rngs::ThreadRng;

use crate::pmf::PMF;

#[derive(Debug)]
pub struct InteractivePMFVerifier {
    checksum_only: bool,
    p: u64,
    poly: PMF,
    asserted_sum: u64,
    rng: ThreadRng,
    active: bool,
    convinced: bool,
    points: Vec<u64>,
    round: usize,
    expect: u64,
}

// modular inverse (https://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/)
//     :param a: a
//     :param m: prime
//     :return: a^-1 mod p
fn mod_inverse(a: u64, m: u64) -> u64 {
    let m0: i128 = m.into();
    let mut y: i128 = 0;
    let mut x: i128 = 1;

    if m == 1 {
        return 0;
    }

    let mut m: i128 = m.into();
    let mut a: i128 = a.into();
    while a > 1 {
        let q = a / m;
        let mut t: i128 = m;
        m = a % m;
        a = t;
        t = y;
        y = x - q * y;
        x = t;
    }
    if x < 0 {
        x += m0;
    }
    x.try_into().unwrap()
}

///Interpolate and evaluate a PMF.
// Adapted from https://www.geeksforgeeks.org/lagranges-interpolation/
// :param points: P_0, P_1, P_2, ..., P_m where P is the product of m multilinear polynomials
// :param r: The point we want to evaluate at. In this scenario, the verifier wants to evaluate p_r.
// :param p: Field size.
// :return: P_r
fn interpolate(points: &[u64], r: u64, p: u64) -> u64 {
    let mut result = 0;
    for i in 0..points.len() {
        let mut term = points[i] % p;
        for j in 0..points.len() {
            if j != i {
                // TODO: `(r - j) % p` & `(i - j) % p` could be not always in the range [0, p-1]
                let r_min_j = if r >= j as u64 {
                    r - j as u64
                } else {
                    r + p - j as u64
                };
                let i_min_j = if i >= j {
                    (i - j) as u64
                } else {
                    i as u64 + p - j as u64
                };
                term = (term * (r_min_j % p) * mod_inverse(i_min_j % p, p)) % p;
            }
        }
        result = (result + term) % p;
    }
    result
}

#[test]
fn test_mod_inverse() {
    assert_eq!(mod_inverse(3, 11), 4);
}
