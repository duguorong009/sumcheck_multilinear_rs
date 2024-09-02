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

#[test]
fn test_mod_inverse() {
    assert_eq!(mod_inverse(3, 11), 4);
}
