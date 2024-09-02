use rand::{
    rngs::{StdRng, ThreadRng},
    Rng, SeedableRng,
};

use crate::pmf::PMF;

pub const MAX_ALLOWED_SOUNDNESS_ERROR: f64 = 2e-64;

trait RandomGen {
    fn get_random_element(&mut self) -> u64 {
        unimplemented!("get_random_element")
    }
}

#[derive(Debug)]
pub struct TrueRandomGen {
    rng: StdRng,
    p: u64,
}

impl TrueRandomGen {
    fn new(seed: u64, p: u64) -> Self {
        let rng = rand::rngs::StdRng::seed_from_u64(seed);
        Self { rng, p }
    }
}

impl RandomGen for TrueRandomGen {
    fn get_random_element(&mut self) -> u64 {
        self.rng.gen_range(0..self.p)
    }
}

/// An interactive verifier that verifies the sum of the polynomial which is the product of multilinear functions
#[derive(Debug)]
pub struct InteractivePMFVerifier {
    checksum_only: bool,
    p: u64,
    poly: PMF,
    asserted_sum: u64,
    rng: TrueRandomGen,
    active: bool,
    convinced: bool,
    points: Vec<u64>,
    round: usize,
    expect: u64,
}

impl InteractivePMFVerifier {
    pub fn new(
        poly: PMF,
        asserted_sum: u64,
        max_allowed_soundness_err: Option<f64>,
        checksum_only: Option<bool>,
        rng: Option<TrueRandomGen>,
    ) -> Self {
        let poly_num_variables = poly.num_variables;
        let p = poly.p;
        let seed = rand::thread_rng().gen_range(0..0xFFFFFFFF);

        // check soundness
        let err = {
            let deg = poly.num_variables * poly.num_multiplicands();
            (poly.num_variables * deg) as f64 / poly.p as f64
        };
        if err > max_allowed_soundness_err.unwrap_or(MAX_ALLOWED_SOUNDNESS_ERROR) {
            panic!("Soundness error {} exceeds maximum allowed soundness error {}. Try to have a prime with size >= ??? bits", err, max_allowed_soundness_err.unwrap_or(MAX_ALLOWED_SOUNDNESS_ERROR));
        }

        // some edge cases: if univariate or constant: no need to be interactive
        if poly.num_variables == 0 {
            if (asserted_sum - poly.eval(&[])) % poly.p == 0 {
                todo!("implement `_convince_and_close` & return")
            } else {
                todo!("implement `_reject_and_close` & return")
            }
        }
        if poly.num_variables == 1 {
            let result = (poly.eval(&[0]) + poly.eval(&[1])) % poly.p;
            if (asserted_sum - result) % poly.p == 0  {
                todo!("implement `_convince_and_close` & return")
            } else {
                todo!("implement `_reject_and_close` & return")
            }
        }

        Self {
            checksum_only: checksum_only.unwrap_or(false),
            p,
            poly,
            asserted_sum: asserted_sum % p,
            rng: rng.unwrap_or_else(|| TrueRandomGen::new(seed, p)),
            active: true,
            convinced: false,
            points: vec![0; poly_num_variables],
            round: 0,
            expect: asserted_sum,
        }
    }
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
