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

#[derive(Debug, Clone)]
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
    pub(crate) active: bool,
    pub(crate) convinced: bool,
    pub(crate) points: Vec<u64>,
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
            if (asserted_sum - result) % poly.p == 0 {
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

    fn random_r(&mut self) -> u64 {
        self.rng.get_random_element()
    }

    fn soundness_error(&self) -> f64 {
        let deg = self.poly.num_variables * self.poly.num_multiplicands();
        (self.poly.num_variables * deg) as f64 / self.p as f64
    }

    /// The minimum size of prime required to meet the soundness error constraint.
    fn required_field_length_bit(&self, e: f64) -> usize {
        let deg = self.poly.num_variables * self.poly.num_multiplicands();
        let min_p = (self.poly.num_variables as f64 * deg as f64) / e;
        min_p.log2().ceil() as usize
    }

    /// Send this verifier the univariate polynomial P(x). P(x) has degree at most the number of multiplicands.
    /// :param msgs: [P(0), P(1), ..., P(m)] where m is the number of multiplicands
    /// :return: accepted, r
    pub fn talk(&mut self, msgs: &[u64]) -> (bool, u64) {
        if !self.active {
            panic!("Unable to prove: the protocol is not active.");
        }
        if msgs.len() != self.poly.num_multiplicands() + 1 {
            panic!(
                "Malformed message: Expect {} points, but got {}",
                self.poly.num_multiplicands() + 1,
                msgs.len()
            );
        }

        let p0 = msgs[0] % self.p;
        let p1 = msgs[1] % self.p;
        if (p0 + p1) % self.p != self.expect % self.p {
            self._reject_and_close();
            return (false, 0);
        }

        // pick r at random
        let r = self.random_r();
        let pr = interpolate(msgs, r, self.p);

        self.expect = pr;
        self.points[self.round] = r;

        // if not final step, end here
        if !(self.round + 1 == self.poly.num_variables) {
            self.round += 1;
            return (true, r);
        }

        // final step: check all
        if self.checksum_only {
            // When checksum_only is on, the verifier do not access the polynomial. It only
            // verifies that the sum of a polynomial is correct.
            // User often use this verifier as a subroutine, and uses self.subclaim() to get a sub-claim for
            // the polynomial.
            self._convince_and_close();
            return (true, r);
        }

        let final_sum = self.poly.eval(&self.points);
        if pr != final_sum {
            self._reject_and_close();
            return (false, r);
        }
        self._convince_and_close();
        (true, r)
    }

    /// The verifier should already checks the sum of the polynomial. If the sum is indeed the sum of polynomial, then
    /// the sub claim should be correct.
    /// The sub claim is in the following form:
    ///  - one point of the polynomial
    ///  - the expected evaluation at this point
    /// :return: (point: Vec<u64>, expected: u64)
    pub fn subclaim(&self) -> (Vec<u64>, u64) {
        if !self.convinced {
            panic!("The verifier is not convinced, and cannot make a sub claim.");
        }
        (self.points.clone(), self.expect)
    }

    /// Accept the sum. Close the protocol.
    fn _convince_and_close(&mut self) {
        self.convinced = true;
        self.active = false;
    }

    /// Reject the sum. Close the protocol.
    fn _reject_and_close(&mut self) {
        self.convinced = false;
        self.active = false;
    }

    // TODO: Should implement the "_repr_html_" method ???
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
