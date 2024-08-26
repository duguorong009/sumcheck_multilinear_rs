use std::{
    collections::HashMap,
    ops::{Add, Mul, MulAssign, Neg, Sub},
};

use num_bigint::{BigUint, RandBigInt};
use num_traits::{One, Zero};
use rand::distributions::{Distribution, Uniform};

use crate::pmf::PMF;

#[derive(Debug, Clone)]
/// A Sparse Representation of a multi-linear polynomial.
pub struct MVLinear {
    pub(crate) num_variables: usize,
    pub(crate) terms: HashMap<usize, BigUint>,
    pub(crate) p: BigUint,
}

impl MVLinear {
    pub fn new(num_variables: usize, term: Vec<(usize, BigUint)>, p: BigUint) -> MVLinear {
        let mut terms = HashMap::new();
        for (k, v) in term {
            if k >> num_variables > 0 {
                panic!("Term is out of range.");
            }
            if (v.clone() % p.clone()).is_zero() {
                continue;
            }
            if terms.contains_key(&k) {
                terms.insert(k, (terms.get(&k).unwrap() + v) % p.clone());
            } else {
                terms.insert(k, v % p.clone());
            }
        }

        MVLinear {
            num_variables,
            terms,
            p,
        }
    }

    fn assert_same_type(&self, other: &MVLinear) {
        assert_eq!(
            self.p, other.p,
            "The function being added is not in the same field."
        );
    }

    pub fn eval(&self, at: &[BigUint]) -> BigUint {
        let mut s = BigUint::zero();
        for &term in self.terms.keys() {
            let mut term = term;
            let mut i = 0;
            let mut val = self.terms.get(&term).unwrap().clone();
            while term != 0 {
                if term & 1 == 1 {
                    val = (val * (at[i].clone() % self.p.clone())) % self.p.clone();
                }
                if val.is_zero() {
                    break;
                }
                term >>= 1;
                i += 1;
            }
            s = (s + val) % self.p.clone();
        }
        s
    }

    /// Evaluate the polynomial where the arguments are in {0, 1}. The ith argument is the ith bit of the polynomial.
    ///
    /// `at`: polynomial argument in binary form
    pub fn eval_bin(&self, at: usize) -> BigUint {
        if at > 2usize.pow(self.num_variables.try_into().unwrap()) {
            panic!("Number of varialbes is larger than expected")
        }
        let mut args = vec![BigUint::zero(); self.num_variables];
        for i in 0..self.num_variables {
            if at & (1 << i) > 0 {
                args[i] = BigUint::one();
            }
        }
        self.eval(&args)
    }

    /// Evaluate part of the arguments of the multilinear polynomial.
    ///
    /// `args`: arguments at beginning
    fn eval_part(&self, args: &[BigUint]) -> MVLinear {
        let s = args.len();
        if s > self.num_variables {
            panic!("len(args) > self.num_variables");
        }
        let mut new_terms = HashMap::new();
        for (t, v) in self.terms.iter() {
            let mut t = *t;
            let mut v = v.clone();
            for k in 0..s {
                if t & (1 << k) > 0 {
                    v = v * (args[k].clone() % self.p.clone()) % self.p.clone();
                    t &= !(1 << k);
                }
            }
            let t_shifted = t >> s;
            new_terms.insert(
                t_shifted,
                (new_terms.get(&t_shifted).unwrap_or(&0u64.into()) + v) % self.p.clone(),
            );
        }
        MVLinear::new(
            self.num_variables - args.len(),
            new_terms.into_iter().collect(),
            self.p.clone(),
        )
    }

    /// Remove redundant unused variable from left.
    ///
    /// `n`: number of variables to collapse
    fn collapse_left(&self, n: usize) -> MVLinear {
        let mut new_terms = HashMap::new();
        let mask = (1 << n) - 1;
        for (t, v) in self.terms.iter() {
            if t & mask > 0 {
                panic!("Cannot collapse: Variable exist.");
            }
            new_terms.insert(t >> n, v.clone());
        }
        MVLinear::new(
            self.num_variables - n,
            new_terms.into_iter().collect(),
            self.p.clone(),
        )
    }

    /// Remove redundant unused variable from right.
    ///
    /// `n`: number of variables to collapse
    fn collapse_right(&self, n: usize) -> MVLinear {
        let mut new_terms = HashMap::new();
        let mask = ((1 << n) - 1) << (self.num_variables - n);
        let anti_mask = (1 << (self.num_variables - n)) - 1;
        for (t, v) in self.terms.iter() {
            if t & mask > 0 {
                panic!("Cannot collapse: Variable exist.");
            }
            new_terms.insert(t & anti_mask, v.clone());
        }
        MVLinear::new(
            self.num_variables - n,
            new_terms.into_iter().collect(),
            self.p.clone(),
        )
    }
}

impl Add for MVLinear {
    type Output = MVLinear;
    fn add(self, other: MVLinear) -> MVLinear {
        self.assert_same_type(&other);
        let mut ans = self.clone();
        ans.num_variables = self.num_variables.max(other.num_variables);
        for (k, v) in other.terms {
            if self.terms.contains_key(&k) {
                ans.terms
                    .insert(k, (self.terms.get(&k).unwrap() + v) % self.p.clone());
                if ans.terms.get(&k).unwrap().is_zero() {
                    ans.terms.remove(&k);
                }
            } else {
                ans.terms.insert(k, v % self.p.clone());
            }
        }
        ans
    }
}

impl Sub for MVLinear {
    type Output = MVLinear;
    fn sub(self, other: MVLinear) -> MVLinear {
        self.assert_same_type(&other);
        let mut ans = self.clone();
        ans.num_variables = self.num_variables.max(other.num_variables);
        for (k, v) in other.terms {
            if self.terms.contains_key(&k) {
                ans.terms
                    .insert(k, (self.terms.get(&k).unwrap() - v) % self.p.clone());
                if ans.terms.get(&k).unwrap().is_zero() {
                    ans.terms.remove(&k);
                }
            } else {
                ans.terms.insert(k, self.p.clone() - v);
            }
        }
        ans
    }
}

impl Neg for MVLinear {
    type Output = MVLinear;
    fn neg(self) -> MVLinear {
        let zero = MVLinear::new(
            self.num_variables,
            vec![(0b0, BigUint::zero())],
            self.p.clone(),
        );
        zero - self
    }
}

impl Mul for MVLinear {
    type Output = MVLinear;
    fn mul(self, other: MVLinear) -> MVLinear {
        self.assert_same_type(&other);
        let mut terms = HashMap::new();
        // native n^2 poly multiplication where n is number of terms
        for sk in self.terms.keys() {
            for ok in other.terms.keys() {
                if sk & ok > 0 {
                    panic!("The product is no longer multi-linear function.")
                }
                let nk = sk + ok;
                if terms.contains_key(&nk) {
                    terms.insert(
                        nk,
                        (terms.get(&nk).unwrap()
                            + self.terms.get(sk).unwrap() * other.terms.get(ok).unwrap())
                            % self.p.clone(),
                    );
                } else {
                    terms.insert(
                        nk,
                        self.terms.get(sk).unwrap() * other.terms.get(ok).unwrap() % self.p.clone(),
                    );
                }
                if terms.get(&nk).unwrap().is_zero() {
                    terms.remove(&nk);
                }
            }
        }
        let terms = terms.into_iter().collect();
        MVLinear::new(self.num_variables.max(other.num_variables), terms, self.p)
    }
}

impl Sub<MVLinear> for u64 {
    type Output = MVLinear;
    fn sub(self, other: MVLinear) -> MVLinear {
        let t = MVLinear::new(
            other.num_variables,
            vec![(0b0, BigUint::from(self))],
            other.p.clone(),
        );
        t - other
    }
}

impl MulAssign for MVLinear {
    fn mul_assign(&mut self, other: MVLinear) {
        *self = self.clone() * other;
    }
}

impl Mul<PMF> for MVLinear {
    type Output = PMF;
    fn mul(self, rhs: PMF) -> PMF {
        let mut multiplicands = rhs.multiplicands.clone();
        multiplicands.push(self);
        PMF::new(multiplicands)
    }
}

impl PartialEq for MVLinear {
    fn eq(&self, other: &MVLinear) -> bool {
        let diff = self.clone() - other.clone();
        diff.terms.is_empty()
    }
}

impl std::fmt::Display for MVLinear {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut s = String::new();
        s += "MVLinear(";
        let mut limit = 8;
        for k in self.terms.keys() {
            if limit != 8 {
                s += " + ";
            }
            if limit == 0 {
                s += "...";
                break;
            }
            s += self.terms[k].to_string().as_str();
            if *k != 0 {
                s += "*"
            }

            let mut i = 0;
            let mut k = *k;
            while k != 0 {
                if k & 1 == 1 {
                    s += format!("x{}", i).as_str();
                }
                i += 1;
                k >>= 1;
            }
            limit -= 1;
        }
        s += ")";
        write!(f, "{}", s)
    }
}

/// Return a function that outputs MVLinear.
pub fn make_mvlinear_constructor(
    num_variables: usize,
    p: BigUint,
) -> impl Fn(Vec<(usize, BigUint)>) -> MVLinear {
    move |terms: Vec<(usize, BigUint)>| MVLinear::new(num_variables, terms, p.clone())
}

// Function to generate a random prime number of a given bit length
pub fn random_prime(bit_length: usize) -> BigUint {
    let mut rng = rand::thread_rng();
    let mut prime_candidate = rng.gen_biguint(bit_length.try_into().unwrap());
    while !is_prime::is_prime(&prime_candidate.to_str_radix(10)) {
        prime_candidate += 2u64;
    }
    prime_candidate
}

// Function to create a random MVLinear
fn random_mvlinear(
    num_variables: usize,
    prime: Option<BigUint>,
    prime_bit_length: Option<usize>,
) -> MVLinear {
    let prime_bit_length = prime_bit_length.unwrap_or(128);
    let num_terms = 2_usize.pow(num_variables as u32);
    let p = match prime {
        Some(p) => p,
        None => random_prime(prime_bit_length), // Ensure the prime fits in usize
    };

    let mv_linear_constructor = make_mvlinear_constructor(num_variables, p.clone());
    let mut terms = HashMap::new();
    let mut rng = rand::thread_rng();
    let range = Uniform::from(0..num_terms);

    for _ in 0..num_terms {
        terms.insert(
            range.sample(&mut rng),
            rng.gen_biguint_range(&0u64.into(), &p),
        );
    }

    mv_linear_constructor(terms.into_iter().collect())
}

#[test]
fn test_mvlinear_new() {
    // P(x0, x1, x2, x3) = 15 + x0 + 4 * x3 + x1 * x2 + 5 * x2 * x3 (in Z_37)
    let p = MVLinear::new(
        4,
        vec![
            (0b0000, 15u64.into()),
            (0b0001, 1u64.into()),
            (0b1000, 4u64.into()),
            (0b0110, 1u64.into()),
            (0b1100, 5u64.into()),
        ],
        37u64.into(),
    );
    assert_eq!(
        p.terms,
        HashMap::from([
            (0, 15u64.into()),
            (1, 1u64.into()),
            (6, 1u64.into()),
            (8, 4u64.into()),
            (12, 5u64.into()),
        ])
    );
    // println!("p: {}", p);
}

#[test]
fn test_mvlinear_add() {
    let p1 = MVLinear::new(4, vec![(0b0000, 15u64.into())], 37u64.into());
    let p2 = MVLinear::new(4, vec![(0b0001, 1u64.into())], 37u64.into());
    let p3 = p1 + p2;
    let expected = MVLinear::new(
        4,
        vec![(0b0000, 15u64.into()), (0b0001, 1u64.into())],
        37u64.into(),
    );
    assert_eq!(p3, expected);
}

#[test]
fn test_mvlinear_sub() {
    let p1 = MVLinear::new(4, vec![(0b0000, 15u64.into())], 37u64.into());
    let p2 = MVLinear::new(4, vec![(0b0001, 1u64.into())], 37u64.into());
    let p3 = p1 - p2;
    let expected = MVLinear::new(
        4,
        vec![(0b0000, 15u64.into()), (0b0001, 36u64.into())],
        37u64.into(),
    );
    assert_eq!(p3, expected);
}

#[test]
fn test_mvlinear_neg() {
    let p1 = MVLinear::new(4, vec![(0b0000, 15u64.into())], 37u64.into());
    let p2 = -p1;
    let expected = MVLinear::new(4, vec![(0b0000, 22u64.into())], 37u64.into());
    assert_eq!(p2, expected);
}

#[test]
fn test_mvlinear_mul() {
    let p1 = MVLinear::new(4, vec![(0b0000, 15u64.into())], 37u64.into());
    let p2 = MVLinear::new(4, vec![(0b0001, 1u64.into())], 37u64.into());
    let p3 = p1 * p2;
    let expected = MVLinear::new(4, vec![(0b0001, 15u64.into())], 37u64.into());
    assert_eq!(p3, expected);
}

#[test]
fn test_mvlinear_eval() {
    let p1 = MVLinear::new(
        4,
        vec![(0b0000, 15u64.into()), (0b0001, 1u64.into())],
        37u64.into(),
    );
    let at = vec![1u64.into(), 1u64.into(), 1u64.into(), 1u64.into()];
    let expected = (15u64 + 1u64).into();
    assert_eq!(p1.eval(&at), expected);
}

#[test]
fn test_mvlinear_eval_bin() {
    let p1 = MVLinear::new(
        4,
        vec![(0b0000, 15u64.into()), (0b0001, 1u64.into())],
        37u64.into(),
    );
    let expected = (15u64 + 1u64).into();
    assert_eq!(p1.eval_bin(0b0001), expected);
    let expected = 15u64.into();
    assert_eq!(p1.eval_bin(0b1110), expected);
}

#[test]
fn test_mvlinear_collpase() {
    // right collpase: remove x3, x2, x1
    let p1 = MVLinear::new(
        4,
        vec![(0b0000, 15u64.into()), (0b0001, 1u64.into())],
        37u64.into(),
    );
    let expected = MVLinear::new(
        1,
        vec![(0b0000, 15u64.into()), (0b0001, 1u64.into())],
        37u64.into(),
    );
    assert_eq!(p1.collapse_right(3), expected);

    // left collpase: remove x0, x1, x2
    let p1 = MVLinear::new(4, vec![(0b1000, 15u64.into())], 37u64.into());
    let expected = MVLinear::new(1, vec![(0b0001, 15u64.into())], 37u64.into());
    assert_eq!(p1.collapse_left(3), expected);
}
