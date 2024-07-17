use std::{
    collections::HashMap,
    ops::{Add, Mul, Neg, Sub},
};

use num_bigint::BigUint;
use num_traits::{One, Zero};

#[derive(Debug, Clone)]
/// A Sparse Representation of a multi-linear polynomial.
struct MVLinear {
    num_variables: usize,
    terms: HashMap<usize, BigUint>,
    p: BigUint,
}

impl MVLinear {
    fn new(num_variables: usize, term: Vec<(usize, BigUint)>, p: BigUint) -> MVLinear {
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

    fn eval(&self, at: Vec<BigUint>) -> BigUint {
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
    fn eval_bin(&self, at: usize) -> BigUint {
        if at > 2usize.pow(self.num_variables.try_into().unwrap()) {
            panic!("Number of varialbes is larger than expected")
        }
        let mut args = vec![BigUint::zero(); self.num_variables];
        for i in 0..self.num_variables {
            if at & (1 << i) > 0 {
                args[i] = BigUint::one();
            }
        }
        self.eval(args)
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
    assert_eq!(p1.eval(at), expected);
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
