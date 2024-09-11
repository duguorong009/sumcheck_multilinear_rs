use std::{
    collections::HashMap,
    ops::{Add, Mul, MulAssign, Neg, Sub},
};

use halo2curves::ff::PrimeField;
use rand::{
    distributions::{Distribution, Uniform},
    Rng,
};

use crate::pmf::PMF;

#[derive(Debug, Clone)]
/// A Sparse Representation of a multi-linear polynomial.
pub struct MVLinear<F>
where
    F: PrimeField + Clone,
{
    pub(crate) num_variables: usize,
    pub(crate) terms: HashMap<usize, F>,
}

impl<F> MVLinear<F>
where
    F: PrimeField + Clone,
{
    pub fn new(num_variables: usize, term: Vec<(usize, F)>) -> MVLinear<F> {
        let mut terms = HashMap::new();
        for (k, v) in term {
            if k >> num_variables > 0 {
                panic!("Term is out of range.");
            }
            if v == F::ZERO {
                continue;
            }
            if terms.contains_key(&k) {
                terms.insert(k, *terms.get(&k).unwrap() + v);
            } else {
                terms.insert(k, v);
            }
        }

        MVLinear {
            num_variables,
            terms,
        }
    }

    // fn assert_same_type(&self, other: &MVLinear<F>) {
    //     assert_eq!(
    //         self.p, other.p,
    //         "The function being added is not in the same field."
    //     );
    // }

    pub fn eval(&self, at: &[F]) -> F {
        let mut s = F::ZERO;
        for &term in self.terms.keys() {
            let mut term = term;
            let mut i = 0;
            let mut val = *self.terms.get(&term).unwrap();
            while term != 0 {
                if term & 1 == 1 {
                    val *= at[i];
                }
                if val == F::ZERO {
                    break;
                }
                term >>= 1;
                i += 1;
            }
            s += val;
        }
        s
    }

    /// Evaluate the polynomial where the arguments are in {0, 1}. The ith argument is the ith bit of the polynomial.
    ///
    /// `at`: polynomial argument in binary form
    pub fn eval_bin(&self, at: usize) -> F {
        if at > 2usize.pow(self.num_variables.try_into().unwrap()) {
            panic!("Number of varialbes is larger than expected")
        }
        let mut args = vec![F::ZERO; self.num_variables];
        for i in 0..self.num_variables {
            if at & (1 << i) > 0 {
                args[i] = F::ONE;
            }
        }
        self.eval(&args)
    }

    /// Evaluate part of the arguments of the multilinear polynomial.
    ///
    /// `args`: arguments at beginning
    pub fn eval_part(&self, args: &[F]) -> MVLinear<F> {
        let s = args.len();
        if s > self.num_variables {
            panic!("len(args) > self.num_variables");
        }
        let mut new_terms: HashMap<usize, F> = HashMap::new();
        for (t, v) in self.terms.iter() {
            let mut t = *t;
            let mut v = *v;
            for k in 0..s {
                if t & (1 << k) > 0 {
                    v *= args[k];
                    t &= !(1 << k);
                }
            }
            let t_shifted = t >> s;
            new_terms.insert(
                t_shifted,
                *new_terms.get(&t_shifted).unwrap_or(&F::ZERO) + v,
            );
        }
        MVLinear::new(
            self.num_variables - args.len(),
            new_terms.into_iter().collect(),
        )
    }

    /// Remove redundant unused variable from left.
    ///
    /// `n`: number of variables to collapse
    fn collapse_left(&self, n: usize) -> MVLinear<F> {
        let mut new_terms = HashMap::new();
        let mask = (1 << n) - 1;
        for (t, v) in self.terms.iter() {
            if t & mask > 0 {
                panic!("Cannot collapse: Variable exist.");
            }
            new_terms.insert(t >> n, *v);
        }
        MVLinear::new(self.num_variables - n, new_terms.into_iter().collect())
    }

    /// Remove redundant unused variable from right.
    ///
    /// `n`: number of variables to collapse
    fn collapse_right(&self, n: usize) -> MVLinear<F> {
        let mut new_terms = HashMap::new();
        let mask = ((1 << n) - 1) << (self.num_variables - n);
        let anti_mask = (1 << (self.num_variables - n)) - 1;
        for (t, v) in self.terms.iter() {
            if t & mask > 0 {
                panic!("Cannot collapse: Variable exist.");
            }
            new_terms.insert(t & anti_mask, *v);
        }
        MVLinear::new(self.num_variables - n, new_terms.into_iter().collect())
    }
}

impl<F> Add for MVLinear<F>
where
    F: PrimeField + Clone,
{
    type Output = MVLinear<F>;
    fn add(self, other: MVLinear<F>) -> MVLinear<F> {
        // self.assert_same_type(&other);
        let mut ans = self.clone();
        ans.num_variables = self.num_variables.max(other.num_variables);
        for (k, v) in other.terms {
            if self.terms.contains_key(&k) {
                ans.terms.insert(k, *self.terms.get(&k).unwrap() + v);
                if ans.terms.get(&k).unwrap() == &F::ZERO {
                    ans.terms.remove(&k);
                }
            } else {
                ans.terms.insert(k, v);
            }
        }
        ans
    }
}

impl<F> Sub for MVLinear<F>
where
    F: PrimeField + Clone,
{
    type Output = MVLinear<F>;
    fn sub(self, other: MVLinear<F>) -> MVLinear<F> {
        // self.assert_same_type(&other);
        let mut ans = self.clone();
        ans.num_variables = self.num_variables.max(other.num_variables);
        for (k, v) in other.terms {
            if self.terms.contains_key(&k) {
                ans.terms.insert(k, *self.terms.get(&k).unwrap() - v);
                if ans.terms.get(&k).unwrap() == &F::ZERO {
                    ans.terms.remove(&k);
                }
            } else {
                ans.terms.insert(k, -v);
            }
        }
        ans
    }
}

impl<F> Neg for MVLinear<F>
where
    F: PrimeField + Clone,
{
    type Output = MVLinear<F>;
    fn neg(self) -> MVLinear<F> {
        let zero = MVLinear::new(self.num_variables, vec![(0b0, F::ZERO)]);
        zero - self
    }
}

impl<F> Mul for MVLinear<F>
where
    F: PrimeField + Clone,
{
    type Output = MVLinear<F>;
    fn mul(self, other: MVLinear<F>) -> MVLinear<F> {
        // self.assert_same_type(&other);
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
                        *terms.get(&nk).unwrap()
                            + *self.terms.get(sk).unwrap() * *other.terms.get(ok).unwrap(),
                    );
                } else {
                    terms.insert(
                        nk,
                        *self.terms.get(sk).unwrap() * *other.terms.get(ok).unwrap(),
                    );
                }
                if terms.get(&nk).unwrap() == &F::ZERO {
                    terms.remove(&nk);
                }
            }
        }
        let terms = terms.into_iter().collect();
        MVLinear::new(self.num_variables.max(other.num_variables), terms)
    }
}

impl<F> From<F> for MVLinear<F>
where
    F: PrimeField + Clone,
{
    fn from(x: F) -> MVLinear<F> {
        MVLinear::new(0, vec![(0b0, x)])
    }
}

impl<F> MulAssign for MVLinear<F>
where
    F: PrimeField + Clone,
{
    fn mul_assign(&mut self, other: MVLinear<F>) {
        *self = self.clone() * other;
    }
}

impl<F> Mul<PMF<F>> for MVLinear<F> where F: PrimeField + Clone {
    type Output = PMF<F>;
    fn mul(self, rhs: PMF<F>) -> PMF<F> {
        let mut multiplicands = rhs.multiplicands.clone();
        multiplicands.push(self);
        PMF::new(multiplicands)
    }
}

impl<F> MulAssign<F> for MVLinear<F>
where
    F: PrimeField + Clone,
{
    fn mul_assign(&mut self, rhs: F) {
        let rhs: MVLinear<F> = rhs.into();
        *self = self.clone() * rhs;
    }
}

impl<F> PartialEq for MVLinear<F>
where
    F: PrimeField + Clone,
{
    fn eq(&self, other: &MVLinear<F>) -> bool {
        let diff = self.clone() - other.clone();
        diff.terms.is_empty()
    }
}

impl<F> std::fmt::Display for MVLinear<F>
where
    F: PrimeField + Clone,
{
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
            s += format!("{:?}", self.terms.get(k)).as_str();
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
pub fn make_mvlinear_constructor<F: PrimeField + Clone>(
    num_variables: usize,
) -> impl Fn(Vec<(usize, F)>) -> MVLinear<F> {
    move |terms: Vec<(usize, F)>| MVLinear::new(num_variables, terms)
}

// Function to generate a random prime number of a given bit length
// NOTE: Currently only works for bit_length <= 64
pub fn random_prime(bit_length: usize) -> u64 {
    let mut rng = rand::thread_rng();
    let lower_bound = 1u64 << (bit_length - 1); // 2 ^ (bit_length - 1)
    let upper_bound = (1u64 << bit_length) - 1; // 2 ^ bit_length - 1
    let mut prime_candidate: u64 = rng.gen_range(lower_bound..=upper_bound) | 1;
    while !is_prime::is_prime(&prime_candidate.to_string()) {
        prime_candidate += 2u64;
    }
    prime_candidate
}

// Function to create a random MVLinear
pub fn random_mvlinear<F>(
    num_variables: usize,
    prime: Option<u64>,
    prime_bit_length: Option<usize>,
) -> MVLinear<F>
where
    F: PrimeField + Clone,
{
    let prime_bit_length = prime_bit_length.unwrap_or(64);
    let num_terms = 2_usize.pow(num_variables as u32);
    let p = match prime {
        Some(p) => p,
        None => random_prime(prime_bit_length), // Ensure the prime fits in usize
    };

    let mv_linear_constructor = make_mvlinear_constructor(num_variables);
    let mut terms = HashMap::new();
    let mut rng = rand::rngs::OsRng;
    let range = Uniform::from(0..num_terms);

    for _ in 0..num_terms {
        terms.insert(range.sample(&mut rng), F::random(rng));
    }

    mv_linear_constructor(terms.into_iter().collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use halo2curves::bn256::Fr;

    #[test]
    fn test_mvlinear_new() {
        // P(x0, x1, x2, x3) = 15 + x0 + 4 * x3 + x1 * x2 + 5 * x2 * x3 (in Z_37)
        let p: MVLinear<Fr> = MVLinear::new(
            4,
            vec![
                (0b0000, 15u64.into()),
                (0b0001, 1u64.into()),
                (0b1000, 4u64.into()),
                (0b0110, 1u64.into()),
                (0b1100, 5u64.into()),
            ],
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
        let p1: MVLinear<Fr> = MVLinear::new(4, vec![(0b0000, 15u64.into())]);
        let p2: MVLinear<Fr> = MVLinear::new(4, vec![(0b0001, 1u64.into())]);
        let p3: MVLinear<Fr> = p1 + p2;
        let expected = MVLinear::new(4, vec![(0b0000, 15u64.into()), (0b0001, 1u64.into())]);
        assert_eq!(p3, expected);
    }

    #[test]
    fn test_mvlinear_sub() {
        let p1: MVLinear<Fr> = MVLinear::new(4, vec![(0b0000, 15u64.into())]);
        let p2: MVLinear<Fr> = MVLinear::new(4, vec![(0b0001, 1u64.into())]);
        let p3: MVLinear<Fr> = p1 - p2;
        let expected = MVLinear::new(4, vec![(0b0000, 15u64.into()), (0b0001, 36u64.into())]);
        assert_eq!(p3, expected);
    }

    #[test]
    fn test_mvlinear_neg() {
        let p1: MVLinear<Fr> = MVLinear::new(4, vec![(0b0000, 15u64.into())]);
        let p2: MVLinear<Fr> = -p1;
        let expected = MVLinear::new(4, vec![(0b0000, 22u64.into())]);
        assert_eq!(p2, expected);
    }

    #[test]
    fn test_mvlinear_mul() {
        let p1: MVLinear<Fr> = MVLinear::new(4, vec![(0b0000, 15u64.into())]);
        let p2: MVLinear<Fr> = MVLinear::new(4, vec![(0b0001, 1u64.into())]);
        let p3: MVLinear<Fr> = p1 * p2;
        let expected = MVLinear::new(4, vec![(0b0001, 15u64.into())]);
        assert_eq!(p3, expected);
    }

    #[test]
    fn test_mvlinear_eval() {
        let p1: MVLinear<Fr> = MVLinear::new(4, vec![(0b0000, 15u64.into()), (0b0001, 1u64.into())]);
        let at = vec![1u64.into(), 1u64.into(), 1u64.into(), 1u64.into()];
        let expected = (15u64 + 1u64).into();
        assert_eq!(p1.eval(&at), expected);
    }

    #[test]
    fn test_mvlinear_eval_bin() {
        let p1: MVLinear<Fr> = MVLinear::new(4, vec![(0b0000, 15u64.into()), (0b0001, 1u64.into())]);
        let expected = (15u64 + 1u64).into();
        assert_eq!(p1.eval_bin(0b0001), expected);
        let expected = 15u64.into();
        assert_eq!(p1.eval_bin(0b1110), expected);
    }

    #[test]
    fn test_mvlinear_collpase() {
        // right collpase: remove x3, x2, x1
        let p1: MVLinear<Fr> = MVLinear::new(4, vec![(0b0000, 15u64.into()), (0b0001, 1u64.into())]);
        let expected = MVLinear::new(1, vec![(0b0000, 15u64.into()), (0b0001, 1u64.into())]);
        assert_eq!(p1.collapse_right(3), expected);

        // left collpase: remove x0, x1, x2
        let p1: MVLinear<Fr> = MVLinear::new(4, vec![(0b1000, 15u64.into())]);
        let expected = MVLinear::new(1, vec![(0b0001, 15u64.into())]);
        assert_eq!(p1.collapse_left(3), expected);
    }
}
