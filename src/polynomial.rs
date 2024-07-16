use std::{collections::HashMap, ops::Add, str::FromStr};

use num_bigint::BigUint;
use num_traits::Zero;

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
        assert_eq!(self.p, other.p, "The function being added is not in the same field.");
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
