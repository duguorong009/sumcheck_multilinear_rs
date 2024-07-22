use std::{clone, fmt::Display, ops::Mul};

use num_bigint::BigUint;
use num_traits::One;

use crate::polynomial::MVLinear;

/// Product of multilinear functions.
#[derive(Debug, Clone)]
pub struct PMF {
    num_variables: usize,
    pub(crate) p: BigUint,
    pub(crate) multiplicands: Vec<MVLinear>,
}

impl PMF {
    pub fn new(multiplicands: Vec<MVLinear>) -> PMF {
        if multiplicands.is_empty() {
            panic!("Multiplicands are empty.");
        }
        let mut num_variables = multiplicands[0].num_variables;
        let p = multiplicands[0].p.clone();
        for poly in &multiplicands {
            num_variables = num_variables.max(poly.num_variables);
            if poly.p != p {
                panic!("Field size mismatch.");
            }
        }
        PMF {
            num_variables,
            p,
            multiplicands,
        }
    }

    fn num_multiplicands(&self) -> usize {
        self.multiplicands.len()
    }

    fn eval(&self, at: Vec<BigUint>) -> BigUint {
        let mut s = BigUint::one();
        for poly in &self.multiplicands {
            s = (s * poly.eval(at.clone())) % self.p.clone();
        }
        s
    }
}

impl Mul<MVLinear> for PMF {
    type Output = PMF;
    fn mul(self, rhs: MVLinear) -> PMF {
        let mut multiplicands = self.multiplicands.clone();
        multiplicands.push(rhs);
        PMF::new(multiplicands)
    }
}

impl Display for PMF {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut s = String::new();
        s += "Product[";
        for poly in &self.multiplicands {
            s += format!("{}", poly).as_str();
            s += ",";
        }
        s += "]";
        write!(f, "{}", s)
    }
}

#[derive(Debug, Clone)]
pub struct DummyPMF {
    num_multiplicands: usize,
    num_variables: usize,
    pub(crate) p: BigUint,
}

impl DummyPMF {
    pub fn new(num_multiplicands: usize, num_variables: usize, p: BigUint) -> DummyPMF {
        // let _pmf: PMF = PMF::new(vec![MVLinear::new(num_variables, vec![], p.clone())]);
        DummyPMF {
            num_multiplicands,
            num_variables,
            p,
        }
    }

    pub fn num_multiplicands(&self) -> usize {
        self.num_multiplicands
    }

    pub fn eval(&self, at: Vec<BigUint>) -> BigUint {
        unimplemented!("Dummy PMF Evaluated.")
    }
}

impl Mul<MVLinear> for DummyPMF {
    type Output = DummyPMF;
    fn mul(self, _rhs: MVLinear) -> DummyPMF {
        unimplemented!("Dummy PMF");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_PMF_eval() {
        let mut pmf = PMF::new(vec![MVLinear::new(
            3,
            vec![(0b0, 1u64.into())],
            2u64.into(),
        )]);
        pmf.p = 7u64.into();
        assert_eq!(pmf.eval(vec![0u64.into()]), 1u64.into());
    }
}
