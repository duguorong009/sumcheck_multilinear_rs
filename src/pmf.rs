use std::{fmt::Display, ops::Mul};

use crate::polynomial::MVLinear;

/// Product of multilinear functions.
#[derive(Debug, Clone)]
pub struct PMF {
    pub(crate) num_variables: usize,
    pub(crate) p: u64,
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

    pub fn num_multiplicands(&self) -> usize {
        self.multiplicands.len()
    }

    pub fn eval(&self, at: &[u64]) -> u64 {
        let mut s = 1;
        for poly in &self.multiplicands {
            s = (s * poly.eval(at)) % self.p.clone();
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
    pub(crate) p: u64,
}

impl DummyPMF {
    pub fn new(num_multiplicands: usize, num_variables: usize, p: u64) -> DummyPMF {
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

    pub fn eval(&self, _at: &[u64]) -> u64 {
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
        assert_eq!(pmf.eval(&[0u64.into()]), 1u64.into());
    }
}
