use std::{fmt::Display, ops::Mul};

use halo2curves::ff::PrimeField;
use serde::{Deserialize, Serialize};

use crate::polynomial::MVLinear;

/// Product of multilinear functions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PMF<F>
where
    F: PrimeField + Clone,
{
    pub(crate) num_variables: usize,
    pub(crate) multiplicands: Vec<MVLinear<F>>,
}

impl<F> PMF<F>
where
    F: PrimeField + Clone,
{
    pub fn new(multiplicands: Vec<MVLinear<F>>) -> PMF<F> {
        if multiplicands.is_empty() {
            panic!("Multiplicands are empty.");
        }
        let mut num_variables = multiplicands[0].num_variables;
        for poly in &multiplicands {
            num_variables = num_variables.max(poly.num_variables);
        }
        PMF {
            num_variables,
            multiplicands,
        }
    }

    pub fn num_multiplicands(&self) -> usize {
        self.multiplicands.len()
    }

    pub fn eval(&self, at: &[F]) -> F {
        let mut s = F::ONE;
        for poly in &self.multiplicands {
            s *= poly.eval(at);
        }
        s
    }
}

impl<F> Mul<MVLinear<F>> for PMF<F>
where
    F: PrimeField + Clone,
{
    type Output = PMF<F>;
    fn mul(self, rhs: MVLinear<F>) -> PMF<F> {
        let mut multiplicands = self.multiplicands.clone();
        multiplicands.push(rhs);
        PMF::new(multiplicands)
    }
}

impl<F> Display for PMF<F>
where
    F: PrimeField + Clone,
{
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

// #[derive(Debug, Clone)]
// pub struct DummyPMF {
//     num_multiplicands: usize,
//     num_variables: usize,
//     pub(crate) p: u64,
// }

// impl DummyPMF {
//     pub fn new(num_multiplicands: usize, num_variables: usize, p: u64) -> DummyPMF {
//         let _pmf: PMF = PMF::new(vec![MVLinear::new(num_variables, vec![], p)]);
//         DummyPMF {
//             num_multiplicands,
//             num_variables,
//             p,
//         }
//     }

//     pub fn num_multiplicands(&self) -> usize {
//         self.num_multiplicands
//     }

//     pub fn eval(&self, _at: &[u64]) -> u64 {
//         unimplemented!("Dummy PMF Evaluated.")
//     }
// }

// impl Mul<MVLinear> for DummyPMF {
//     type Output = DummyPMF;
//     fn mul(self, _rhs: MVLinear) -> DummyPMF {
//         unimplemented!("Dummy PMF");
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;
    use halo2curves::bn256::Fr;

    #[test]
    fn test_PMF_eval() {
        let pmf: PMF<Fr> = PMF::new(vec![MVLinear::new(3, vec![(0b0, 1u64.into())])]);
        assert_eq!(pmf.eval(&[0u64.into()]), 1u64.into());
    }
}
