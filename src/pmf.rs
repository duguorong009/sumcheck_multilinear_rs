use num_bigint::BigUint;
use num_traits::One;

use crate::polynomial::MVLinear;

/// Product of multilinear functions.
struct PMF {
    num_variables: usize,
    pub(crate) p: BigUint,
    multiplicands: Vec<MVLinear>,
}

impl PMF {
    fn new(multiplicands: Vec<MVLinear>) -> PMF {
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
