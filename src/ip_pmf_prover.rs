use halo2curves::ff::PrimeField;

use crate::{ip_pmf_verifier::InteractivePMFVerifier, pmf::PMF};

/// A linear honest prover of sum-check protocol for product of multilinear polynomials using dynamic programming.
#[derive(Debug)]
struct InteractivePMFProver<F>
where
    F: PrimeField + Clone,
{
    poly: PMF<F>,
}

impl<F> InteractivePMFProver<F>
where
    F: PrimeField + Clone,
{
    pub fn new(poly: PMF<F>) -> InteractivePMFProver<F> {
        InteractivePMFProver { poly }
    }

    /// Attempt to prove the sum.
    ///  :param As: The bookkeeping table for each MVLinear in the PMF
    ///  :param verifier: the active interactive PMF verifier instance
    ///  :param gen: in FS mode, attemptProve will provide source of randomness for pseudorandom generator
    ///  :return: the prover message
    pub fn attempt_prove(
        &mut self,
        As: Vec<Vec<F>>,
        verifier: InteractivePMFVerifier<F>,
        // gen: Option<PseudoRandomGen>,
    ) -> Vec<Vec<F>> {
        let l = self.poly.num_variables;
        let mut msgs: Vec<Vec<F>> = vec![];
        for i in 1..l + 1 {
            let mut products_sum: Vec<F> = vec![F::ZERO; self.poly.num_multiplicands() + 1];
            for b in 0..2usize.pow((l - i) as u32) {
                for t in 0..self.poly.num_multiplicands() + 1 {
                    let mut product = F::ONE;
                    for j in 0..self.poly.num_multiplicands() {
                        let A = As[j];
                        product *= (A[b << 1] * (F::ONE - F::from_u128(t as u128))) + (A[(b << 1) + 1] * F::from_u128(t as u128));
                    }
                    products_sum[t] += product;
                }
            }
            if let Some(gen) = gen {
                gen.message.push(products_sum.clone());
            } else {
                msgs.push(products_sum);
            }
            let (result, r) = verifier.talk(&products_sum);

            assert!(result);

            for j in 0..self.poly.num_multiplicands() {
                for b in 0..2usize.pow((l - i) as u32) {
                    As[j][b] = As[j][b << 1] * (F::ONE - r) + As[(j << 1) + 1][b << 1] * r;
                }
            }
        }
        if let Some(gen) = gen {
            gen.message
        } else {
            msgs
        }
    }

    pub fn calculate_single_table(&mut self) {
        todo!()
    }

    pub fn calculate_all_bookkeeping_tables(&mut self) {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_completeness() {
        todo!()
    }
}
