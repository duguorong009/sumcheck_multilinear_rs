use halo2curves::ff::PrimeField;

use crate::{
    fs_pmf_verifier::PseudoRandomGen, gkr_prover::binary_to_list,
    ip_pmf_verifier::InteractivePMFVerifier, pmf::PMF,
};

/// A linear honest prover of sum-check protocol for product of multilinear polynomials using dynamic programming.
#[derive(Debug)]
pub struct InteractivePMFProver<F>
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
        &self,
        mut As: Vec<Vec<F>>,
        verifier: &mut InteractivePMFVerifier<F>,
        gen: &mut Option<PseudoRandomGen<F>>,
    ) -> Vec<Vec<F>> {
        let l = self.poly.num_variables;
        let mut msgs: Vec<Vec<F>> = vec![];
        for i in 1..l + 1 {
            let mut products_sum: Vec<F> = vec![F::ZERO; self.poly.num_multiplicands() + 1];
            for b in 0..2usize.pow((l - i) as u32) {
                for t in 0..self.poly.num_multiplicands() + 1 {
                    let mut product = F::ONE;
                    for j in 0..self.poly.num_multiplicands() {
                        let A = As[j].clone();
                        product *= (A[b << 1] * (F::ONE - F::from_u128(t as u128)))
                            + (A[(b << 1) + 1] * F::from_u128(t as u128));
                    }
                    products_sum[t] += product;
                }
            }
            if let Some(gen) = gen {
                gen.message.push(products_sum.clone());
            } else {
                msgs.push(products_sum.clone());
            }
            let (result, r) = verifier.talk(&products_sum);

            assert!(result);

            for j in 0..self.poly.num_multiplicands() {
                for b in 0..2usize.pow((l - i) as u32) {
                    As[j][b] = As[j][b << 1] * (F::ONE - r) + As[j][(b << 1) + 1] * r;
                }
            }
        }
        if let Some(gen) = gen {
            gen.message.clone()
        } else {
            msgs
        }
    }

    /// Calculate the bookkeeping table of a single MVLinear in the PMF.
    ///
    ///  :param index: the index of the MVLinear in PMF
    ///  :return: A bookkeeping table where the index is the binary form of argument of polynomial and value is the
    ///    evaluated value
    pub fn calculate_single_table(&self, index: usize) -> Vec<F> {
        if index >= self.poly.num_multiplicands() {
            panic!(
                "PMF has onlyl {} multiplicands. index = {}",
                self.poly.num_multiplicands(),
                index
            );
        }

        let mut A: Vec<F> = vec![F::ZERO; 2usize.pow(self.poly.num_variables as u32)];
        for p in 0..2usize.pow(self.poly.num_variables as u32) {
            A[p] = self.poly.multiplicands[index].eval(&binary_to_list(
                F::from_u128(p as u128),
                self.poly.num_variables,
            ));
        }
        A
    }

    /// For all multiplicands of the PMF, calculate its bookkeeping table. The function all calculates the sum.
    ///
    ///  :return: All bookkeeping table. The sum of the PMF.
    pub fn calculate_all_bookkeeping_tables(&self) -> (Vec<Vec<F>>, F) {
        let mut S: Vec<F> = vec![F::ONE; 2usize.pow(self.poly.num_variables as u32)];
        let mut As: Vec<Vec<F>> = vec![];
        for i in 0..self.poly.num_multiplicands() {
            let A = self.calculate_single_table(i);
            for j in 0..2usize.pow(self.poly.num_variables as u32) {
                S[j] *= A[j];
            }
            As.push(A);
        }
        (As, S.iter().sum())
    }
}

#[cfg(test)]
mod tests {
    use halo2curves::bn256::Fr;

    use crate::polynomial::random_mvlinear;

    use super::*;

    #[test]
    fn test_completeness() {
        for _ in 0..100 {
            let p = PMF::<Fr>::new(vec![random_mvlinear(7); 5]);
            let pv = InteractivePMFProver::new(p.clone());
            let (As, s) = pv.calculate_all_bookkeeping_tables();
            let mut v = InteractivePMFVerifier::new(p, s, None, None, None);
            let _ = pv.attempt_prove(As, &mut v, &mut None);
        }
    }
}
