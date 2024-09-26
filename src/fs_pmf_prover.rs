use std::{cell::RefCell, rc::Rc};

use crate::{
    fs_pmf_verifier::{Proof, PseudoRandomGen, Theorem},
    ip_pmf_prover::InteractivePMFProver,
    ip_pmf_verifier::InteractivePMFVerifier,
    pmf::PMF,
};
use halo2curves::ff::PrimeField;

/// Generate the theorem(poly itself and the asserted sum) and its proof.
///
///  :param maxAllowedSoundnessError:
///  :param poly: The PMF polynomial
///  :return: theorem, proof, and the (hopefully) convinced pseudorandom verifier
pub fn generate_theorem_and_proof<F>(
    poly: PMF<F>,
) -> (Theorem<F>, Proof<F>, InteractivePMFVerifier<F>)
where
    F: PrimeField + Clone,
{
    let pv = InteractivePMFProver::new(poly.clone());
    let (As, s) = pv.calculate_all_bookkeeping_tables();

    let gen = PseudoRandomGen::new(poly.clone());
    let mut v = InteractivePMFVerifier::new(
        poly.clone(),
        s,
        None,
        None,
        Some(Rc::new(RefCell::new(gen.clone()))),
    );
    let msgs = pv.attempt_prove(As, &mut v, &mut Some(gen));

    let theorem = Theorem::new(poly, s);
    let proof = Proof::new(msgs);
    (theorem, proof, v)
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::{fs_pmf_verifier::verify_proof, polynomial::random_mvlinear};

    use halo2curves::bn256::Fr;
    
    #[test]
    fn test_completeness() {
        for _ in 0..100 {
            let p = PMF::<Fr>::new(vec![random_mvlinear(7); 5]);
            let (theorem, proof, _) = generate_theorem_and_proof(p);
            assert!(verify_proof(&theorem, &proof, None));
        }
    }
}