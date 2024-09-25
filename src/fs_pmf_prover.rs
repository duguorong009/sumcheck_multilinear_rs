use crate::{
    fs_pmf_verifier::{Proof, Theorem},
    ip_pmf_verifier::InteractivePMFVerifier,
    pmf::PMF,
};
use halo2curves::ff::PrimeField;

pub fn generate_theorem_and_proof<F>(
    poly: PMF<F>,
) -> (Theorem<F>, Proof<F>, InteractivePMFVerifier<F>)
where
    F: PrimeField + Clone,
{
    // unimplemented!("generate_theorem_and_proof")
}
