use std::{cell::RefCell, rc::Rc};

use halo2curves::ff::{FromUniformBytes, PrimeField};
use rand::rngs::OsRng;

use crate::{
    ip_pmf_verifier::{InteractivePMFVerifier, RandomGen},
    pmf::PMF,
};

const MAX_ALLOWED_SOUNDNESS_ERROR: f64 = 2e-64;

// /// Sample a random element in the field using hash function which takes the polynomial and prover message as input.
// ///
// ///  :param poly: The polynomial
// ///  :param proverMessage: List of messages [P(0), P(1), P(2), ..., P(m)]
// ///  :return:
// fn random_element<F>(poly: &PMF<F>, prover_message: &Vec<Vec<F>>) -> F where F: PrimeField + Clone {
//     let hash_size = (F::NUM_BITS as usize + 7) / 8;

//     let mut sha = blake2b_simd::Params::new().hash_length(hash_size).to_state();
//     sha.update(&serde_json::to_vec(poly).unwrap());
//     for msg in prover_message {
//         for &point in msg {
//             sha.update(b"N");
//             sha.update(&point.to_repr().as_ref());
//         }
//         sha.update(b"X");
//     }
//     let mut result = F::from_repr(sha.finalize().as_bytes().into()).unwrap();
//     // while result >= poly.p {
//     //     sha.update(b"\xFF");
//     //     result = i32::from_le_bytes(sha.finalize().try_into().unwrap());
//     // }
//     result
// }

#[derive(Debug, Clone)]
pub struct PseudoRandomGen<F>
where
    F: PrimeField + Clone,
{
    pub(crate) poly: PMF<F>,
    pub(crate) message: Vec<Vec<F>>,
}

impl<F> PseudoRandomGen<F>
where
    F: PrimeField + Clone,
{
    pub fn new(poly: PMF<F>) -> PseudoRandomGen<F> {
        PseudoRandomGen {
            poly,
            message: vec![],
        }
    }
}

impl<F> RandomGen<F> for PseudoRandomGen<F>
where
    F: PrimeField + Clone,
{
    fn get_random_element(&mut self) -> F {
        // TODO: implement `random_element` since it is Fiat-Shamir transform which should absorb the message.
        // random_element(&self.poly, &self.message)
        F::random(OsRng)
    }
}

/// A data structure representing offline theorem of PMF sum.
pub struct Theorem<F>
where
    F: PrimeField + Clone,
{
    pub(crate) poly: PMF<F>,
    pub(crate) asserted_sum: F,
}

impl<F> Theorem<F>
where
    F: PrimeField + Clone,
{
    pub fn new(poly: PMF<F>, asserted_sum: F) -> Theorem<F> {
        Theorem { poly, asserted_sum }
    }
}

/// A data structure representing proof of a theorem.
pub struct Proof<F>
where
    F: PrimeField + Clone,
{
    pub(crate) prover_message: Vec<Vec<F>>,
}

impl<F> Proof<F>
where
    F: PrimeField + Clone,
{
    pub fn new(prover_message: Vec<Vec<F>>) -> Proof<F> {
        Proof { prover_message }
    }
}

pub fn verify_proof<F>(
    theorem: &Theorem<F>,
    proof: &Proof<F>,
    max_allowed_soundness_error: Option<f64>,
) -> bool
where
    F: PrimeField + Clone,
{
    let max_allowed_soundness_error =
        max_allowed_soundness_error.unwrap_or(MAX_ALLOWED_SOUNDNESS_ERROR);
    let mut gen = PseudoRandomGen::new(theorem.poly.clone());
    let mut v = InteractivePMFVerifier::new(
        theorem.poly.clone(),
        theorem.asserted_sum.clone(),
        Some(max_allowed_soundness_error),
        None,
        Some(Rc::new(RefCell::new(gen.clone()))),
    );

    for msg in proof.prover_message.iter() {
        gen.message.push(msg.clone());
        v.talk(&msg);
    }
    v.convinced
}
