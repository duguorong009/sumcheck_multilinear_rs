use std::{cell::RefCell, rc::Rc};

use halo2curves::ff::PrimeField;

use crate::{ip_pmf_verifier::{InteractivePMFVerifier, RandomGen}, pmf::PMF};

const MAX_ALLOWED_SOUNDNESS_ERROR: f64 = 2e-64;

/// Sample a random element in the field using hash function which takes the polynomial and prover message as input.
/// 
///  :param poly: The polynomial
///  :param proverMessage: List of messages [P(0), P(1), P(2), ..., P(m)]
///  :return:
fn random_element<F>(poly: &PMF<F>, prover_message: &Vec<Vec<F>>) -> F where F: PrimeField + Clone {
    // let hash_size = (F::NUM_BITS + 7) / 8;
    // let byte_length = hash_size;

    todo!()
}

#[derive(Debug, Clone)]
pub struct PseudoRandomGen<F>
where
    F: PrimeField + Clone,
{
    pub(crate) poly: PMF<F>,
    pub(crate) message: Vec<Vec<F>>,
}

impl<F> PseudoRandomGen<F> where F: PrimeField + Clone {
    pub fn new(poly: PMF<F>) -> PseudoRandomGen<F> {
        PseudoRandomGen { poly, message: vec![] }
    }
}

impl<F> RandomGen<F> for PseudoRandomGen<F> where F: PrimeField + Clone {
    fn get_random_element(&mut self) -> F {
        random_element(&self.poly, &self.message)
    }
}

/// A data structure representing offline theorem of PMF sum.
struct Theorem<F> where F: PrimeField + Clone {
    poly: PMF<F>,
    asserted_sum: F,
}

impl<F> Theorem<F> where F: PrimeField + Clone {
    pub fn new(poly: PMF<F>, asserted_sum: F) -> Theorem<F> {
        Theorem { poly, asserted_sum }
    }
}

/// A data structure representing proof of a theorem.
struct Proof<F> where F: PrimeField + Clone {
    prover_message: Vec<Vec<F>>,
}

impl<F> Proof<F> where F: PrimeField + Clone {
    pub fn new(prover_message: Vec<Vec<F>>) -> Proof<F> {
        Proof { prover_message }
    }
}

fn verify_proof<F>(theorem: &Theorem<F>, proof: &Proof<F>, max_allowed_soundness_error: Option<f64>) -> bool where F: PrimeField + Clone {
    let max_allowed_soundness_error = max_allowed_soundness_error.unwrap_or(MAX_ALLOWED_SOUNDNESS_ERROR);
    let mut gen = PseudoRandomGen::new(theorem.poly.clone());
    let mut v = InteractivePMFVerifier::new(theorem.poly.clone(), theorem.asserted_sum.clone(), Some(max_allowed_soundness_error), None, Some(Rc::new(RefCell::new(gen.clone()))));

    for msg in proof.prover_message.iter() {
        gen.message.push(msg.clone());
        v.talk(&msg);
    }
    v.convinced
}
