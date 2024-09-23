use halo2curves::ff::PrimeField;

use crate::pmf::PMF;

#[derive(Debug)]
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

    pub fn get_random_element(&self) -> F {
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
    let gen = PseudoRandomGen::new(theorem.poly.clone());
    let v = InteractivePMFVerifier::new(theorem.poly.clone(), theorem.asserted_sum.clone(), Some(max_allowed_soundness_error), None, Some(gen));

    for msg in proof.prover_message.iter() {
        gen.message.push(msg.clone());
        v.talk(&msg);
    }
    v.convinced
}
