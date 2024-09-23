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
    pub fn new(poly: PMF<F>, message: Vec<Vec<F>>) -> PseudoRandomGen<F> {
        PseudoRandomGen { poly, message }
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
