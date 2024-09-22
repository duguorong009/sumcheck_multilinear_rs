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
