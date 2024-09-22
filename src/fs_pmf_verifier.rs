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
