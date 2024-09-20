use halo2curves::ff::PrimeField;

use crate::pmf::PMF;


#[derive(Debug)]
struct InteractivePMFProver<F> where F: PrimeField + Clone {
    poly: PMF<F>
}

impl<F> InteractivePMFProver<F>
where
    F: PrimeField + Clone
{
    pub fn new(poly: PMF<F>) -> InteractivePMFProver<F> {
        InteractivePMFProver { poly }
    }

    pub fn attempt_prove(&mut self) {
        todo!()
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
