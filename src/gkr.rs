use std::collections::HashMap;

use halo2curves::ff::PrimeField;

/// f1: Sparse polynomial f1(g, x, y) represented by a map of argument and its evaluation. Argument is little
/// endian binary form. For example, 0b10111 means f(1, 1, 1, 0, 1)
/// f2: Dense polynomial represented by a map of argument (index) and its evaluation (value).
/// f3: Dense polynomial represented by a map of argument (index) and its evaluation (value).
/// p: field size
/// l: number of variables in f2 and f3
#[derive(Debug, Clone)]
pub struct GKR<F: PrimeField + Clone> {
    pub(crate) f1: HashMap<usize, F>,
    pub(crate) f2: Vec<F>,
    pub(crate) f3: Vec<F>,
    pub(crate) l: usize,
}

impl<F> GKR<F>
where
    F: PrimeField + Clone,
{
    pub fn new(f1: HashMap<usize, F>, f2: Vec<F>, f3: Vec<F>, l: usize) -> Self {
        assert!(f2.len() == 1 << l); // f2(x) should have size 2^L
        assert!(f3.len() == 1 << l); // f3(y) should have size 2^L

        for k in f1.keys() {
            if *k >= (1 << (3 * l)) {
                panic!(
                    "f1 has invalid term {k} cannot be represented by {} variables",
                    3 * l
                );
            }
        }

        Self { f1, f2, f3, l }
    }
}
