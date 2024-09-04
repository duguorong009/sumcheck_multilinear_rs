use std::collections::HashMap;

/// f1: Sparse polynomial f1(g, x, y) represented by a map of argument and its evaluation. Argument is little
/// endian binary form. For example, 0b10111 means f(1, 1, 1, 0, 1)
/// f2: Dense polynomial represented by a map of argument (index) and its evaluation (value).
/// f3: Dense polynomial represented by a map of argument (index) and its evaluation (value).
/// p: field size
/// l: number of variables in f2 and f3
#[derive(Debug, Clone)]
pub struct GKR {
    pub(crate) f1: HashMap<usize, u64>,
    pub(crate) f2: Vec<u64>,
    pub(crate) f3: Vec<u64>,
    pub(crate) p: u64,
    pub(crate) l: usize,
}

impl GKR {
    pub fn new(f1: HashMap<usize, u64>, f2: Vec<u64>, f3: Vec<u64>, p: u64, l: usize) -> Self {
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

        Self { f1, f2, f3, p, l }
    }
}
