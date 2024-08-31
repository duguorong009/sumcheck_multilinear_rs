use rand::rngs::ThreadRng;

use crate::pmf::PMF;

#[derive(Debug)]
pub struct InteractivePMFVerifier {
    checksum_only: bool,
    p: u64,
    poly: PMF,
    asserted_sum: u64,
    rng: ThreadRng,
    active: bool,
    convinced: bool,
    points: Vec<u64>,
    round: usize,
    expect: u64,
}
