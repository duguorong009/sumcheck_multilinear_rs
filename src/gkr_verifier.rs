use rand::rngs::ThreadRng;
use std::collections::HashMap;

use crate::ip_pmf_verifier::InteractivePMFVerifier;

#[derive(Debug)]
pub enum GKRVerifierState {
    PhaseOneListening,
    PhaseTwoListening,
    ACCEPT,
    REJECT,
}

#[derive(Debug)]
pub struct GKRVerifier {
    state: GKRVerifierState,
    rng: ThreadRng,
    f1: HashMap<usize, u64>,
    f2: Vec<u64>,
    f3: Vec<u64>,
    g: Vec<u64>,
    p: u64,
    l: usize,
    asserted_sum: u64,
    phase1_verifier: InteractivePMFVerifier,
    phase2_verifier: Option<InteractivePMFVerifier>,
}

impl GKRVerifier {
    pub fn new() -> Self {
        todo!()
    }

    pub fn talk_phase1(&mut self, msgs: &[u64]) -> (bool, u64) {
        todo!()
    }

    pub fn talk_phase2(&mut self, msgs: &[u64]) -> (bool, u64) {
        todo!()
    }

    pub fn _verdict(&mut self) -> bool {
        todo!()
    }

    fn get_randomness_u(&self) -> Vec<u64> {
        todo!()
    }

    fn get_randomness_v(&self) -> Vec<u64> {
        todo!()
    }
}
