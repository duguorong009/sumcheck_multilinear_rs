use std::collections::HashMap;
use rand::rngs::ThreadRng;

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

