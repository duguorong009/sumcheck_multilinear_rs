use rand::rngs::ThreadRng;
use std::collections::HashMap;

use crate::{
    gkr::GKR,
    ip_pmf_verifier::{InteractivePMFVerifier, TrueRandomGen},
    pmf::{DummyPMF, PMF}, polynomial::MVLinear,
};

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
    rng: Option<TrueRandomGen>,
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
    pub fn new(gkr: GKR, g: Vec<u64>, asserted_sum: u64, rng: Option<TrueRandomGen>) -> Self {
        Self {
            state: GKRVerifierState::PhaseOneListening,
            rng: rng.clone(),
            f1: gkr.f1.clone(),
            f2: gkr.f2.clone(),
            f3: gkr.f3.clone(),
            g,
            p: gkr.p,
            l: gkr.l,
            asserted_sum,
            phase1_verifier: InteractivePMFVerifier::new(
                // DummyPMF::new(2, gkr.l, gkr.p),   // TODO: SHOULD enable this code after handling the inheritance case
                PMF::new(vec![MVLinear::new(gkr.l, vec![], gkr.p), MVLinear::new(gkr.l, vec![], gkr.p)]),
                asserted_sum,
                None,
                Some(true),
                rng,
            ),
            phase2_verifier: None,
        }
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
