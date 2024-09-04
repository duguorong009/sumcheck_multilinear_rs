use std::collections::HashMap;

use crate::{
    gkr::GKR,
    ip_pmf_verifier::{InteractivePMFVerifier, TrueRandomGen},
    multilinear_extension::{evaluate, evaluate_sparse},
    pmf::PMF,
    polynomial::MVLinear,
};

#[derive(Debug, PartialEq, PartialOrd)]
pub enum GKRVerifierState {
    PhaseOneListening,
    PhaseTwoListening,
    ACCEPT,
    REJECT,
}

#[derive(Debug)]
pub struct GKRVerifier {
    pub(crate) state: GKRVerifierState,
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
                PMF::new(vec![
                    MVLinear::new(gkr.l, vec![], gkr.p),
                    MVLinear::new(gkr.l, vec![], gkr.p),
                ]),
                asserted_sum,
                None,
                Some(true),
                rng,
            ),
            phase2_verifier: None,
        }
    }

    pub fn talk_phase1(&mut self, msgs: &[u64]) -> (bool, u64) {
        if self.state != GKRVerifierState::PhaseOneListening {
            panic!("Verifier is not in phase 1.");
        }

        let (_, r) = self.phase1_verifier.talk(msgs);
        if self.phase1_verifier.convinced {
            let l = self.l;
            self.phase2_verifier = Some(InteractivePMFVerifier::new(
                PMF::new(vec![
                    MVLinear::new(l, vec![], self.p),
                    MVLinear::new(l, vec![], self.p),
                ]),
                self.phase1_verifier.subclaim().1,
                None,
                Some(true),
                self.rng.clone(),
            ));
            self.state = GKRVerifierState::PhaseTwoListening;
            return (true, r);
        }

        if !self.phase1_verifier.active && !self.phase1_verifier.convinced {
            self.state = GKRVerifierState::REJECT;
            return (false, r);
        }

        (true, r)
    }

    pub fn talk_phase2(&mut self, msgs: &[u64]) -> (bool, u64) {
        if self.state != GKRVerifierState::PhaseTwoListening {
            panic!("Verifier is not in phase 2.");
        }

        let (_, r) = self.phase2_verifier.as_mut().unwrap().talk(msgs);
        if self.phase2_verifier.as_mut().unwrap().convinced {
            return (self._verdict(), r);
        }

        if !self.phase2_verifier.as_mut().unwrap().active
            && !self.phase2_verifier.as_mut().unwrap().convinced
        {
            self.state = GKRVerifierState::REJECT;
            return (false, r);
        }

        (true, r)
    }

    /// Verify the sub claim of verifier 2, using the u from sub claim 1 and v from sub claim 2.
    /// This requires three polynomial evaluation.
    pub fn _verdict(&mut self) -> bool {
        if self.state != GKRVerifierState::PhaseTwoListening {
            panic!("Verifier is not in phase 2.");
        }
        if !self.phase2_verifier.as_mut().unwrap().convinced {
            panic!("Phase 2 verifier is not convinced.");
        }
        let u = self.phase1_verifier.subclaim().0;
        let v = self.phase2_verifier.as_mut().unwrap().subclaim().0;

        // verify phase 2 verifier's claim
        let self_f1: Vec<(usize, u64)> = self.f1.clone().into_iter().collect();
        let mut args = Vec::with_capacity(self.g.len() + u.len() + v.len());
        args.extend_from_slice(&self.g);
        args.extend_from_slice(&u);
        args.extend_from_slice(&v);
        let m1 = evaluate_sparse(&self_f1, &args, self.p);

        let m2 = evaluate(&self.f3, &v, self.p) * evaluate(&self.f2, &u, self.p) % self.p;

        let expected = m1 * m2 % self.p;

        if (self.phase2_verifier.as_mut().unwrap().subclaim().1 - expected) % self.p != 0 {
            self.state = GKRVerifierState::REJECT;
            return false;
        }

        self.state = GKRVerifierState::ACCEPT;
        true
    }

    fn get_randomness_u(&self) -> Vec<u64> {
        if self.state == GKRVerifierState::PhaseOneListening
            || self.state == GKRVerifierState::REJECT
        {
            panic!("Not in correct phase.");
        }
        self.phase1_verifier.points.clone()
    }

    pub fn get_randomness_v(&self) -> Vec<u64> {
        if self.state != GKRVerifierState::ACCEPT {
            panic!("Not in correct phase.");
        }
        self.phase2_verifier.as_ref().unwrap().points.clone()
    }
}
