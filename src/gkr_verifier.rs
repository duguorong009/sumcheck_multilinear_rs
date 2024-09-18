use std::collections::HashMap;

use halo2curves::ff::PrimeField;

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
pub struct GKRVerifier<F: PrimeField + Clone> {
    pub(crate) state: GKRVerifierState,
    rng: Option<TrueRandomGen>,
    f1: HashMap<usize, F>,
    f2: Vec<F>,
    f3: Vec<F>,
    g: Vec<F>,
    l: usize,
    pub(crate) asserted_sum: F,
    phase1_verifier: InteractivePMFVerifier<F>,
    phase2_verifier: Option<InteractivePMFVerifier<F>>,
}

impl<F> GKRVerifier<F>
where
    F: PrimeField + Clone,
{
    pub fn new(gkr: GKR<F>, g: Vec<F>, asserted_sum: F, rng: Option<TrueRandomGen>) -> Self {
        Self {
            state: GKRVerifierState::PhaseOneListening,
            rng: rng.clone(),
            f1: gkr.f1.clone(),
            f2: gkr.f2.clone(),
            f3: gkr.f3.clone(),
            g,
            l: gkr.l,
            asserted_sum,
            phase1_verifier: InteractivePMFVerifier::new(
                // DummyPMF::new(2, gkr.l, gkr.p),   // TODO: SHOULD enable this code after handling the inheritance case
                PMF::new(vec![
                    MVLinear::new(gkr.l, vec![]),
                    MVLinear::new(gkr.l, vec![]),
                ]),
                asserted_sum,
                None,
                Some(true),
                rng,
            ),
            phase2_verifier: None,
        }
    }

    pub fn talk_phase1(&mut self, msgs: &[F]) -> (bool, F) {
        if self.state != GKRVerifierState::PhaseOneListening {
            panic!("Verifier is not in phase 1.");
        }

        let (_, r) = self.phase1_verifier.talk(msgs);
        if self.phase1_verifier.convinced {
            let l = self.l;
            self.phase2_verifier = Some(InteractivePMFVerifier::new(
                PMF::new(vec![MVLinear::new(l, vec![]), MVLinear::new(l, vec![])]),
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

    pub fn talk_phase2(&mut self, msgs: &[F]) -> (bool, F) {
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
        let self_f1: Vec<(usize, F)> = self.f1.clone().into_iter().collect();
        let mut args = Vec::with_capacity(self.g.len() + u.len() + v.len());
        args.extend_from_slice(&self.g);
        args.extend_from_slice(&u);
        args.extend_from_slice(&v);
        let m1 = evaluate_sparse(&self_f1, &args);

        let m2 = evaluate(&self.f3, &v) * evaluate(&self.f2, &u);

        let expected = m1 * m2;

        if (self.phase2_verifier.as_mut().unwrap().subclaim().1 - expected) != F::ZERO {
            self.state = GKRVerifierState::REJECT;
            return false;
        }

        self.state = GKRVerifierState::ACCEPT;
        true
    }

    pub fn get_randomness_u(&self) -> Vec<F> {
        if self.state == GKRVerifierState::PhaseOneListening
            || self.state == GKRVerifierState::REJECT
        {
            panic!("Not in correct phase.");
        }
        self.phase1_verifier.points.clone()
    }

    pub fn get_randomness_v(&self) -> Vec<F> {
        if self.state != GKRVerifierState::ACCEPT {
            panic!("Not in correct phase.");
        }
        self.phase2_verifier.as_ref().unwrap().points.clone()
    }
}
