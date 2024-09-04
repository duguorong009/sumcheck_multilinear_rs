use std::collections::HashMap;

use crate::{
    gkr::GKR,
    gkr_verifier::{GKRVerifier, GKRVerifierState},
};

/// Change binary form to list of arguments.
///
/// `b`: The binary form in little endian encoding. For example, 0b1011 means g(x0=1, x1=1, x2=0, x3=1)
/// `num_variables`: number of variables
pub fn binary_to_list(mut b: usize, num_variables: usize) -> Vec<usize> {
    let mut lst = vec![0; num_variables];
    let mut i = 0;
    while b != 0 {
        lst[i] = b & 1;
        b >>= 1;
        i += 1;
    }
    lst
}

pub fn precompute(g: Vec<u64>, p: u64) -> Vec<u64> {
    let l = g.len();
    let mut g = vec![0; 1 << l];
    g[0] = 1 - g[0];
    g[1] = g[0];
    for i in 1..l {
        let old_g = g.clone();
        for b in 0..1 << i {
            g[b] = old_g[b] * (1 - g[i]) % p;
            g[b + (1 << i)] = old_g[b] * g[i] % p;
        }
    }
    g
}

/// Split the argument into three parts.
///
/// `arg`: The argument of length 3L
/// `l`: The number of variables
/// `z`: (first L bits), `x`: (second L bits), `y`: (last L bits)
pub fn _three_split(arg: usize, l: usize) -> (usize, usize, usize) {
    let z = arg & ((1 << l) - 1);
    let x = (arg & (((1 << l) - 1) << l)) >> l;
    let y = (arg & (((1 << l) - 1) << (2 * l))) >> (2 * l);
    (z, x, y)
}

///
///
/// `f1`: Sparse polynomial f1(z, x, y) Sparse MVLinear represented by a map of [argument in little endian binary form, evaluation].
/// `l`: number of variables of f3
/// `p`: field size
/// `a_f3`: Bookkeeping table of f3 (where f3 is the multilinear extension of f3)
/// `g`: fixed parameter g of f1
/// return: Bookkeeping table of h_g = sum over y: f1(g, x, y)*f3(y). It has size 2 ** l. It also returns G, which is precompute(g, p), that is useful for phase two.
pub fn initialize_phase_one(
    f1: HashMap<usize, u64>,
    l: usize,
    p: u64,
    a_f3: Vec<u64>,
    g: Vec<u64>,
) -> (Vec<u64>, Vec<u64>) {
    assert!(a_f3.len() == 1 << l);
    assert!(g.len() == l);

    let mut a_hg = vec![0; 1 << l];
    let g = precompute(g, p);

    // rely on sparsity
    for (arg, ev) in f1.into_iter() {
        let (z, x, y) = _three_split(arg, l);
        a_hg[x] += g[z] * ev * a_f3[y] % p;
    }
    (a_hg, g)
}

/// calculate the sum of the GKR.
pub fn sum_of_gkr(a_hg: &[u64], f2: &[u64], p: u64) -> u64 {
    assert!(a_hg.len() == f2.len());
    let mut s = 0;
    for i in 0..a_hg.len() {
        s += a_hg[i] * f2[i] % p;
    }
    s
}

/// phase two
///
/// `f1`: Sparse polynomial f1(z, x, y) Sparse MVLinear represented by a map of [argument in little endian binary form, evaluation].
/// `g`: precompute(g, p), which is outputted in phase one. It has size 2 ** l.
/// `u`: randomness of previous phase sum check protocol. It has size l (#variables in f2, f3).
/// `p`: field size
/// return: Bookkeeping table of f1(g, u, y) over y. It has size 2 ** l.
pub fn initialize_phase_two(f1: HashMap<usize, u64>, g: &[u64], u: &[u64], p: u64) -> Vec<u64> {
    let l = u.len();
    let u = precompute(g.to_vec(), p);
    assert!(u.len() == g.len());
    let mut a_f1 = vec![0; 1 << l];
    for (arg, ev) in f1.into_iter() {
        let (z, x, y) = _three_split(arg, l);
        a_f1[y] = (a_f1[y] + g[z] * u[x] * ev) % p;
    }
    a_f1
}

// type Talker = dyn Fn(&Vec<u64>) -> (bool, u64);

fn talk_process<Talker>(
    mut as_vec: (Vec<u64>, Vec<u64>),
    l: usize,
    p: u64,
    talker: &Talker,
    msg_recorder: &mut Option<Vec<Vec<u64>>>,
) where
    Talker: Fn(&Vec<u64>) -> (bool, u64),
{
    let num_multiplicands = 2;
    for i in 1..=l {
        let mut product_sum: Vec<u64> = vec![0u64.into(); num_multiplicands as usize + 1];
        for b in 0..(1 << (l - i)) {
            for t in 0..=num_multiplicands {
                let mut product = 1;
                for j in 0..num_multiplicands {
                    let a = match j {
                        0 => &as_vec.0.clone(),
                        1 => &as_vec.1.clone(),
                        _ => unreachable!(),
                    };
                    product = product
                        * (((a[b << 1] * ((1 - t as u64) % p))
                            + (a[((b << 1) + 1) as usize] * (t as u64)) % p)
                            % p)
                        % p;
                }
                product_sum[t as usize] = (product_sum[t as usize] + product) % p;
            }
        }

        if let Some(msg_recorder) = msg_recorder {
            msg_recorder.push(product_sum.clone());
        }
        let (result, r) = talker(&product_sum);

        assert!(result);
        for j in 0..num_multiplicands {
            for b in 0..(1 << (l - i)) {
                match j {
                    0 => {
                        as_vec.0[b as usize] = (as_vec.0[(b << 1) as usize] * (1 - r)
                            + as_vec.0[((b << 1) + 1) as usize] * r)
                            % p;
                    }
                    1 => {
                        as_vec.1[b as usize] = (as_vec.1[(b << 1) as usize] * (1 - r)
                            + as_vec.1[((b << 1) + 1) as usize] * r)
                            % p;
                    }
                    _ => unreachable!(),
                }
            }
        }
    }
}

/// Attempt to prove to GKR verifier.
///
/// :param randomGen: add randomness
/// :param A_hg: Bookkeeping table of hg. A_hg will be modified in-place. Do not reuse it!
/// :param gkr: The GKR function
/// :return: randomness, f2(u)
fn talk_to_verifier_phase_one(
    a_hg: &[u64],
    gkr: GKR,
    verifier: GKRVerifier,
    msg_recorder: &mut Option<Vec<Vec<u64>>>,
) -> (Vec<u64>, u64) {
    // sanity check
    let l = gkr.l;
    let p = gkr.p;
    assert!(
        verifier.state == GKRVerifierState::PhaseOneListening,
        "Verifier is not in phase one."
    );
    assert!(a_hg.len() == 1 << l, "Mismatch A_hg size and L");

    let As: (Vec<u64>, Vec<u64>) = (a_hg.clone().to_vec(), gkr.f2.clone());
    talk_process(As, l, p, &verifier.talk_phase1, msg_recorder);

    (verifier.get_randomness_v(), As.1[0])
}

fn talk_to_verifier_phase_two(
    a_f1: &[u64],
    gkr: GKR,
    f2u: u64,
    verifier: GKRVerifier,
    msg_recorder: &mut Option<Vec<Vec<u64>>>,
) {
    let l = gkr.l;
    let p = gkr.p;
    let a_f3_f2u: Vec<u64> = gkr.f3.clone().into_iter().map(|x| x * f2u % p).collect();

    assert!(
        verifier.state == GKRVerifierState::PhaseTwoListening,
        "Verifier is not in phase two."
    );
    assert!(a_f1.len() == 1 << l, "Mismatch A_f1 size and L");

    let As: (Vec<u64>, Vec<u64>) = (a_f1.clone().to_vec(), a_f3_f2u.clone());
    talk_process(As, l, p, &verifier.talk_phase2, msg_recorder);
}

#[derive(Debug)]
pub struct GKRProver {
    gkr: GKR,
}

impl GKRProver {
    pub fn new(gkr: GKR) -> Self {
        Self { gkr }
    }

    pub fn init_and_get_sum(&self, g: &[u64]) -> (Vec<u64>, Vec<u64>, u64) {
        todo!("come back after GKRVerifier is ready")
    }

    pub fn prove_to_verifier(
        &self,
        a_hg: &[u64],
        g: &[u64],
        s: u64,
        verifier: GKRVerifier,
        msg_recorder_phase_1: &mut Option<Vec<Vec<u64>>>,
        msg_recorder_phase_2: &mut Option<Vec<Vec<u64>>>,
    ) {
        todo!("come back after GKRVerifier is ready")
    }
}
