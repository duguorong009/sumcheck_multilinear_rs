use std::collections::HashMap;

use halo2curves::ff::PrimeField;

use crate::{
    gkr::GKR,
    gkr_verifier::{GKRVerifier, GKRVerifierState},
};

/// Change binary form to list of arguments.
///
/// `b`: The binary form in little endian encoding. For example, 0b1011 means g(x0=1, x1=1, x2=0, x3=1)
/// `num_variables`: number of variables
pub fn binary_to_list<F>(b: F, num_variables: usize) -> Vec<F>
where
    F: PrimeField + Clone,
{
    let bits = to_bits(b.to_repr().as_ref());
    let sliced_bits = bits[..num_variables].to_vec();
    let vec: Vec<F> = sliced_bits.iter().map(|&x| F::from(u64::from(x))).collect();
    vec
}

/// Converts given bytes to the bits.
pub fn to_bits(num: &[u8]) -> Vec<bool> {
    let len = num.len() * 8;
    let mut bits = Vec::new();
    for i in 0..len {
        let bit = num[i / 8] & (1 << (i % 8)) != 0;
        bits.push(bit);
    }
    bits
}

pub fn precompute<F>(g: Vec<F>) -> Vec<F>
where
    F: PrimeField + Clone,
{
    let l = g.len();
    let mut G: Vec<F> = vec![F::ZERO; 1 << l];
    G[0] = F::ONE - g[0];
    G[1] = g[0];
    for i in 1..l {
        let old_G = G.clone();
        for b in 0..1 << i {
            G[b] = old_G[b] * (F::ONE - g[i]);
            G[b + (1 << i)] = old_G[b] * g[i];
        }
    }
    G
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
pub fn initialize_phase_one<F>(
    f1: HashMap<usize, F>,
    l: usize,
    a_f3: Vec<F>,
    g: Vec<F>,
) -> (Vec<F>, Vec<F>)
where
    F: PrimeField + Clone,
{
    assert!(a_f3.len() == 1 << l);
    assert!(g.len() == l);

    let mut a_hg = vec![F::ZERO; 1 << l];
    let G = precompute(g);

    // rely on sparsity
    for (arg, ev) in f1.into_iter() {
        let (z, x, y) = _three_split(arg, l);
        a_hg[x] += G[z] * ev * a_f3[y];
    }
    (a_hg, G)
}

/// calculate the sum of the GKR.
pub fn sum_of_gkr<F>(a_hg: &[F], f2: &[F]) -> F
where
    F: PrimeField + Clone,
{
    assert!(a_hg.len() == f2.len());
    let mut s = F::ZERO;
    for i in 0..a_hg.len() {
        s += a_hg[i] * f2[i];
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
pub fn initialize_phase_two<F>(f1: &HashMap<usize, F>, g: &[F], u: &[F]) -> Vec<F>
where
    F: PrimeField + Clone,
{
    let l = u.len();
    let U = precompute(u.to_vec());
    assert!(U.len() == g.len());
    let mut a_f1 = vec![F::ZERO; 1 << l];
    for (arg, ev) in f1.into_iter() {
        let (z, x, y) = _three_split(*arg, l);
        a_f1[y] += g[z] * U[x] * ev;
    }
    a_f1
}

// type Talker = dyn Fn(&Vec<u64>) -> (bool, u64);

fn talk_process<F, Talker>(
    mut as_vec: (Vec<F>, Vec<F>),
    l: usize,
    talker: &mut Talker,
    msg_recorder: &mut Option<Vec<Vec<F>>>,
) where
    F: PrimeField + Clone,
    Talker: FnMut(&Vec<F>) -> (bool, F),
{
    let num_multiplicands: usize = 2;
    for i in 1..=l {
        let mut product_sum: Vec<F> = vec![F::ZERO; num_multiplicands + 1];
        for b in 0..(1 << (l - i)) {
            for t in 0..=num_multiplicands {
                let mut product = F::ONE;
                for j in 0..num_multiplicands {
                    let a = match j {
                        0 => &as_vec.0.clone(),
                        1 => &as_vec.1.clone(),
                        _ => unreachable!(),
                    };
                    product *= (a[b << 1] * (F::ONE - F::from_u128(t as u128)))
                            + (a[(b << 1) + 1] * F::from_u128(t as u128));
                }
                product_sum[t] += product;
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
                        as_vec.0[b as usize] = as_vec.0[(b << 1) as usize] * (F::ONE - r)
                            + as_vec.0[((b << 1) + 1) as usize] * r;
                    }
                    1 => {
                        as_vec.1[b as usize] = as_vec.1[(b << 1) as usize] * (F::ONE - r)
                            + as_vec.1[((b << 1) + 1) as usize] * r;
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
fn talk_to_verifier_phase_one<F>(
    a_hg: &[F],
    gkr: GKR<F>,
    verifier: &mut GKRVerifier<F>,
    msg_recorder: &mut Option<Vec<Vec<F>>>,
) -> (Vec<F>, F)
where
    F: PrimeField + Clone,
{
    // sanity check
    let l = gkr.l;
    assert!(
        verifier.state == GKRVerifierState::PhaseOneListening,
        "Verifier is not in phase one."
    );
    assert!(a_hg.len() == 1 << l, "Mismatch A_hg size and L");

    let r#as: (Vec<F>, Vec<F>) = (a_hg.to_vec(), gkr.f2.clone());
    talk_process(
        r#as.clone(),
        l,
        &mut |arg| verifier.talk_phase1(arg),
        msg_recorder,
    );

    (verifier.get_randomness_v(), r#as.1[0])
}

fn talk_to_verifier_phase_two<F>(
    a_f1: &[F],
    gkr: GKR<F>,
    f2u: F,
    verifier: &mut GKRVerifier<F>,
    msg_recorder: &mut Option<Vec<Vec<F>>>,
) where
    F: PrimeField + Clone,
{
    let l = gkr.l;
    let a_f3_f2u: Vec<F> = gkr.f3.clone().into_iter().map(|x| x * f2u).collect();

    assert!(
        verifier.state == GKRVerifierState::PhaseTwoListening,
        "Verifier is not in phase two."
    );
    assert!(a_f1.len() == 1 << l, "Mismatch A_f1 size and L");

    let r#as: (Vec<F>, Vec<F>) = (a_f1.to_vec(), a_f3_f2u.clone());
    talk_process(r#as, l, &mut |arg| verifier.talk_phase2(arg), msg_recorder);
}

#[derive(Debug)]
pub struct GKRProver<F: PrimeField + Clone> {
    gkr: GKR<F>,
}

impl<F> GKRProver<F>
where
    F: PrimeField + Clone,
{
    pub fn new(gkr: GKR<F>) -> Self {
        Self { gkr }
    }

    /// :param g: fixed g
    /// :return: Bookkeeping table h_g, G: precompute cache, sum
    pub fn init_and_get_sum(&self, g: &[F]) -> (Vec<F>, Vec<F>, F) {
        assert!(g.len() == self.gkr.l, "Size of g is incorrect");
        let (a_hg, g) = initialize_phase_one(
            self.gkr.f1.clone(),
            self.gkr.l,
            self.gkr.f3.clone(),
            g.to_vec(),
        );
        let s = sum_of_gkr(&a_hg, &self.gkr.f2);
        (a_hg, g, s)
    }

    /// :param A_hg: bookkeeping table h_g
    /// :param G: precompute cache
    /// :param s: sum
    /// :param verifier: GKR verifier
    pub fn prove_to_verifier(
        &self,
        a_hg: &[F],
        g: &[F],
        s: F,
        verifier: &mut GKRVerifier<F>,
        msg_recorder_phase_1: &mut Option<Vec<Vec<F>>>,
        msg_recorder_phase_2: &mut Option<Vec<Vec<F>>>,
    ) {
        assert!(verifier.asserted_sum == s, "Asserted sum mismatch");

        let (u, f2u) =
            talk_to_verifier_phase_one(a_hg, self.gkr.clone(), verifier, msg_recorder_phase_1);

        assert!(
            verifier.state == GKRVerifierState::PhaseTwoListening,
            "Verifier does not accept phase 1 proof"
        );

        let a_f1 = initialize_phase_two(&self.gkr.f1, g, &u);
        talk_to_verifier_phase_two(&a_f1, self.gkr.clone(), f2u, verifier, msg_recorder_phase_2);
    }
}

#[cfg(test)]
mod tests {
    use halo2curves::{bn256::Fr, ff::Field};
    use rand::{rngs::OsRng, Rng};

    use crate::{
        multilinear_extension::{evaluate, extend_sparse},
        polynomial::{random_mvlinear, random_prime, MVLinear},
    };

    use super::*;

    fn generate_random_f1<F>(l: usize) -> HashMap<usize, F>
    where
        F: PrimeField + Clone,
    {
        let mut rng = rand::thread_rng();
        let n = ((1 << (3 * l)) as f64).sqrt() as usize;
        let mut ans = HashMap::new();
        for _ in 0..n {
            let term = rng.gen_range(0..(1 << 3 * l));
            let ev = F::random(OsRng);
            ans.insert(term, ev);
        }
        ans
    }

    fn random_gkr<F>(l: usize) -> GKR<F>
    where
        F: PrimeField + Clone,
    {
        let f1 = generate_random_f1(l);
        let f2 = (0..(1 << l)).map(|_| F::random(OsRng)).collect();
        let f3 = (0..(1 << l)).map(|_| F::random(OsRng)).collect();
        GKR { l, f1, f2, f3 }
    }

    /// :return: A bookkeeping table where the index is the binary form of argument of polynomial and value is the
    /// evaluated value; the sum
    fn calculate_bookkeeping_table<F>(poly: MVLinear<F>) -> (Vec<F>, F)
    where
        F: PrimeField + Clone,
    {
        let mut a: Vec<F> = vec![F::ZERO; 1 << poly.num_variables];
        let mut s = F::ZERO;
        for p in 0..(1 << poly.num_variables) {
            let p_f = F::from_u128(p as u128);
            a[p as usize] = poly.eval(&binary_to_list(p_f, poly.num_variables));
            s = s + a[p as usize];
        }
        (a, s)
    }

    #[test]
    fn test_initialize_phase_one_two() {
        let l = 5;
        println!(
            "Testing GKR Prover Bookkeeping table generator functions... Use L = {l}, p = ???"
        );

        // generate random sparse f1, random f3, g
        let d_f1 = generate_random_f1::<Fr>(l);
        let f3 = random_mvlinear::<Fr>(l);
        let g: Vec<Fr> = (0..l).map(|_| Fr::random(OsRng)).collect();

        // get bookkeeping table for f3
        let (a_f3, _) = calculate_bookkeeping_table(f3);

        // get poly form for f1 (only for test checking)
        let mut a_f1 = vec![Fr::zero(); 1 << (3 * l)];
        for (k, v) in d_f1.iter() {
            a_f1[*k] = *v;
        }

        let d_f1_: Vec<(usize, Fr)> = d_f1.clone().into_iter().collect();
        let f1 = extend_sparse::<Fr>(&d_f1_, 3 * l);

        let f1_fix_g = f1.eval_part(&g);
        assert!(f1_fix_g.num_variables == 2 * l);

        let mut a_hg_expected: Vec<Fr> = vec![Fr::zero(); 1 << l];
        for i in 0..(1 << (2 * l)) {
            let x = i & ((1 << l) - 1);
            let y = (i & (((1 << l) - 1) << l)) >> l;
            a_hg_expected[x] = (a_hg_expected[x] + f1_fix_g.eval_bin(i) * a_f3[y]);
        }
        let (a_hg_actual, G) = initialize_phase_one(d_f1.clone(), l, a_f3, g);
        for i in 0..(1 << l) {
            assert!(a_hg_expected[i] == a_hg_actual[i]);
        }
        println!("PASS: initialize_PhaseOne");

        // phase 2
        let u: Vec<Fr> = (0..l).map(|_| Fr::random(OsRng)).collect();
        let f1_fix_gu = f1_fix_g.eval_part(&u);
        assert!(f1_fix_gu.num_variables == l);

        let mut a_f1_expected: Vec<Fr> = vec![Fr::zero(); 1 << l];
        for i in 0..(1 << l) {
            let y = i & ((1 << l) - 1);
            a_f1_expected[y] = f1_fix_gu.eval_bin(y);
        }

        let a_f1_actual = initialize_phase_two(&d_f1, &G, &u);
        for i in 0..(1 << l) {
            assert!(a_f1_expected[i] == a_f1_actual[i]);
        }
        println!("PASS: initialize_PhaseTwo");
    }

    #[test]
    fn test_completeness_sanity() {
        println!("Test Completeness of GKR protocol (test individual functions");
        let L = 7;
        let gkr = random_gkr(L);
        let g: Vec<Fr> = (0..L).map(|_| Fr::random(OsRng)).collect();

        let (A_hg, G) = initialize_phase_one(gkr.f1.clone(), L, gkr.f3.clone(), g.clone());
        let s = sum_of_gkr(&A_hg, &gkr.f2);
        let mut v = GKRVerifier::new(gkr.clone(), g, s, None);
        assert!(
            v.state == GKRVerifierState::PhaseOneListening,
            "Wrong verifier state"
        );

        let (u, f2u) = talk_to_verifier_phase_one(&A_hg, gkr.clone(), &mut v, &mut None);
        assert!(u.len() == L, "wrong randomness size");
        assert!(
            v.state == GKRVerifierState::PhaseTwoListening,
            "verifier should be in phase two"
        );
        assert!(
            f2u == evaluate(&gkr.f2, &u),
            "f2(u) returned by talk_to_verifier_phase_one is incorrect"
        );

        println!("initialize_phase_one, sum_of_gkr, talk_to_verifier_phase_one looks good. ");
        let A_f1 = initialize_phase_two(&gkr.f1, &G, &u);
        talk_to_verifier_phase_two(&A_f1, gkr.clone(), f2u, &mut v, &mut None);
        assert!(
            v.state == GKRVerifierState::ACCEPT,
            "verifier does not accept this proof. "
        );
        println!("initialize_phase_two, talk_to_verifier_phase_two looks good. ");
        println!("Completeness test PASS!");
    }

    #[test]
    fn test_protocol_comprehensive() {
        let num_tests = 10;
        let Ls = vec![10, 11, 12];
        println!("Performing GKR interactive protocol comprehensive test...");
        for i in 0..num_tests {
            let p = random_prime(32);
            let L = Ls[i % Ls.len()];
            println!("test #{}: |g|=|x|=|y| = {L}: TESTING", i + 1);
            let gkr = random_gkr(L);
            let g: Vec<Fr> = (0..L).map(|_| Fr::random(OsRng)).collect();

            let pv = GKRProver::new(gkr.clone());
            let (A_hg, G, s) = pv.init_and_get_sum(&g);

            let mut v = GKRVerifier::new(gkr, g, s, None);
            assert!(
                v.state == GKRVerifierState::PhaseOneListening,
                "Verifier sanity check failed"
            );
            pv.prove_to_verifier(&A_hg, &G, s, &mut v, &mut None, &mut None);
            assert!(v.state == GKRVerifierState::ACCEPT);
            println!("PASS");
        }
        println!("ALL GKR interactive protocol comprehensive tests PASSED!");
    }
}
