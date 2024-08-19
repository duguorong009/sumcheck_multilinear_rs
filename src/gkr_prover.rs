use std::collections::HashMap;

use num_bigint::BigUint;
use num_traits::{One, Zero};

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

pub fn precompute(g: Vec<BigUint>, p: BigUint) -> Vec<BigUint> {
    let l = g.len();
    let mut g = vec![BigUint::ZERO; 1 << l];
    g[0] = BigUint::one() - g[0].clone();
    g[1] = g[0].clone();
    for i in 1..l {
        let old_g = g.clone();
        for b in 0..1 << i {
            g[b] = old_g[b].clone() * (BigUint::one() - g[i].clone()) % p.clone();
            g[b + (1 << i)] = old_g[b].clone() * g[i].clone() % p.clone();
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
    f1: HashMap<usize, BigUint>,
    l: usize,
    p: BigUint,
    a_f3: Vec<BigUint>,
    g: Vec<BigUint>,
) -> (Vec<BigUint>, Vec<BigUint>) {
    assert!(a_f3.len() == 1 << l);
    assert!(g.len() == l);

    let mut a_hg = vec![BigUint::ZERO; 1 << l];
    let g = precompute(g, p.clone());

    // rely on sparsity
    for (arg, ev) in f1.into_iter() {
        let (z, x, y) = _three_split(arg, l);
        a_hg[x] += g[z].clone() * ev * a_f3[y].clone() % p.clone();
    }
    (a_hg, g)
}

/// calculate the sum of the GKR.
pub fn sum_of_gkr(a_hg: &[BigUint], f2: &[BigUint], p: BigUint) -> BigUint {
    assert!(a_hg.len() == f2.len());
    let mut s = BigUint::zero();
    for i in 0..a_hg.len() {
        s += a_hg[i].clone() * f2[i].clone() % p.clone();
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
pub fn initialize_phase_two(
    f1: HashMap<usize, BigUint>,
    g: &[BigUint],
    u: &[BigUint],
    p: BigUint,
) -> Vec<BigUint> {
    let l = u.len();
    let u = precompute(g.to_vec(), p.clone());
    assert!(u.len() == g.len());
    let mut a_f1 = vec![BigUint::ZERO; 1 << l];
    for (arg, ev) in f1.into_iter() {
        let (z, x, y) = _three_split(arg, l);
        a_f1[y] = (a_f1[y].clone() + g[z].clone() * u[x].clone() * ev) % p.clone();
    }
    a_f1
}

type Talker = dyn Fn(&Vec<BigUint>) -> (bool, BigUint);

fn talk_process(
    mut as_vec: Vec<Vec<BigUint>>,
    l: usize,
    p: BigUint,
    talker: &Talker,
    msg_recorder: Option<&mut Vec<Vec<BigUint>>>,
) {
    let num_multiplicands = 2;
    for i in 1..=l {
        let mut product_sum: Vec<BigUint> = vec![0u64.into(); num_multiplicands as usize + 1];
        for b in 0..(1 << (l - i)) {
            for t in 0..=num_multiplicands {
                let mut product = BigUint::one();
                for j in 0..num_multiplicands {
                    let a = &as_vec[j].clone();
                    product = product
                        * (((a[(b << 1) as usize].clone() * ((1 - t) % p.clone()))
                            + (a[((b << 1) + 1) as usize].clone() * t) % p.clone())
                            % p.clone())
                        % p.clone();
                }
                product_sum[t as usize] = (product_sum[t as usize].clone() + product) % p.clone();
            }
        }
        todo!("resolve the error in following commented code");
        // if let Some(ref msg_recorder) = msg_recorder {
        //     msg_recorder.push(product_sum.clone());
        // }
        let (result, r) = talker(&product_sum);

        assert!(result);
        for j in 0..num_multiplicands {
            for b in 0..(1 << (l - i)) {
                as_vec[j][b as usize] = (as_vec[j][(b << 1) as usize].clone()
                    * (BigUint::one() - r.clone())
                    + as_vec[j][((b << 1) + 1) as usize].clone() * r.clone())
                    % p.clone();
            }
        }
    }
}
