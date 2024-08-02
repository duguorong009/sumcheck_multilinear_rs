use num_bigint::BigUint;
use num_traits::One;

/// Change binary form to list of arguments.
///
/// `b`: The binary form in little endian encoding. For example, 0b1011 means g(x0=1, x1=1, x2=0, x3=1)
/// `num_variables`: number of variables
fn binary_to_list(mut b: usize, num_variables: usize) -> Vec<usize> {
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
pub fn _three_split(arg: BigUint, l: usize) -> (BigUint, BigUint, BigUint) {
    let z = arg.clone() & ((BigUint::one() << l) - BigUint::one());
    let x = (arg.clone() & (((BigUint::one() << l) - BigUint::one()) << l)) >> l;
    let y = (arg & (((BigUint::one() << l) - BigUint::one()) << (2 * l))) >> (2 * l);
    (z, x, y)
}
