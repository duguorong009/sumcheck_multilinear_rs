use std::collections::HashMap;

use num_bigint::BigUint;

use crate::polynomial::{make_MVLinear_constructor, MVLinear};

/// Convert an array to a polynomial where the argument is the binary form of index.
///
/// `data`: Array of size 2^l. If the size of array is not power of 2, out-of-range part will be arbitrary.
/// `field_size`: The size of the finite field that the array value belongs.
///
pub fn extend(data: Vec<usize>, field_size: BigUint) -> MVLinear {
    let l = (data.len() as f64).log2().ceil() as usize;
    let p = field_size;
    let gen = make_MVLinear_constructor(l, p.clone());
    let x: Vec<MVLinear> = (0..l)
        .into_iter()
        .map(|i| gen(vec![(1 << i, 1u64.into())]))
        .collect();

    let mut poly_terms: HashMap<usize, BigUint> = (0..2usize.pow(l.try_into().unwrap()))
        .into_iter()
        .map(|i| (i, 0u64.into()))
        .collect();

    for b in 0..data.len() {
        let mut sub_poly = gen(vec![(b, data[b].into())]);
        let xi0 = {
            let mut xi0 = vec![];
            for i in 0..l {
                if b & (1 << i) > 0 {
                    xi0.push(x[i].clone());
                }
            }
            xi0
        };
        sub_poly *= _product1mx(&xi0, 0, xi0.len() - 1);
        for (t, v) in sub_poly.terms {
            poly_terms.insert(t, (poly_terms.get(&t).unwrap() + v) % p.clone());
        }
    }

    gen(poly_terms.into_iter().collect())
}

/// Convert a sparse map to a polynomial where the argument is the binary form of index.
///
/// `data`: sparse map if index < 2^L. If the size of the array is not power of 2, out-of-range part will be arbitrary.
/// `num_var`: number of variables
/// `field_size`: The size of the finite field that the array value belongs.
pub fn extend_sparse(data: Vec<(usize, BigUint)>, num_var: usize, field_size: BigUint) -> MVLinear {
    let l = num_var;
    let p = field_size;
    let gen = make_MVLinear_constructor(l, p.clone());
    let x: Vec<MVLinear> = (0..l)
        .into_iter()
        .map(|i| gen(vec![(1 << i, 1u64.into())]))
        .collect();

    let mut poly_terms: HashMap<usize, BigUint> = (0..2usize.pow(l.try_into().unwrap()))
        .into_iter()
        .map(|i| (i, 0u64.into()))
        .collect();

    for (b, vb) in data {
        let sub_poly = gen(vec![(b, vb)]);
        let xi0 = {
            let mut xi0 = vec![];
            for i in 0..l {
                if b & (1 << i) > 0 {
                    xi0.push(x[i].clone());
                }
            }
            xi0
        };
        let sub_poly = _product1mx(&xi0, 0, xi0.len() - 1) * sub_poly;
        for (t, v) in sub_poly.terms {
            poly_terms.insert(t, (poly_terms.get(&t).unwrap() + v) % p.clone());
        }
    }
    gen(poly_terms.into_iter().collect())
}

/// Divide and conquer algorithm for calculating product of (1 - xi)
fn _product1mx(xs: &[MVLinear], lo: usize, hi: usize) -> MVLinear {
    if lo == hi {
        return 1u64 - xs[lo].clone();
    }
    if hi > lo {
        let left = _product1mx(xs, lo, lo + (hi - lo) / 2);
        let right = _product1mx(xs, lo + (hi - lo) / 2, hi);
        return left * right;
    }
    MVLinear::new(0, vec![(0b0000, 1u64.into())], xs[0].p.clone())
}
