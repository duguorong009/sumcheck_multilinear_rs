use std::collections::HashMap;

use crate::polynomial::{make_mvlinear_constructor, MVLinear};

/// Convert an array to a polynomial where the argument is the binary form of index.
///
/// `data`: Array of size 2^l. If the size of array is not power of 2, out-of-range part will be arbitrary.
/// `field_size`: The size of the finite field that the array value belongs.
///
pub fn extend(data: &[u64], field_size: u64) -> MVLinear {
    let l = (data.len() as f64).log2().ceil() as usize;
    let p = field_size;
    let gen = make_mvlinear_constructor(l, p);
    let x: Vec<MVLinear> = (0..l).map(|i| gen(vec![(1 << i, 1u64.into())])).collect();

    let mut poly_terms: HashMap<usize, u64> = (0..2usize.pow(l.try_into().unwrap()))
        .map(|i| (i, 0u64.into()))
        .collect();

    for b in 0..data.len() {
        let mut sub_poly = gen(vec![(b, data[b])]);
        let xi0 = {
            let mut xi0 = vec![];
            for i in 0..l {
                if (b >> i) & 1 == 0 {
                    xi0.push(x[i].clone());
                }
            }
            xi0
        };
        sub_poly *= _product1mx(&xi0, 0, xi0.len() - 1);
        for (t, v) in sub_poly.terms {
            poly_terms.insert(t, (poly_terms.get(&t).unwrap() + v) % p);
        }
    }

    gen(poly_terms.into_iter().collect())
}

/// Convert a sparse map to a polynomial where the argument is the binary form of index.
///
/// `data`: sparse map if index < 2^L. If the size of the array is not power of 2, out-of-range part will be arbitrary.
/// `num_var`: number of variables
/// `field_size`: The size of the finite field that the array value belongs.
pub fn extend_sparse(data: &[(usize, u64)], num_var: usize, field_size: u64) -> MVLinear {
    let l = num_var;
    let p = field_size;
    let gen = make_mvlinear_constructor(l, p);
    let x: Vec<MVLinear> = (0..l).map(|i| gen(vec![(1 << i, 1u64.into())])).collect();

    let mut poly_terms: HashMap<usize, u64> = (0..2usize.pow(l.try_into().unwrap()))
        .map(|i| (i, 0u64.into()))
        .collect();

    for (b, vb) in data {
        let sub_poly = gen(vec![(*b, *vb)]);
        let xi0 = {
            let mut xi0 = vec![];
            for i in 0..l {
                if (b >> i) & 1 == 0 {
                    xi0.push(x[i].clone());
                }
            }
            xi0
        };
        let sub_poly = _product1mx(&xi0, 0, xi0.len() - 1) * sub_poly;
        for (t, v) in sub_poly.terms {
            poly_terms.insert(t, (poly_terms.get(&t).unwrap() + v) % p);
        }
    }
    gen(poly_terms.into_iter().collect())
}

/// Divide and conquer algorithm for calculating product of (1 - xi)
///
/// **NOTE**: When `xs` is empty, the result will be 1.
/// This should be handled by the caller.
fn _product1mx(xs: &[MVLinear], lo: usize, hi: usize) -> MVLinear {
    assert!(hi >= lo);
    assert!(!xs.is_empty());

    if lo == hi {
        return 1u64 - xs[lo].clone();
    }

    let left = _product1mx(xs, lo, lo + (hi - lo) / 2);
    let right = _product1mx(xs, lo + (hi - lo) / 2 + 1, hi);
    left * right
}

/// Directly evaluate a polynomial based on multilinear extension. The function takes linear time to the size of the data.
///
/// `data`: The bookkeeping table (where the multilinear extension is based on)
/// `arguments`: Input argument
/// `field_size`: The size of the finite field that the array value belongs.
pub fn evaluate(data: &[u64], arguments: &[u64], field_size: u64) -> u64 {
    let l = arguments.len();
    let p = field_size;
    assert!(data.len() <= (1 << l));

    let mut a = data.to_vec();
    if a.len() < (1 << l) {
        a.resize(1 << l, 0u64.into());
    }

    for i in 1..l + 1 {
        let r = arguments[i - 1];
        for b in 0..2usize.pow((l - i).try_into().unwrap()) {
            a[b] = (a[b << 1] * (1u64 - r) + a[(b << 1) + 1] * r) % p;
        }
    }
    a[0]
}

/// Sparse version of the function `evaluate`. The function also takes linear time to the size of the data.
///
/// `data`: dictionary indicating a map between binary argument and its value (sparse bookkeeping table)
/// `arguments`: Input argument
/// `field_size`: The size of the finite field that the array value belongs.
pub fn evaluate_sparse(data: &[(usize, u64)], arguments: &[u64], field_size: u64) -> u64 {
    let l = arguments.len();
    let p = field_size;

    let mut dp0 = data.to_vec();
    let mut dp1: HashMap<usize, u64> = HashMap::new();
    for i in 0..l {
        let r = arguments[i];
        for (k, v) in dp0 {
            dp1.entry(k >> 1).or_insert_with(|| 0u64.into());
            if k & 1 == 0 {
                dp1.insert(k >> 1, (dp1.get(&(k >> 1)).unwrap() + v * (1u64 - r)) % p);
            } else {
                dp1.insert(k >> 1, (dp1.get(&(k >> 1)).unwrap() + v * r) % p);
            }
        }
        dp0 = dp1.into_iter().collect();
        dp1 = HashMap::new();
    }
    dp0.into_iter()
        .find(|(k, _)| *k == 0)
        .map(|(_, v)| v)
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use rand::Rng;

    use crate::polynomial::random_prime;

    use super::*;

    #[test]
    fn test_extend() {
        // Create a random number generator
        let mut rng = rand::thread_rng();

        for t in 0..1 {
            let p = random_prime(32);
            let arr: Vec<u64> = (0..1024).map(|_| rng.gen_range(0..p)).collect();
            let poly = extend(&arr, p);
            for _ in 0..poly.num_variables.pow(2) {
                let i = rng.gen_range(0..arr.len());
                assert!(arr[i] == poly.eval_bin(i));
            }
            println!("Test #{} passed", t);
        }
    }

    #[test]
    fn test_extend_sparse() {
        let mut rng = rand::thread_rng();

        for t in 0..1 {
            let p = random_prime(32);
            let l = 10;
            let data: HashMap<usize, u64> = (0..256)
                .map(|_| (rng.gen_range(0..1 << l), rng.gen_range(0..p)))
                .collect();
            let data_vec: Vec<(usize, u64)> = data.clone().into_iter().collect();
            let poly = extend_sparse(&data_vec, l, p);
            for k in 0..1 << l {
                let zero = 0u64.into();
                let expected = data.get(&k).unwrap_or(&zero);
                let actual = poly.eval_bin(k);
                assert!(*expected == actual);
            }
            println!("Test #{} passed", t);
        }
    }

    #[test]
    fn test_evaluate() {
        let mut rng = rand::thread_rng();

        for _ in 0..1 {
            let p = random_prime(32);
            let l = 8;
            let arr: Vec<u64> = (0..1 << l).map(|_| rng.gen_range(0..p)).collect();
            let poly = extend(&arr, p);
            let args: Vec<u64> = (0..l).map(|_| rng.gen_range(0..p)).collect();
            assert!(poly.eval(&args) == evaluate(&arr, &args, p));
        }
    }

    #[test]
    fn test_evaluate_sparse() {
        let mut rng = rand::thread_rng();
        for _ in 0..1 {
            let p = random_prime(32);
            let l = 9;
            let data: HashMap<usize, u64> = (0..1 << 3)
                .map(|_| (rng.gen_range(0..1 << l), rng.gen_range(0..p)))
                .collect();
            let data_vec: Vec<(usize, u64)> = data.clone().into_iter().collect();
            let poly = extend_sparse(&data_vec, l, p);
            let args: Vec<u64> = (0..l).map(|_| rng.gen_range(0..p)).collect();
            assert!(poly.eval(&args) == evaluate_sparse(&data_vec, &args, p));
        }
    }
}
