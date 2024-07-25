use std::collections::HashMap;

use num_bigint::BigUint;

use crate::polynomial::{make_MVLinear_constructor, MVLinear};

/// Convert an array to a polynomial where the argument is the binary form of index.
///
/// `data`: Array of size 2^l. If the size of array is not power of 2, out-of-range part will be arbitrary.
/// `field_size`: The size of the finite field that the array value belongs.
///
pub fn extend(data: &[BigUint], field_size: BigUint) -> MVLinear {
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
        let mut sub_poly = gen(vec![(b, data[b].clone())]);
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
pub fn extend_sparse(data: &[(usize, BigUint)], num_var: usize, field_size: BigUint) -> MVLinear {
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
        let sub_poly = gen(vec![(*b, vb.clone())]);
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
    MVLinear::new(
        xs[0].num_variables,
        vec![(0b0000, 1u64.into())],
        xs[0].p.clone(),
    )
}

/// Directly evaluate a polynomial based on multilinear extension. The function takes linear time to the size of the data.
///
/// `data`: The bookkeeping table (where the multilinear extension is based on)
/// `arguments`: Input argument
/// `field_size`: The size of the finite field that the array value belongs.
pub fn evaluate(data: &[BigUint], arguments: &[BigUint], field_size: BigUint) -> BigUint {
    let l = arguments.len();
    let p = field_size;
    assert!(data.len() <= (1 << l));

    let mut a = data.to_vec();
    if a.len() < (1 << l) {
        a.resize(1 << l, 0u64.into());
    }

    for i in 1..l + 1 {
        let r = arguments[i - 1].clone();
        for b in 0..2usize.pow((l - i).try_into().unwrap()) {
            a[b] = (a[b << 1].clone() * (1u64 - r.clone()) + a[(b << 1) + 1].clone() * r.clone())
                % p.clone();
        }
    }
    a[0].clone()
}

/// Sparse version of the function `evaluate`. The function also takes linear time to the size of the data.
///
/// `data`: dictionary indicating a map between binary argument and its value (sparse bookkeeping table)
/// `arguments`: Input argument
/// `field_size`: The size of the finite field that the array value belongs.
pub fn evaluate_sparse(
    data: &[(usize, BigUint)],
    arguments: &[BigUint],
    field_size: BigUint,
) -> BigUint {
    let l = arguments.len();
    let p = field_size;

    let mut dp0 = data.to_vec();
    let mut dp1: HashMap<usize, BigUint> = HashMap::new();
    for i in 0..l {
        let r = arguments[i].clone();
        for (k, v) in dp0 {
            if !dp1.contains_key(&(k >> 1)) {
                dp1.insert(k >> 1, 0u64.into());
            }
            if k & 1 == 0 {
                dp1.insert(
                    k >> 1,
                    (dp1.get(&(k >> 1)).unwrap() + v * (1u64 - r.clone())) % p.clone(),
                );
            } else {
                dp1.insert(
                    k >> 1,
                    (dp1.get(&(k >> 1)).unwrap() + v * r.clone()) % p.clone(),
                );
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
    use num_bigint::RandBigInt;
    use rand::Rng;

    use crate::polynomial::random_prime;

    use super::*;

    #[test]
    fn test_extend() {
        // Create a random number generator
        let mut rng = rand::thread_rng();

        for t in 0..1 {
            let p = random_prime(64);
            let arr: Vec<BigUint> = (0..1024)
                .map(|_| rng.gen_biguint_range(&0u64.into(), &p))
                .collect();
            let poly = extend(&arr, p.clone());
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
            let p = random_prime(64);
            let l = 10;
            let data: HashMap<usize, BigUint> = (0..256)
                .map(|_| {
                    (
                        rng.gen_range(0..1 << l),
                        rng.gen_biguint_range(&0u64.into(), &p),
                    )
                })
                .collect();
            let data_vec: Vec<(usize, BigUint)> = data.clone().into_iter().collect();
            let poly = extend_sparse(&data_vec, l, p.clone());
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
            let p = random_prime(64);
            let l = 8;
            let arr: Vec<BigUint> = (0..1 << l)
                .map(|_| rng.gen_biguint_range(&0u64.into(), &p))
                .collect();
            let poly = extend(&arr, p.clone());
            let args: Vec<BigUint> = (0..l)
                .map(|_| rng.gen_biguint_range(&0u64.into(), &p))
                .collect();
            assert!(poly.eval(&args) == evaluate(&arr, &args, p));
        }
    }

    #[test]
    fn test_evaluate_sparse() {
        let mut rng = rand::thread_rng();
        for _ in 0..1 {
            let p = random_prime(64);
            let l = 9;
            let data: HashMap<usize, BigUint> = (0..1 << 3)
                .map(|_| {
                    (
                        rng.gen_range(0..1 << l),
                        rng.gen_biguint_range(&0u64.into(), &p),
                    )
                })
                .collect();
            let data_vec: Vec<(usize, BigUint)> = data.clone().into_iter().collect();
            let poly = extend_sparse(&data_vec, l, p.clone());
            let args: Vec<BigUint> = (0..l)
                .map(|_| rng.gen_biguint_range(&0u64.into(), &p))
                .collect();
            assert!(poly.eval(&args) == evaluate_sparse(&data_vec, &args, p));
        }
    }
}
