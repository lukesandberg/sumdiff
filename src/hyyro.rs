use rug::ops::SubFrom;
use rug::{Assign, Integer};

use crate::{
    lcs_utils::{common_prefix, common_suffix},
    token::{Token, Tokens},
};

/// Computes the length of an LCS between two slices of tokens using the Hyyro algorithm.
/// Bit-Parallel LCS-length Computation Revisited
/// https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=7b1385ba60875b219ce76d5dc0fb343f664c6d6a
///
/// The algorithm leverages bit-parallelism and vector instructions to speed up computation
///
/// This is mostly for benchmarking purposes, as it consistenly loses to the `meyers` algorithm.
pub fn hyyro_lcs_len(tokens: &Tokens, mut left: &[Token], mut right: &[Token]) -> usize {
    let prefix = common_prefix(left, right);
    if prefix > 0 {
        left = &left[prefix..];
        right = &right[prefix..];
    }
    let suffix = common_suffix(left, right);
    if suffix > 0 {
        left = &left[..left.len() - suffix];
        right = &right[..right.len() - suffix];
    }
    let n = left.len();
    let m = right.len();
    let trimmed_length = prefix + suffix;
    if n == 0 || m == 0 {
        return trimmed_length;
    }
    // TODO: arrange for n <= m
    let num_bits = n;
    // For every token, we store a bitset of all the positions it appears in the left string
    let mut matches = vec![Integer::with_capacity(num_bits); tokens.token_upper_bound()];
    for (i, a) in left.iter().enumerate() {
        let pm: &mut Integer = matches.get_mut(*a as usize).unwrap();
        pm.set_bit(i as u32, true);
    }
    let mut v = Integer::with_capacity(num_bits + 1);
    v.set_bit(num_bits.try_into().unwrap(), true);
    v -= 1;
    v.keep_bits_mut(num_bits as u32);
    // use temporaries to avoid re-allocating integers inside the loop
    let mut u = Integer::with_capacity(num_bits);
    // t1 needs an extra bit for the addition to avoid reallocation.
    // we don't actually care about the extra bit
    let mut t1 = Integer::with_capacity(num_bits + 1);
    for b in right {
        let pm = matches.get(*b as usize).unwrap();
        // From the paper
        // u = V' &PM[b[j]]
        // v' = (v' + u) | (v' - u)
        u.assign(&v & pm);
        t1.assign(&v + &u);
        // The addition may have carried a bit that we don't want
        t1.keep_bits_mut(num_bits as u32);
        u.sub_from(&v);
        v.assign(&t1 | &u);
    }
    // confusingly count_zeros only works on negative numbers, so we count ones and subtract from the total to get
    // the number of zeros as required.
    (num_bits - v.count_ones().unwrap() as usize) + trimmed_length
}

#[cfg(test)]
mod tests {

    mod length {
        use crate::{hyyro::hyyro_lcs_len, lex::lex_characters, token::Tokens};

        #[test]
        fn test_fencepost_both_empty() {
            let mut tokens = Tokens::new();

            let z = lex_characters(&mut tokens, "");
            let length = hyyro_lcs_len(&tokens, &z, &z);
            assert_eq!(length, 0);
        }

        #[test]
        fn test_fencepost_one_empty() {
            let mut tokens = Tokens::new();

            let a = lex_characters(&mut tokens, "a");
            let z = lex_characters(&mut tokens, "");
            let length = hyyro_lcs_len(&tokens, &a, &z);
            assert_eq!(length, 0);
            let length = hyyro_lcs_len(&tokens, &z, &a);
            assert_eq!(length, 0);
        }

        #[test]
        fn test_fencepost_no_overlap() {
            let mut tokens = Tokens::new();

            let a = lex_characters(&mut tokens, "a");
            let z = lex_characters(&mut tokens, "bb");
            let length = hyyro_lcs_len(&tokens, &a, &z);
            assert_eq!(length, 0);
            let length = hyyro_lcs_len(&tokens, &z, &a);
            assert_eq!(length, 0);
        }

        #[test]
        fn test_example() {
            let mut tokens = Tokens::new();

            let abc = lex_characters(&mut tokens, "abc");
            let bac = lex_characters(&mut tokens, "bac");
            let length = hyyro_lcs_len(&tokens, &abc, &bac);
            assert_eq!(length, 2);
        }

        #[test]
        fn test_example_longer() {
            let mut tokens = Tokens::new();
            let abacad = lex_characters(&mut tokens, "abacad");
            let bacada = lex_characters(&mut tokens, "bacada");
            let length = hyyro_lcs_len(&tokens, &abacad, &bacada);
            assert_eq!(length, 5);
        }

        #[test]
        fn test_full_diff() {
            let mut tokens = Tokens::new();
            let abcde = lex_characters(&mut tokens, "abcde");
            let fghijk = lex_characters(&mut tokens, "fghijk");
            let length = hyyro_lcs_len(&tokens, &abcde, &fghijk);
            assert_eq!(length, 0);
        }
    }
}
