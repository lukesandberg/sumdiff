use crate::token::Token;
use std::mem::MaybeUninit;

/// An implementation of the Myers LCS algorithm for calculating the length of an LCS
/// - http://www.xmailserver.org/diff2.pdf is the original paper that is quite readable.
/// - https://blog.jcoglan.com/2017/02/12/the-myers-diff-algorithm-part-1/ was a secondary resource
///   with more examples.
///
/// This implementation is very true to the paper.
pub fn meyers_lcs_length(left: &[Token], right: &[Token]) -> usize {
    let n = left.len();
    let m = right.len();
    if n == 0 || m == 0 {
        return 0;
    }
    let max_distance = n + m;
    // Save time by not initializing the entire vector
    // SAFETY: we will only access values that have been initialized.  By induction
    // when d==0, on the loop over k onjly accesses `v[k+1]` which is exacrly the value we initialize below.
    // When d>0 and even, k only accesses the odd values for a smaller range of k values and assigns the even ones.
    // When d>0 and odd, k only accesses the even values for a smaller range of k values and assigns the odd ones.
    // In this way we only ever access values of v that were initialized by a prior iteration of the loop.
    let mut v = vec![MaybeUninit::<usize>::uninit(); 2 * max_distance + 1].into_boxed_slice();
    v[max_distance + 1] = MaybeUninit::new(0);
    for d in 0..=max_distance {
        // In the algorithm, k is modeled as a signed integer, but we can use an unsigned integer by offsetting
        // by max distance
        let lower = max_distance - d;
        let upper = max_distance + d;
        // k is our diagonal index
        for k in (lower..=upper).step_by(2) {
            let mut x = unsafe {
                // SAFETY: we only access values of v that have been initialized
                if k == lower || (k != upper && v[k - 1].assume_init() < v[k + 1].assume_init()) {
                    v[k + 1].assume_init()
                } else {
                    v[k - 1].assume_init() + 1
                }
            };

            // The corresponding y value is x - k but our k is offset by max_distance so correct for that
            let mut y = max_distance + x - k;

            // Follow the diagonal as far as we can.
            // Using suffix arrays and LCP arrays, we can skip over common prefixes in O(1) time, but construction is
            // complex.
            while x < n && y < m && left[x] == right[y] {
                x += 1;
                y += 1;
            }

            v[k] = MaybeUninit::new(x);
            if x == n && y == m {
                return (max_distance - d) / 2;
            }
        }
    }
    panic!("Should have found a solution by now");
}

#[cfg(test)]
mod tests {
    use crate::{lex::lex_characters, token::Tokens};

    use super::*;

    #[test]
    fn test_fencepost_both_empty() {
        let mut tokens = Tokens::new();

        let z = lex_characters(&mut tokens, "");
        let length = meyers_lcs_length(&z, &z);
        assert_eq!(length, 0);
    }

    #[test]
    fn test_fencepost_one_empty() {
        let mut tokens = Tokens::new();

        let a = lex_characters(&mut tokens, "a");
        let z = lex_characters(&mut tokens, "");
        let length = meyers_lcs_length(&a, &z);
        assert_eq!(length, 0);
        let length = meyers_lcs_length(&z, &a);
        assert_eq!(length, 0);
    }

    #[test]
    fn test_fencepost_no_overlap() {
        let mut tokens = Tokens::new();

        let a = lex_characters(&mut tokens, "a");
        let z = lex_characters(&mut tokens, "bb");
        let length = meyers_lcs_length(&a, &z);
        assert_eq!(length, 0);
        let length = meyers_lcs_length(&z, &a);
        assert_eq!(length, 0);
    }

    #[test]
    fn test_example() {
        let mut tokens = Tokens::new();

        let abc = lex_characters(&mut tokens, "abc");
        let bac = lex_characters(&mut tokens, "bac");
        let length = meyers_lcs_length(&abc, &bac);
        assert_eq!(length, 2);
    }
    #[test]
    fn test_example_longer() {
        let mut tokens = Tokens::new();
        let abacad = lex_characters(&mut tokens, "abacad");
        let bacada = lex_characters(&mut tokens, "bacada");
        let length = meyers_lcs_length(&abacad, &bacada);
        assert_eq!(length, 5);
    }
}
