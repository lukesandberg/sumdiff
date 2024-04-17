use assume::assume;

use crate::{
    lcs_utils::{common_prefix, common_suffix},
    token::Token,
};

/// An implementation of the Wu-Manber-Meyers algorithm for computing the length of the longest common subsequence
/// See https://publications.mpi-cbg.de/Wu_1990_6334.pdf
///
/// This algorithm is faster than the Meyers algorithm with worst cast O(NP) time and exectped O(N+PD) time
/// Where `D` is the edit distance and `P` is number of deletions in an edit script.
///
pub fn wu_manber_meyers_lcs_length(mut left: &[Token], mut right: &[Token]) -> usize {
    // Start by trimming suffixes and prefixes to reduce the size of the superlinear algorithm.
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
    // The algorithm requires m>=n so swap if needed.
    (if n > m {
        do_wu_manber_meyers_lcs_length(m, right, n, left)
    } else {
        do_wu_manber_meyers_lcs_length(n, left, m, right)
    }) + trimmed_length
}

fn do_wu_manber_meyers_lcs_length(n: usize, left: &[Token], m: usize, right: &[Token]) -> usize {
    assume!(unsafe: n>0);
    assume!(unsafe: m>0);
    assume!(unsafe: n<=m);
    assume!(unsafe: n == left.len());
    assume!(unsafe: m == right.len());
    // The largest possiple `p` value is `n` and we need to be able to access a sentinel value at offset-p-1
    // so offset has to be at least n+1;
    let offset = n + 1;
    // Unlike the paper we initialize to 0 instead of -1.
    // This is because rust0 based indexing and usize cannot take on negative numbers, so
    // we store 1-based indices in here and simply adjust them in `snake` to 0-based indices.
    // The largest value we access is when `p=n``, in which case we access `p+delta_index+1`
    // in the second loop.
    // By simplifying that we see that the largest value is `n+m+2`, so an array containing that
    // value must have size n+m+3
    let mut fp = vec![0usize; n + m + 3];
    // This is the diagonal that (n,m) is on.  A simplification of `offset + (m-n)`
    let delta_index = m + 1;

    // At this point we know that there is no snake along diagonal 0, so at least some calls to
    // `snake` in the first iteration will be worthless.  However, we do access other diagonals
    // based on delta so skipping/simplifying the first iteration is non-trivial unless delta is 0.
    // Don't bother.
    let mut p: usize = 0;
    loop {
        // The paper defines this as [-p, delta-1], which is [-p,delta) as a half exclusive interval.
        // NOTE: `p` is bounded by `n` because it tracks the number of deletations and we cannot
        // delete more than `n` tokens.
        assume!(unsafe: p<=n);
        let mut lower = 0;
        for k_index in (offset - p)..delta_index {
            // Because p <=n and offset = n+1, we know that k_index >= 1, make sure the compiler can tell.
            assume!(unsafe: k_index >= 1);
            let y = std::cmp::max(fp[k_index - 1] + 1, fp[k_index + 1]);
            let x = offset + y - k_index;
            lower = snake(x, y, n, left, m, right);
            fp[k_index] = lower;
        }

        let mut upper = 0;
        // The paper defines this as [delta+p-> delta+1]
        for k_index in ((delta_index + 1)..=(delta_index + p)).rev() {
            // NOTE: again that `p` is bounded by `n`, so `k_index` is bounded by `n+m+1`
            assume!(unsafe: k_index <= n + m + 1);
            let y = std::cmp::max(fp[k_index - 1] + 1, fp[k_index + 1]);
            let x = offset + y - k_index;
            upper = snake(x, y, n, left, m, right);
            fp[k_index] = upper;
        }
        let y = std::cmp::max(lower + 1, upper);
        let x = offset + y - delta_index;
        let fp_delta = snake(x, y, n, left, m, right);
        // `snake` returns 1 based indices and we are looking for when we hit the end of the right string
        // so our exit condition is `fp_delta == m+1` which is conveniently equivalent to
        // `fp_delta == delta_index` due to how the offsets are set up.
        if fp_delta == delta_index {
            // The LCS length is (max_distance - edit_distance) / 2
            // max_distance = n + m
            // edit_distance = delta + 2p
            // delta = m - n
            // so by the power of algebra we get n-p
            return n - p;
        }
        fp[delta_index] = fp_delta;
        p += 1;
    }
}

fn snake(mut x: usize, mut y: usize, n: usize, left: &[Token], m: usize, right: &[Token]) -> usize {
    // x and y are stored as 1 based indices
    assume!(unsafe: x > 0);
    assume!(unsafe: y > 0);
    x -= 1;
    y -= 1;
    while x < n && y < m && left[x] == right[y] {
        x += 1;
        y += 1;
    }
    y + 1
}

#[cfg(test)]
mod tests {

    mod length {
        use crate::{
            lex::lex_characters, token::Tokens, wu_manber_meyers::wu_manber_meyers_lcs_length,
        };

        #[test]
        fn test_fencepost_both_empty() {
            let mut tokens = Tokens::new();

            let z = lex_characters(&mut tokens, "");
            let length = wu_manber_meyers_lcs_length(&z, &z);
            assert_eq!(length, 0);
        }

        #[test]
        fn test_fencepost_one_empty() {
            let mut tokens = Tokens::new();

            let a = lex_characters(&mut tokens, "a");
            let z = lex_characters(&mut tokens, "");
            let length = wu_manber_meyers_lcs_length(&a, &z);
            assert_eq!(length, 0);
            let length = wu_manber_meyers_lcs_length(&z, &a);
            assert_eq!(length, 0);
        }

        #[test]
        fn test_fencepost_no_overlap() {
            let mut tokens = Tokens::new();

            let a = lex_characters(&mut tokens, "aa");
            let z = lex_characters(&mut tokens, "bbb");
            let length = wu_manber_meyers_lcs_length(&a, &z);
            assert_eq!(length, 0);
            let length = wu_manber_meyers_lcs_length(&z, &a);
            assert_eq!(length, 0);
        }

        #[test]
        fn test_example() {
            let mut tokens = Tokens::new();

            let abc = lex_characters(&mut tokens, "abc");
            let bac = lex_characters(&mut tokens, "bac");
            let length = wu_manber_meyers_lcs_length(&abc, &bac);
            assert_eq!(length, 2);
        }

        #[test]
        fn test_example_longer() {
            let mut tokens = Tokens::new();
            let abacad = lex_characters(&mut tokens, "abacad");
            let bacada = lex_characters(&mut tokens, "bacada");
            let length = wu_manber_meyers_lcs_length(&abacad, &bacada);
            assert_eq!(length, 5);
        }

        #[test]
        fn test_full_diff() {
            let mut tokens = Tokens::new();
            let abcde = lex_characters(&mut tokens, "abcde");
            let fghijk = lex_characters(&mut tokens, "fghij");
            let length = wu_manber_meyers_lcs_length(&abcde, &fghijk);
            assert_eq!(length, 0);
        }

        #[test]
        fn test_crash_uneven() {
            let mut tokens = Tokens::new();
            let a = lex_characters(&mut tokens, "a");
            let bab = lex_characters(&mut tokens, "bab");

            let lcs = wu_manber_meyers_lcs_length(&a, &bab);
            assert_eq!(lcs, 1);
        }
    }
}
