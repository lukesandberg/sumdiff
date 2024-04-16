use crate::{
    lcs_utils::{common_prefix, common_suffix, remove_suffixes_and_prefixes, Trimmed},
    token::{CommonRange, Token},
};
use assume::assume;
use std::mem::MaybeUninit;

// TODO: in addition to returning the split point we could also hint at what direction the diagonal would be in
// The common_prefix/common_suffix scans in do_meyers_lcs are at least 50% wasted effort (one or both is guaranteed
// to be empty).

fn middle_snake(
    left: &[Token],
    right: &[Token],
    v_forward: &mut [usize],
    v_backward: &mut [usize],
) -> Option<(usize, usize)> {
    let n = left.len();
    let m = right.len();
    // If either is empty there is no snake to be found.
    // if either has length one and the other has length 2 or 1 then we are done since
    // prefix/suffixes have been removed so the comparisons are done.
    if m == 0 || n == 0 || (n == 1 && m <= 2) || (n <= 2 && m == 1) {
        return None;
    }
    if m > n {
        return do_middle_snake(m, right, n, left, v_backward, v_forward).map(|(x, y)| (y, x));
    }
    return do_middle_snake(n, left, m, right, v_forward, v_backward);
}
fn do_middle_snake(
    n: usize,
    left: &[Token],
    m: usize,
    right: &[Token],
    v_forward: &mut [usize],
    v_backward: &mut [usize],
) -> Option<(usize, usize)> {
    assume!(unsafe: n>=m);
    assume!(unsafe: n>0);
    assume!(unsafe: m>0);
    assume!(unsafe: n == left.len());
    assume!(unsafe: m == right.len());
    let delta = n - m;
    let delta_odd = delta % 2 == 1;
    let num_diagonals = n + m + 1;
    assume!(unsafe: v_forward.len()>=num_diagonals+2);
    assume!(unsafe: v_backward.len()>=num_diagonals+2);

    // Our iteration in the forward direction starts at diagonal 0
    // In the backward direction it starts at delta
    //
    //  0    1
    //   +---+---+
    //   |   |   |
    //-1 +---+---+
    //   |   |   |
    //-2 +---+---+
    //   |   |   |
    //   +---+---+
    //
    // The grid above assume `n=2` and `m=3` so `delta=-1`, The diagonals are numbered
    // by their upper left corner.
    //
    // The forward iteration starts from the upper left (diagonal `0`) and the backward
    // iteration starts from the lower right (diagonal `-1`, aka delta).

    // Our diagonals are (logically) numbered from `-m..n` so to shift them into positive numbers our offset needs to be `m`+1 for a boundary value.
    let v_offset = m + 1;
    // We store 1 based indexes of left in v_forward to allow for a 0 sentinel value.
    let mut f_lower = v_offset;
    let mut f_upper = v_offset;
    v_forward[f_upper] = 1;
    let mut b_lower = v_offset + delta;
    let mut b_upper = v_offset + delta;
    v_backward[b_upper] = n;
    let min_diagonal = v_offset - m; //  could be simplified to `1`
    let max_diagonal = v_offset + n; // aka `num_diagonals

    loop {
        if f_lower > min_diagonal {
            f_lower -= 1;
            v_forward[f_lower - 1] = 0;
        } else {
            f_lower += 1;
        }
        if f_upper < max_diagonal {
            f_upper += 1;
            v_forward[f_upper + 1] = 0;
        } else {
            f_upper -= 1;
        }
        debug_assert!((f_upper - f_lower) % 2 == 0);
        // NOTE: we can choose to iterate in increasing _or_ decreasing K, which will change the reported LCS.
        // Git apparently iterates in decreasing order.
        for k in (f_lower..=f_upper).step_by(2) {
            let l = v_forward[k - 1];
            let u = v_forward[k + 1];
            let mut x = if l < u {
                // Implicitly we are keeping the prior x value and incrementing 'y'
                u - 1
            } else {
                // Extend from the prior x and increment it.
                l
            };
            let mut y = v_offset + x - k;
            // Follow the diagonal as far as we can.
            // Using suffix arrays and LCP arrays, we can skip over common prefixes in O(1) time, but construction is
            // complex to say the least.

            assume!(unsafe: x <= n);
            // NOTE: y might be greater than m, but we will never access right[y] in that case.
            while x < n && y < m && left[x] == right[y] {
                x += 1;
                y += 1;
            }
            // TODO: In Git there is a heuristic to to exit the outer loop early whenever we find a long snake here
            // (where long is >20), this accelerates the algorithm when the diff is large since it allows us to
            // divide the space sooner.
            // TODO: consider an alternative heuristic to simply choose this snake when it crosses the middle of the
            // the space.  This would save an extra scan of the same snake in reverse.
            v_forward[k] = x + 1;
            if delta_odd && k >= b_lower && k <= b_upper && x >= v_backward[k] {
                return Some((x, y));
            }
        }

        if b_lower > min_diagonal {
            b_lower -= 1;
            v_backward[b_lower - 1] = usize::MAX;
        } else {
            b_lower += 1;
        }
        if b_upper < max_diagonal {
            b_upper += 1;
            v_backward[b_upper + 1] = usize::MAX;
        } else {
            b_upper -= 1;
        }
        debug_assert!((b_upper - b_lower) % 2 == 0);
        // Now follow the reverse diagonals
        for k in (b_lower..=b_upper).step_by(2) {
            let l = v_backward[k - 1];
            let u = v_backward[k + 1];
            let mut x = if l < u {
                l // keep x and decrement y
            } else {
                u - 1 // decrement x
            };

            // The corresponding y value is x - k but our k is offset by v_offset so correct for that.
            // NOTE: it is possible for this operation to underflow since the 'decrement' in the first
            // branch above might simply be out of range, in that case we just keep the `x` value but
            // we need to not under flow here.
            let mut y = (v_offset + x) as isize - k as isize;

            // Follow the diagonal as far as we can.
            assume!(unsafe: x <= n && y <= m as isize);
            while x > 0 && y > 0 && left[x - 1] == right[y as usize - 1] {
                x -= 1;
                y -= 1;
            }

            v_backward[k] = x;
            if !delta_odd && k >= f_lower && k <= f_upper && x < v_forward[k] {
                return Some((x, y as usize));
            }
        }

        // TODO: in git there is a further heuristic to exit early when the edit cost (`d`) gets too large.
        // When that happens we just recover the longest snake from the current round.  This allows us to avoid
        // full quadratic behavior in the worst case of large diffs, at the expensive of choosing non-optimal split points.
        // Setting this location to be `max_distance`/2 should preserve the overall complexity of the algorithm
    }
}

// Returns the longest common subsequence of `right` and `left` using the Meyers LCS algorithm.
//
// This uses this linear space variation that executes in `O((N+M)D)` time and `O(N+M)` space.
// See http://www.xmailserver.org/diff2.pdf for the original paper
// The discussion at https://blog.jcoglan.com/2017/02/12/the-myers-diff-algorithm-part-1/ is also helpful for a more detailed
// perspective and https://blog.robertelder.org/diff-algorithm/ offers some insight as well.
// The xdiff implementation in side of git, uses some of the advice from the last link (independently discovered afaict)
//
// The implementation here is straightforward to the paper with some care taken to
// * deal with using unsigned indices as rust requires
// * reduce memory usage and initialization overhead
//
pub fn meyers_lcs(left: &[Token], right: &[Token]) -> Vec<CommonRange> {
    let Trimmed {
        left,
        right,
        mut matches,
        prefix,
        suffix,
    } = match remove_suffixes_and_prefixes(left, right) {
        Ok(matches) => {
            return matches;
        }
        Err(t) => t,
    };
    let n = left.len();
    let m = right.len();
    // In `do_middle_snake` we need to traverse all the diagonals from `-m..=n`
    // and we also place sentinel values on either side to simplify loop conditions.
    // so `+2` for the sentinels and `+1` because the range of diagonals is inclusive.
    let num_diagonals_and_slack = n + m + 3;
    // This vector is used to store the forward and backward diagonals as we go
    // SAFETY: we will only access values that have been initialized on each call
    // to `middle_snake`.  To reduce allocations we construct one vector and
    // create 2 slices.
    // A better approach might be to use a dynamic array, if the LCS is long we won't end up
    // needinng the full array in practice, however one allocation of uninitialized memory should
    // be very cheap in practice.
    let mut v: Vec<usize> = unsafe {
        std::mem::transmute(vec![
            MaybeUninit::<usize>::uninit();
            2 * num_diagonals_and_slack
        ])
    };
    let (forward, backward) = v.split_at_mut(num_diagonals_and_slack);
    do_meyers_lcs(prefix, left, prefix, right, forward, backward, &mut matches);
    if suffix > 0 {
        matches.push(CommonRange {
            left_start: prefix + left.len(),
            right_start: prefix + right.len(),
            length: suffix,
        });
    }

    matches
}

fn do_meyers_lcs(
    mut left_offset: usize,
    mut left: &[Token],
    mut right_offset: usize,
    mut right: &[Token],
    v_forward: &mut [usize],
    v_backward: &mut [usize],
    out: &mut Vec<CommonRange>,
) {
    // TODO: If the length is small enough we can use a simpler quadratic algorithm which would short circuit
    // the analysis here.  One of the bit-parallel algorithms would be a good choice.
    let prefix = common_prefix(left, right);
    if prefix > 0 {
        out.push(CommonRange {
            left_start: left_offset,
            right_start: right_offset,
            length: prefix,
        });
        left_offset += prefix;
        right_offset += prefix;
        left = &left[prefix..];
        right = &right[prefix..];
    }
    let suffix = common_suffix(left, right);
    let suffix_range: Option<CommonRange>;
    if suffix > 0 {
        suffix_range = Some(CommonRange {
            left_start: left_offset + left.len() - suffix,
            right_start: right_offset + right.len() - suffix,
            length: suffix,
        });
        left = &left[..left.len() - suffix];
        right = &right[..right.len() - suffix];
    } else {
        suffix_range = None;
    }

    if let Some((x, y)) = middle_snake(left, right, v_forward, v_backward) {
        do_meyers_lcs(
            left_offset,
            &left[..x],
            right_offset,
            &right[..y],
            v_forward,
            v_backward,
            out,
        );
        let actual_left_start = left_offset + x;
        let actual_right_start = right_offset + y;

        do_meyers_lcs(
            actual_left_start,
            &left[x..],
            actual_right_start,
            &right[y..],
            v_forward,
            v_backward,
            out,
        );
    }
    if let Some(r) = suffix_range {
        out.push(r);
    }
}

/// An implementation of the Myers LCS algorithm for calculating the length of an LCS
/// - http://www.xmailserver.org/diff2.pdf is the original paper that is quite readable.
/// - https://blog.jcoglan.com/2017/02/12/the-myers-diff-algorithm-part-1/ was a secondary resource
///   with more examples.
///
/// This implementation is very true to the paper. There are some slight differences to avoid dealing
/// with negative numbers as array indices.  The algorithm is O((N+M)D) where D is the length of the
/// diff.
pub fn meyers_lcs_length(mut left: &[Token], mut right: &[Token]) -> usize {
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
    let max_distance = n + m;
    let num_diagonals = n + m + 1;
    // Save time by not initializing the entire vector
    // SAFETY: we will only access values that have been initialized.
    // For the first iteration we only access `k+1` which is pre-initialized to 0.
    // Prior to each iteration of the inner loop we assign values outside the [lower,upper] range.
    // Within the loop we only access values offset by 1 from the current value of k and k increments by 2
    // so these are either the preinitialized values, or the values assigned by the prior iteration (when `d` had
    // a different parity).

    // In this way we only ever access values of v that were initialized by a prior iteration of the loop.
    let mut v: Vec<usize> =
        unsafe { std::mem::transmute(vec![MaybeUninit::<usize>::uninit(); num_diagonals + 2]) };

    // We store 1 based `x` indexes in `v` to avoid negative numbers and reserve 0 as a lowermost value.
    let v_offset = m + 1;
    let mut lower = v_offset;
    let mut upper = v_offset;
    let mut min_diagonal = v_offset - m; // aka `1`
    let mut max_diagonal = v_offset + n; // aka `num_diagonals`

    // Because we trimmed the prefix and suffix we can start at 1
    // but we have to bootstrap the first value of the array.
    v[v_offset] = 1;
    let mut d = 1;
    // This loop is guaranteed to only execute `max_distance` iterations, but the compiler cannot see that
    // so use an infinite loop construct.
    loop {
        // Adjust lower and upper for the next diagonal
        // avoid exceeding the bounds of the array to reduce work by not traversing 'fake' locations
        // Also, put a sentinel value at the ends to avoid some boundary conditions in the inner loop.
        if lower > min_diagonal {
            lower -= 1;
            v[lower - 1] = 0;
        } else {
            // this means we have found the bottom edge, so we can start moving up diagonals
            lower += 1;
            min_diagonal += 1;
            // NOTE: this optimization only helps when `d>m` `d>n` at which point performance is
            // already disaterous, however this is useful to allow us to decrease the
            // the space of `v` that is accessed
        }
        if upper < max_diagonal {
            upper += 1;
            v[upper + 1] = 0;
        } else {
            // this means we have found the right hand side, so we can start decreasing the max
            // diagonal.
            upper -= 1;
            max_diagonal -= 1;
        }
        assume!(unsafe: (upper - lower) % 2 == 0);
        // k is our diagonal index
        for k in (lower..=upper).step_by(2) {
            let l = v[k - 1];
            let u = v[k + 1];
            let mut x = if l < u {
                // Implicitly we are keeping the prior x value and incrementing 'y'
                // NOTE: this might put y out of range
                u - 1
            } else {
                // Extend from the prior x and increment it.
                // because `v` is one indexed, it is already incremented
                l
            };

            // The corresponding y value is x - k but our k is offset by offset so correct for that.
            // The order of operations here is important to avoid underflow for our unsigned variables.
            let mut y = v_offset + x - k;
            // Follow the diagonal as far as we can.
            // Using suffix arrays and LCP arrays, we can skip over common prefixes in O(1) time, but construction is
            // complex to say the least.
            while x < n && y < m && left[x] == right[y] {
                x += 1;
                y += 1;
            }
            v[k] = x + 1;
            if x == n && y == m {
                // d is the number of edits, so the length of the LCS is (n+m-d)/2
                // This is because every edit accounts for one character on either side, but every common character
                // accounts for one character on both sides.
                return trimmed_length + (max_distance - d) / 2;
            }
        }
        d += 1;
    }
}

#[cfg(test)]
mod tests {

    mod length {
        use crate::{lex::lex_characters, meyers::meyers_lcs_length, token::Tokens};

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

        #[test]
        fn test_full_diff() {
            let mut tokens = Tokens::new();
            let abcde = lex_characters(&mut tokens, "abcde");
            let fghijk = lex_characters(&mut tokens, "fghijk");
            let length = meyers_lcs_length(&abcde, &fghijk);
            assert_eq!(length, 0);
        }
    }

    mod lcs {
        use crate::{
            lcs_utils::check_is_lcs,
            lex::lex_characters,
            meyers::meyers_lcs,
            token::{CommonRange, Tokens},
        };

        #[test]
        fn test_fencepost_both_empty() {
            let mut tokens = Tokens::new();

            let z = lex_characters(&mut tokens, "");
            let lcs = meyers_lcs(&z, &z);
            assert_eq!(lcs.len(), 0);
        }

        #[test]
        fn test_fencepost_one_empty() {
            let mut tokens = Tokens::new();

            let a = lex_characters(&mut tokens, "a");
            let z = lex_characters(&mut tokens, "");
            let lcs = meyers_lcs(&a, &z);
            check_is_lcs(&CommonRange::flatten(&lcs), &a, &z).unwrap();
            let lcs = meyers_lcs(&z, &a);
            check_is_lcs(&CommonRange::flatten(&lcs), &z, &a).unwrap();
        }

        #[test]
        fn test_fencepost_no_overlap() {
            let mut tokens = Tokens::new();

            let a = lex_characters(&mut tokens, "a");
            let bb = lex_characters(&mut tokens, "bb");
            let lcs = meyers_lcs(&a, &bb);
            check_is_lcs(&CommonRange::flatten(&lcs), &a, &bb).unwrap();
            let lcs = meyers_lcs(&bb, &a);
            check_is_lcs(&CommonRange::flatten(&lcs), &bb, &a).unwrap();
        }

        #[test]
        fn test_example() {
            let mut tokens = Tokens::new();

            let abc = lex_characters(&mut tokens, "abc");
            let bac = lex_characters(&mut tokens, "bac");
            let lcs = meyers_lcs(&abc, &bac);
            check_is_lcs(&CommonRange::flatten(&lcs), &abc, &bac).unwrap();
        }

        #[test]
        fn test_example_longer() {
            let mut tokens = Tokens::new();
            let abacad = lex_characters(&mut tokens, "abacad");
            let bacada = lex_characters(&mut tokens, "bacada");
            let lcs = meyers_lcs(&abacad, &bacada);
            check_is_lcs(&CommonRange::flatten(&lcs), &abacad, &bacada).unwrap();
        }

        #[test]
        fn test_full_diff_odd_delta() {
            let mut tokens = Tokens::new();
            let abcde = lex_characters(&mut tokens, "abcde");
            let fghijk = lex_characters(&mut tokens, "fghijk");
            let lcs = meyers_lcs(&abcde, &fghijk);
            check_is_lcs(&CommonRange::flatten(&lcs), &abcde, &fghijk).unwrap();
        }

        #[test]
        fn test_full_diff_even_delta() {
            let mut tokens = Tokens::new();
            let abcde = lex_characters(&mut tokens, "abcde");
            let fghij = lex_characters(&mut tokens, "fghij");
            let lcs = meyers_lcs(&abcde, &fghij);
            check_is_lcs(&CommonRange::flatten(&lcs), &abcde, &fghij).unwrap();
        }

        #[test]
        fn test_crash_uneven() {
            let mut tokens = Tokens::new();
            let abcde = lex_characters(&mut tokens, "a");
            let fghij = lex_characters(&mut tokens, "baaab");

            let lcs = meyers_lcs(&abcde, &fghij);
            check_is_lcs(&CommonRange::flatten(&lcs), &abcde, &fghij).unwrap();
        }
    }
}
