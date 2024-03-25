use std::{
    mem::{transmute, MaybeUninit},
    vec,
};

use crate::{
    lcs_utils::{counting_sort, exponential_search_range, remove_suffixes_and_prefixes, Trimmed},
    token::{CommonRange, Token, Tokens},
};

struct MatchList {
    // The ranges are are indexes into `right_indices`
    matches: Vec<(usize, usize)>,
    // N.B. These are 1-indexed to simplify later comparisons and allow for a zero sentinel value
    // This is a sorted list but will be mutated as we go.
    right_indices: Vec<usize>,
}

impl MatchList {
    // Builds a matchlist from two token lists in O(N+M) time with O(N+M) space
    fn build(tokens: &Tokens, left: &[Token], right: &[Token]) -> MatchList {
        debug_assert!(!left.is_empty());
        debug_assert!(!right.is_empty());
        let max_token = tokens.token_upper_bound();
        let left_indices = counting_sort(left, max_token);
        let mut right_indices = counting_sort(right, max_token);
        let mut matches = vec![MaybeUninit::<(usize, usize)>::uninit(); left.len()];
        let mut i = 0;
        let mut ri = 0;
        while i < left_indices.len() {
            let mut li = left_indices[i];
            let token = left[li];
            while ri < right_indices.len() && right[right_indices[ri]] < token {
                ri += 1;
            }
            let match_start = ri;
            while ri < right_indices.len() && right[right_indices[ri]] == token {
                ri += 1;
            }
            loop {
                matches[li] = MaybeUninit::new((match_start, ri));
                i += 1;
                if i >= left_indices.len() {
                    break;
                }
                li = left_indices[i];
                if left[li] != token {
                    break;
                }
            }
        }
        for item in &mut right_indices {
            *item += 1
        }
        MatchList {
            matches: unsafe { transmute(matches) },
            right_indices,
        }
    }
}

struct DMatch {
    i: usize,
    j: usize,
    prev: usize,
    len: usize,
}

/// Returns the longest common subsequence of two sequences of tokens using
/// the algorithm described by [Kuo and Cross](https://dl.acm.org/doi/pdf/10.1145/74697.74702)
/// which improves on the classic algorithm by Hunt and Szymanski.
///
/// This also includes a few minor optimizations over the algorithms described in the paper
/// * Use exponential search to find insertion points in thresh
///   - Saves substantial time in the worst case
/// * Use exponential search to find offsets in matchlist for each row
///   - for larger matchlists this is substantial
/// * Use counting sort to build the matchlist (dropping O(NlogN) -> O(N) for the sort)
///   - This is workable due to our tokenization technique, the universe of tokens is trivially bound by the size of the inputs.
/// * Use a pool of links to amortize allocation overheads
///   - This doesn't affect asymptotic complexity but does reduce the number of allocations and it is a practical
///     enhancement.
///   - Compress DMatch nodes so they can efficiently represent ranges.
///
/// The runtime of this algorithm is complex
///  - The matchlist is built in O(N) time (due to counting sort !)
///  - The readout at the end takes at O(N) time (generally much faster)
///  - Anaylyzing the loop is difficult
///     - For each match we run 2 exponential searches, one is on thresh (thus O(logN)) the other is on the matchlist which is worse case
///       O(logR) (where R is the size of the matchlist), but averaging across all incides of `n` we should expect O(log(sqrt(R)) instead
///       In the worst case R is O(N^2), so again we should expect each match to take O(logN) time.
///     - The number of matches is not easy to bound however. `R` is a trivial upper bound but it is not tight.
///
/// Thus the overall complexity is O(NlogN).
///
/// * TODO(luke): find a way to use the implicit histogram of the matchlist to leverage a patience sort technique
pub fn kc_lcs(tokens: &Tokens, left: &[Token], right: &[Token]) -> Vec<CommonRange> {
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
    debug_assert!((n > 1 && m > 0) || (n > 0 && m > 1), "n = {}, m = {}", n, m);

    // At this point we know that the inputs are non-empty, at least one is longer than 1, and they have no common
    // prefix or suffix.

    // drop Vec methods for extending and memory supporting it
    let mut thresh = vec![n; m + 1];
    thresh[0] = 0;
    let matchlist = MatchList::build(tokens, left, right);
    // TODO: justify capacities
    let mut link_pool: Vec<DMatch> = Vec::with_capacity(n);
    // The zero entry is a dummy value
    link_pool.push(DMatch {
        i: n,
        j: m,
        prev: 0,
        len: 0,
    });
    let mut links: Vec<usize> = vec![0; m + 1];
    let mut max_thresh = 1;
    for i in 0..n {
        let range = &matchlist.matches[i];
        let matches = &matchlist.right_indices[range.0..range.1];
        if matches.is_empty() {
            continue;
        }
        let mut k = 0;
        // These two fields serve as a tiny buffer to delay writes to the links array
        let mut r = 0;
        let mut c = links[0];
        let mut mi = 0;
        // This also reserves the 0th entry for the dummy value above and simplifies read-out.
        let mut j = matches[mi];
        loop {
            k = exponential_search_range(&thresh, k, max_thresh, j);
            let prev_thresh = thresh[k];
            thresh[k] = j;
            max_thresh = std::cmp::max(max_thresh, k + 1);
            let prev = links[k - 1];
            links[r] = c;
            r = k;
            let prev_link = link_pool.get_mut(prev).unwrap();
            // Check if we are simply extending a previous match.
            if prev_link.i + prev_link.len == i && prev_link.j + prev_link.len == j - 1 {
                prev_link.len += 1;
                c = prev;
            } else {
                c = link_pool.len();
                // This is the most expensive thing in the loop, extending the link pool is costly.
                // What we want to know is when it is safe to reuse a prior link.
                link_pool.push(DMatch {
                    i,
                    j: j - 1,
                    prev,
                    len: 1,
                });
            }
            // Search for the next match index, greater than this one and smaller than the previous threshold
            // This allows us to exclude values that couldn't possibly decrease thresh, since all values in thresh
            // at indexes greater than k are strictly greater than prev_thresh
            mi += 1;
            mi = exponential_search_range(matches, mi, matches.len(), prev_thresh + 1);
            if mi >= matches.len() {
                break;
            }
            j = matches[mi]
        }
        links[r] = c;
    }

    if max_thresh > 1 {
        let mut i = links[max_thresh - 1];
        let link = &link_pool[i];
        let mut c: CommonRange = CommonRange {
            left_start: prefix + link.i + link.len,
            right_start: prefix + link.j + link.len,
            length: 0,
        };
        loop {
            let link = &link_pool[i];
            if c.left_start == prefix + link.i + link.len
                && c.right_start == prefix + link.j + link.len
            {
                c.length += link.len;
                c.left_start -= link.len;
                c.right_start -= link.len;
            } else {
                debug_assert!(c.length > 0);
                matches.push(c);

                c = CommonRange {
                    left_start: prefix + link.i,
                    right_start: prefix + link.j,
                    length: link.len,
                }
            };

            if link.prev == 0 {
                debug_assert!(c.length > 0);
                matches.push(c);
                break;
            }
            i = link.prev;
        }
        if prefix > 0 {
            matches[1..].reverse();
        } else {
            matches.reverse();
        }
    }
    if suffix > 0 {
        matches.push(CommonRange {
            left_start: prefix + n,
            right_start: prefix + m,
            length: suffix,
        });
    }
    matches
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_matchlist() {
        let mut tokens = Tokens::new();
        let a = tokens.get_token(b"a".to_vec());
        let b = tokens.get_token(b"b".to_vec());
        let c = tokens.get_token(b"c".to_vec());
        let d = tokens.get_token(b"d".to_vec());
        let matchlist = MatchList::build(&tokens, &vec![a, b, c], &vec![b, c, c, d]);
        assert_eq!(matchlist.matches, vec![(0, 0), (0, 1), (1, 3)]);
    }
}
