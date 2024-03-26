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
    // At this point we know that the inputs are non-empty, at least one is longer than 1, and they have no common
    // prefix or suffix.
    debug_assert!((n > 1 && m > 0) || (n > 0 && m > 1), "n = {}, m = {}", n, m);

    // Thresh stores indexes of right where we could advance the LCS.  To simplify some comparisons we
    // Initialize it with m+1 so all values are greater than any index of right, and we make it as
    // long as the longest possible LCS.  Furthermore we initialize `thresh[0]=0` and use 1 based indexing
    // into right to resolve some fencepost issues.
    // drop Vec methods for extending and memory supporting it
    let mut thresh = vec![m + 1; std::cmp::min(n, m) + 1].into_boxed_slice();
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
    for (i, range) in matchlist.matches.iter().enumerate() {
        if range.0 == range.1 {
            continue;
        }
        let matches = &matchlist.right_indices[range.0..range.1];
        let mut k = 0;
        // These two fields serve as a tiny buffer to delay writes to the links array
        let mut r = 0;
        let mut c = links[0];
        let mut mi = 0;
        // This also reserves the 0th entry for the dummy value above and simplifies read-out.
        let mut j = matches[mi];
        loop {
            let new_k = exponential_search_range(&thresh, k, max_thresh, j);
            debug_assert!(
                thresh[new_k - 1] < j && thresh[new_k] >= j,
                "new_k={}, j={}, thresh[{}..{}]={:?}",
                new_k,
                j,
                k,
                max_thresh + 1,
                &thresh[k..max_thresh + 1]
            );
            let prev_thresh = thresh[new_k];
            k = new_k;
            // Only add a match if we are strictly decreasing the threshold
            if j < prev_thresh {
                thresh[new_k] = j;
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
            }
            // Search for the next match index, greater than this one and smaller than the previous threshold
            // This allows us to exclude values that couldn't possibly decrease thresh, since all values in thresh
            // at indexes greater than k are strictly greater than prev_thresh
            mi += 1;
            if mi >= matches.len() {
                break;
            }
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
    use crate::{lcs_utils::check_is_lcs, lex::lex_characters};

    use super::*;

    #[test]
    fn test_build_matchlist() {
        let mut tokens = Tokens::new();
        let abc = lex_characters(&mut tokens, "abc");
        let bccd = lex_characters(&mut tokens, "bccd");
        let matchlist = MatchList::build(&tokens, &abc, &bccd);
        assert_eq!(matchlist.matches, vec![(0, 0), (0, 1), (1, 3)]);
    }

    #[test]
    fn test_kc_lcs() {
        let mut tokens = Tokens::new();

        let abc = lex_characters(&mut tokens, "abc");
        let bac = lex_characters(&mut tokens, "bac");
        let lcs = CommonRange::flatten(&kc_lcs(&tokens, &abc, &bac));
        assert_eq!(lcs, vec![(1, 0), (2, 2)]);

        let abacad = lex_characters(&mut tokens, "abacad");
        let bacada = lex_characters(&mut tokens, "bacada");
        let lcs = CommonRange::flatten(&kc_lcs(&tokens, &abacad, &bacada));
        assert_eq!(lcs, vec![(1, 0), (2, 1), (3, 2), (4, 3), (5, 4)]);
    }

    #[test]
    fn test_kc_not_subsequence() {
        let mut tokens = Tokens::new();
        let left = ["22", "15", "19", "15", "25", "25"]
            .iter()
            .map(|s| tokens.get_token(s.as_bytes().to_vec()))
            .collect::<Vec<Token>>();
        let right = ["26", "15", "30", "15", "25", "28"]
            .iter()
            .map(|s| tokens.get_token(s.as_bytes().to_vec()))
            .collect::<Vec<Token>>();
        let kc = CommonRange::flatten(&kc_lcs(&tokens, &left, &right));
        check_is_lcs(&kc, &left, &right).unwrap();
    }

    /// This test triggered a fencepost error in the prefix/suffix logic.
    #[test]
    fn test_kc_buggy_suffix_logic() {
        let mut tokens = Tokens::new();

        let left = lex_characters(&mut tokens, "aba");
        let right = lex_characters(&mut tokens, "b");
        let kc = CommonRange::flatten(&kc_lcs(&tokens, &left, &right));
        check_is_lcs(&kc, &left, &right).unwrap();
    }

    #[test]
    fn test_kc_buggy_thresh_update() {
        let mut tokens = Tokens::new();

        let left = lex_characters(&mut tokens, "b");
        let right = lex_characters(&mut tokens, "aba");
        let kc = CommonRange::flatten(&kc_lcs(&tokens, &left, &right));
        check_is_lcs(&kc, &left, &right).unwrap();
    }
}
