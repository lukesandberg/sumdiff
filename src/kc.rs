use assume::assume;
use std::{
    mem::{transmute, MaybeUninit},
    vec,
};

use crate::{
    lcs_utils::{
        binary_search_range, counting_sort, exponential_search_range, remove_suffixes_and_prefixes,
        Trimmed,
    },
    token::{CommonRange, Token, Tokens},
};

struct MatchList {
    // The ranges are are indexes into `right_indices`
    matches: Vec<(usize, usize)>,
    // N.B. These are 1-indexed to simplify later comparisons and allow for a zero sentinel value
    // This is a sorted list but will be mutated as we go.
    right_indices: Vec<usize>,
    // The highest index in thresh that a corresponding index in right_indices could be assigned to
    // This allows us to accelerate matching when matchlists are visited multiple times.
    starting_k: Vec<usize>,
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
            starting_k: vec![0; right_indices.len()],
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
/// It describes an algorithm with a runtime of `O(R + NL + NLog(N))` where:
///   * `R` is the number of matching locations between left and right, i.e.  `R=|{(i,j) s.t. left[i] = right[j]}|`
///      * Note, in the worst case `R` is O(N^2), and these cases are not particularly
///        rare (consider that every blank line will match and if 10% of lines are blank, R is now >= 0.01*N^2).
///        This is particularly likely with small-fixed alphabets (consider a genome comparison with only 4 tokens),
///        we should expect a quadratic number of matches.
///   * `L` is the length of the longest common subsequence
///   * `N` is the length of the input sequences
///
/// The `R` term is due to how the matchlists are iterated, the `NL` term tracks assignments to threshold and the
/// `NLog(N)` term is related to building the matchlist.
///
/// The implementation here improves on that in a few ways.
/// * Use counting sort to build the matchlist (dropping `O(NlogN)` -> `O(N)` for the sort)
///   - This is workable due to our tokenization technique, the universe of tokens is trivially bound by the size of the inputs.
///   - This improvement is useful for mostly theoretical reasons, the other optimizations are more practical.
///
/// * We use `exponential` search to find both insertion points in `thresh` and to iterate the `matchlist`. Despite the
///   faster searches, this doesn't improve worst case performance (still benchmarks measured consistent 30%+ performance
///   improvements).
///
/// * We use a set of `starting_k` values to accelerate the search of `thresh` for a given matchlist. This helps
///   when we visit the same matchlist multiple times since it allows us to converge on possible update locations
///   faster.  This optimization is more impactful for larger inputs, but ultimately doesn't close the performance
///   gap with meyers (up to 20% on larger inputs).
///
/// * We restrict the range of matches we consider on each row to only those that could increase the LCS length
///   - This relies on an insight from Hyyro that if the LCS length is `t` then the only matches that could be on the
///     LCS are those on diagonals (m-t)..(n-t) of the implicit matrix L from the Wagner-Fischer recurrence.
///     Because we refine our estimate of the LLCS as we execute (whenever we extend `thresh`) this means we can
///     narrow the range of matches we consider.  This is both a practical and a theoretical optimization to KC
///     It doesn't affect the number of assignments to `thresh` but it does reduce the number of comparisons we make for
///     each row by `L` (the length of the LCS).  This reduces the `R` term to `R/L`.
///
/// * Use a pool of links to amortize allocation overheads
///   - This doesn't affect asymptotic complexity but does reduce the number of allocations and it is a practical
///     enhancement.
///   - Compress DMatch nodes so they can efficiently represent ranges which are common
///   - TODO(luke): Using a bump allocator would be superior
///
/// Thus the overall complexity is still at least O(NL + R/L). This is more or less informed by the `thresh` datastructure.
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
    assume!(unsafe: (n > 1 && m > 0) || (n > 0 && m > 1), "n = {}, m = {}", n, m);

    // Intuition for `thresh`: Each value `j` in thresh at index `i` is a dominant match between left and right on
    // a common subsequence of length `i`. We store the predecessor of that match in the `links` array, such that
    // each entry in `thresh` is the head of a 'chain' of matches through the matchlist.

    // Thresh stores indexes of right where we could advance the LCS.  To simplify some comparisons we
    // Initialize it with m+1 so all values are greater than any index of right, and we make it as
    // long as the longest possible LCS.  Furthermore we initialize `thresh[0]=0` and use 1 based indexing
    // into right to resolve some fencepost issues.
    // drop Vec methods for extending and memory supporting it
    let mut thresh = vec![m + 1; std::cmp::min(n, m) + 1].into_boxed_slice();
    thresh[0] = 0;
    let mut matchlist = MatchList::build(tokens, left, right);
    // The worst case is that we need L^2 links, our link compression techniques will reduce this in practice
    // but not the theoretical worst case.
    let mut link_pool: Vec<DMatch> = Vec::with_capacity(n);
    // The zero entry is a dummy value
    link_pool.push(DMatch {
        i: n,
        j: m,
        prev: 0,
        len: 0,
    });
    let mut links: Vec<usize> = vec![0; std::cmp::min(n, m) + 1];
    let mut max_thresh = 1;

    // From "Bit Parallel LCS-length Computation Revisited"
    // https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=7b1385ba60875b219ce76d5dc0fb343f664c6d6a
    // We see the observation that if LLCS(left, right) >=t, then the only matches between left and right
    // that are on the LCS must be in diagonals -(n-t)..(m-t) of the implicit matrix L from the Wagner+Fishcer recurrence
    //
    // The implication of this is that we can restrict the range of matches we consider in each iteration to only those within
    // the range of diagonals defined by our lower bound estimate `max_thresh-1`. As we expand `thresh` our LCS estimate is
    // refined and we can further restrict the range of matches we consider.
    // This approach was also described by Bergroth in
    // "Utilizing Dynamically Updated Estimates in Solving the Longest Common Subsequence Problem", though I only discovered that
    // later.
    // TODO: implement the BestNext heuristic from Bergroth to further improve this.

    for (i, (mut mi, mut matches_end)) in matchlist.matches.iter().enumerate() {
        if mi == matches_end {
            continue;
        }
        // TODO: this is very conditional, can we simplify it?
        let lcs_estimate = max_thresh - 1;

        let max_diag = n - lcs_estimate;
        if i > max_diag {
            // we need to +1 because the values stored in matches are 1-indexed
            mi = binary_search_range(&matchlist.right_indices, mi, matches_end, i - max_diag + 1);
            if mi == matches_end {
                continue;
            }
        }

        let negative_min_diag = m - lcs_estimate;
        let max_j = i + negative_min_diag;
        if max_j < m {
            // +2 because we are 1 indexed and we want to be inclusive
            matches_end = binary_search_range(&matchlist.right_indices, mi, matches_end, max_j + 2);
            if mi == matches_end {
                continue;
            }
        }
        let mut k = matchlist.starting_k[mi];
        // These two fields serve as a tiny buffer to delay writes to the links array
        let mut r = 0;
        let mut c = links[0];
        let mut prev_thresh = thresh[k];
        loop {
            // Search for the next match index, greater than this one and smaller than the previous threshold
            // This allows us to exclude values that couldn't possibly decrease thresh, since all values in thresh
            // at indexes greater than k are strictly greater than prev_thresh
            mi = exponential_search_range(
                &matchlist.right_indices,
                mi,
                matches_end,
                prev_thresh + 1,
            );
            if mi >= matches_end {
                break;
            }
            let j = matchlist.right_indices[mi];
            // For a given mi value, see if we have already computed a higher starting location
            // because we iterate over the matchlist multiple times we can accelerate the search
            k = std::cmp::max(k, matchlist.starting_k[mi]);
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
            prev_thresh = thresh[new_k];
            k = new_k;
            // Only add a match if we are strictly decreasing the threshold at this index
            if j < prev_thresh {
                // TODO: when replacing a thresh that is <m+1 we could potentially shift
                // the previous value _up_ in thresh.  This is only valid if that value
                // also matches an index in left that is > i. To do that we would need
                // to store more information in thresh and possibly store matchlists in
                // both directions.
                // If when storing in `thresh` we also stored all the indices in `left` of
                // that value > i, then we could remove that value from the matchlist.
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
            }
            matchlist.starting_k[mi] = k; // Update the starting_k value for the next time we visit this matchlist
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

    #[test]
    fn test_kc_buggy_diagonal_range() {
        let mut tokens = Tokens::new();

        let left = lex_characters(&mut tokens, "b");
        let right = lex_characters(&mut tokens, "aaba");
        let kc = CommonRange::flatten(&kc_lcs(&tokens, &left, &right));
        check_is_lcs(&kc, &left, &right).unwrap();
    }

    #[test]
    fn test_kc_buggy_diagonal_range_2() {
        let mut tokens = Tokens::new();
        let left = lex_characters(&mut tokens, "abc");
        let right = lex_characters(&mut tokens, "cab");
        let kc = CommonRange::flatten(&kc_lcs(&tokens, &left, &right));
        check_is_lcs(&kc, &left, &right).unwrap();
    }
}
