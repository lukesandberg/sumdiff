use std::{
    mem::{transmute, MaybeUninit},
    vec,
};

use crate::token::{CommonRange, Token, Tokens};

/// Peforms a simple 'counting' sort on the token array, but returns the indices of the sorted token locations
/// rather than sorting the tokens themselves.
/// So output[0] is the index of the first occurrence of the smallest token, and so on.
fn counting_sort(arr: &[Token], max_token: usize) -> Vec<usize> {
    let mut counts: Vec<u32> = vec![0; max_token];
    // Count the number of occurrences of each token
    for &i in arr.iter() {
        counts[i as usize] += 1;
    }
    // Tranform to a prefix sum
    let mut accum = counts[0];
    for v in counts.iter_mut().skip(1) {
        accum += *v;
        *v = accum;
    }
    // We don't need to initialize the slice since we will overwrite all the elements
    let mut indices = vec![MaybeUninit::<usize>::uninit(); arr.len()];
    for i in (0..arr.len()).rev() {
        let token = arr[i];
        let count = counts[token as usize] - 1;
        counts[token as usize] = count;
        indices[count as usize] = MaybeUninit::new(i);
    }
    // SAFETY: safe because we are guaranteed to have initialized all the elements of indices
    unsafe { transmute(indices) }
}

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

fn common_prefix<T: PartialEq + PartialOrd>(a: &[T], b: &[T]) -> usize {
    let mut i = 0;
    while i < a.len() && i < b.len() && a[i] == b[i] {
        i += 1;
    }
    i
}
fn common_suffix<T: PartialEq + PartialOrd>(a: &[T], b: &[T]) -> usize {
    let n = a.len();
    let m = b.len();
    let limit = std::cmp::min(n, m);
    let mut i = 0;
    while i < limit && a[n - i - 1] == b[m - i - 1] {
        i += 1;
    }
    i
}

/// A branchfree implementation of binary search over a subrange
/// It is the callers responsibility to ensure that offset and end are in range
///
/// For our usecase this significantly outperforms the standard library binary search.
fn binary_search_range(arr: &[usize], offset: usize, end: usize, x: usize) -> usize {
    let mut low = offset;
    let mut size = end - low;
    // The loop condition is easy to predict but the comparison within the
    // loop is not.  So we take care to ensure that that uses a conditional move
    // instruction.
    while size > 1 {
        let half = size / 2;
        let mid = low + size / 2;
        // SAFETY: mid is always >= offset and < end
        let mid_value = unsafe { *arr.get_unchecked(mid) };
        // This is compiled to conditional moves
        low = if mid_value < x { mid } else { low };
        size -= half;
    }
    let low_value = unsafe { *arr.get_unchecked(low) };
    // This branch is unpredictable but it is only taken once and generally compiled to a conditional increment.
    if low_value < x {
        low + 1
    } else {
        low
    }
}

/// Returns the index of the greatest element of `arr` that is <= `x` using the exponential search algorithm.
///
/// Exponential search is a good fit for our use case because it is very efficient when the target is expected to be near the beginning of the array.
fn exponential_search_range(arr: &[usize], offset: usize, end: usize, x: usize) -> usize {
    let mut bound = 1;
    let mut probe = offset + bound;
    while probe < end && unsafe { *arr.get_unchecked(probe) } < x {
        bound *= 2;
        probe = offset + bound;
    }
    binary_search_range(
        arr,
        // start at the previous search location which we know is <=x
        offset + bound / 2,
        // end needs to be at most the size of the array or one past the last test which is >= x
        std::cmp::min(end, probe + 1),
        x,
    )
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
    let mut n = left.len();
    let mut m = right.len();
    if n == 0 || m == 0 {
        return vec![];
    }
    let prefix = common_prefix(left, right);
    let mut matches = Vec::new();
    if prefix > 0 {
        matches.push(CommonRange {
            left_start: 0,
            right_start: 0,
            length: prefix,
        });
        // Exact match or one sequence is a prefix of the other
        if prefix == n || prefix == m {
            return matches;
        }
    }
    // wait to append the suffix matches until after the main loop
    let suffix = common_suffix(left, right);
    n -= suffix;
    m -= suffix;

    let left = &left[prefix..n];
    let right = &right[prefix..m];
    n -= prefix;
    m -= prefix;
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
        matches[1..].reverse();
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
    fn test_common_prefix() {
        let a = vec![1, 2, 3, 4, 5];
        let b = vec![1, 2, 3, 4, 5];
        assert_eq!(common_prefix(&a, &b), 5);
        let a = vec![1, 2, 3, 4, 5];
        let b = vec![1, 2, 3, 4, 6];
        assert_eq!(common_prefix(&a, &b), 4);
        let a = vec![1, 2, 3, 4, 5];
        let b = vec![1, 2, 3, 4, 5, 6];
        assert_eq!(common_prefix(&a, &b), 5);
    }
    #[test]
    fn test_common_suffix() {
        let a = vec![1, 2, 3, 4, 5];
        let b = vec![1, 2, 3, 4, 5];
        assert_eq!(common_suffix(&a, &b), 5);
        let a = vec![1, 2, 4, 4, 5];
        let b = vec![1, 2, 3, 4, 5];
        assert_eq!(common_suffix(&a, &b), 2);
        let a = vec![1, 2, 3, 4, 5];
        let b = vec![1, 2, 3, 4, 5, 6];
        assert_eq!(common_suffix(&a, &b), 0);
    }
    #[test]
    fn test_counting_sort() {
        let arr = vec![1, 2, 1, 2, 1, 2];
        let indices = counting_sort(&arr, 10);
        assert_eq!(indices, vec![0, 2, 4, 1, 3, 5]);
    }
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

    #[test]
    fn test_binary_search_suffix() {
        let arr = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        assert_eq!(binary_search_range(&arr, 0, arr.len(), 0), 0);
        assert_eq!(binary_search_range(&arr, 0, arr.len(), 5), 5);
        assert_eq!(binary_search_range(&arr, 0, arr.len(), 9), 9);
        assert_eq!(binary_search_range(&arr, 0, arr.len(), 10), 10);
    }

    #[test]
    fn test_exponential_search_suffix() {
        let arr = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        assert_eq!(exponential_search_range(&arr, 0, arr.len(), 0), 0);
        assert_eq!(exponential_search_range(&arr, 0, arr.len(), 5), 5);
        assert_eq!(exponential_search_range(&arr, 0, arr.len(), 9), 9);
        assert_eq!(exponential_search_range(&arr, 0, arr.len(), 10), 10);
    }
}
