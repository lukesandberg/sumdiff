use std::{
    mem::{transmute, MaybeUninit},
    ops::Range,
    vec,
};

use crate::token::{CommonRange, Token, Tokens};

/// Peforms a simple 'counting' sort on the token array, but returns the indices of the sorted token locations
/// rather than sorting the tokens themselves.
/// So output[0] is the index of the first occurrence of the smallest token, and so on.
fn counting_sort(arr: &[Token], max_token: usize) -> Vec<usize> {
    let mut counts: Vec<u32> = vec![0; max_token as usize];
    // Count the number of occurrences of each token
    for &i in arr.iter() {
        counts[i as usize] += 1;
    }
    // Tranform to a prefix sum
    let mut accum = counts[0];
    for i in 1..counts.len() {
        accum += counts[i];
        counts[i] = accum;
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
    matches: Vec<Range<usize>>,
    // N.B. These are 1-indexed to simplify later comparisons and allow for a zero sentinel value
    // This is a sorted list but will be mutated as we go.
    right_indices: Vec<usize>,
}

impl MatchList {
    // Builds a matchlist from two token lists in O(N+M) time with O(N+M) space
    fn build(tokens: &Tokens, left: &[Token], right: &[Token]) -> MatchList {
        debug_assert!(left.len() > 0);
        debug_assert!(right.len() > 0);
        let max_token = tokens.token_upper_bound();
        let left_indices = counting_sort(left, max_token);
        let mut right_indices = counting_sort(right, max_token);
        let mut matches: Vec<Range<usize>> = vec![0..0; left.len()];
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
                matches[li] = match_start..ri;
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
        for i in 0..right_indices.len() {
            right_indices[i] += 1;
        }
        MatchList {
            matches,
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
    return binary_search_range(
        arr,
        // start at the previous search location which we know is <=x
        offset + bound / 2,
        // end needs to be at most the size of the array or one past the last test which is >= x
        std::cmp::min(end, probe + 1),
        x,
    );
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
///   - The biggest opportunity is finding ways to allocate fewer DMatch nodes.
///
/// * TODO(luke): find a way to use the implicit histogram of the matchlist to leverage a patience sort technique
pub fn kc_lcs(tokens: &Tokens, left: &Vec<Token>, right: &Vec<Token>) -> Vec<CommonRange> {
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
        if range.is_empty() {
            continue;
        }
        let matches = &matchlist.right_indices[range.clone()];
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
            mi = exponential_search_range(&matches, mi, matches.len(), prev_thresh + 1);
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
        assert_eq!(matchlist.matches, vec![0..0, 0..1, 1..3]);
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
    #[cfg(test)]
    mod kc_tests {
        use super::*;

        fn naive_lcs_length(a: &Vec<Token>, b: &Vec<Token>) -> usize {
            let n = a.len();
            let m = b.len();
            let mut table = vec![vec![0; m + 1]; n + 1];
            for i in 1..=n {
                for j in 1..=m {
                    if a[i - 1] == b[j - 1] {
                        table[i][j] = table[i - 1][j - 1] + 1;
                    } else {
                        table[i][j] = std::cmp::max(table[i - 1][j], table[i][j - 1]);
                    }
                }
            }
            table[a.len()][b.len()]
        }

        #[test]
        fn test_naive_lsc_length() {
            let mut tokens = Tokens::new();
            let a = tokens.get_token(b"a".to_vec());
            let b = tokens.get_token(b"b".to_vec());
            let c = tokens.get_token(b"c".to_vec());
            let d = tokens.get_token(b"d".to_vec());
            assert_eq!(naive_lcs_length(&vec![a, b, c], &vec![b, a, c]), 2);
            assert_eq!(
                naive_lcs_length(&vec![a, b, a, c, a, d], &vec![b, a, c, a, d, a]),
                5
            );
        }

        fn chars_to_toks(tokens: &mut Tokens, s: &str) -> Vec<Token> {
            let mut toks = Vec::new();
            s.chars().for_each(|c| {
                toks.push(tokens.get_token(vec![c as u8]));
            });
            toks
        }

        struct LcsTest {
            left: &'static str,
            right: &'static str,
            length: Option<usize>,
            value: Option<Vec<(usize, usize)>>,
        }
        fn do_lcs_test(test: LcsTest) {
            let mut tokens = Tokens::new();
            let left_toks = chars_to_toks(&mut tokens, &test.left);
            let right_toks = chars_to_toks(&mut tokens, &test.right);
            let lcs = CommonRange::flatten(&kc_lcs(&tokens, &left_toks, &right_toks));
            assert_eq!(lcs.len(), naive_lcs_length(&left_toks, &right_toks));
            match test.length {
                Some(len) => assert_eq!(lcs.len(), len),
                None => {}
            }
            match test.value {
                Some(value) => assert_eq!(lcs, value),
                None => {}
            }
        }
        macro_rules! lcs_test {
            ($($name:ident: $value:expr,)*) => {
            $(
                #[test]
                fn $name() {
                    do_lcs_test($value);
                }
            )*
            }
        }
        lcs_test! {
            tiny: LcsTest{left: "a", right: "a", length: Some(1), value: Some(vec![(0, 0)])},
            left_empty: LcsTest{left: "", right: "a", length: Some(0), value: Some(vec![])},
            right_empty: LcsTest{left: "a", right: "", length: Some(0), value: Some(vec![])},
            trivial: LcsTest{left: "aaaaaa", right: "aaaaaa", length: Some(6), value: None},
            example: LcsTest{
                left: "abacad",
                right: "bacada",
                length: Some(5),
                value: Some(vec![(1, 0), (2, 1), (3, 2), (4, 3), (5, 4)])},
            only_prefix: LcsTest{left: "ab", right: "ac", length: Some(1), value: Some(vec![(0, 0)])},
            middle_run: LcsTest{left: "abbbba", right: "cbbbbc", length: Some(4), value: Some(vec![(1, 1), (2, 2), (3,3), (4,4)])},
        }
        #[test]
        fn test_kc_lcs() {
            let mut tokens = Tokens::new();
            let s = "abacad";
            s.chars().for_each(|c| {
                tokens.get_token(vec![c as u8]);
            });
            let a = tokens.get_token(b"a".to_vec());
            let b = tokens.get_token(b"b".to_vec());
            let c = tokens.get_token(b"c".to_vec());
            let d = tokens.get_token(b"d".to_vec());
            let lcs = CommonRange::flatten(&kc_lcs(&tokens, &vec![a, b, c], &vec![b, a, c]));
            assert_eq!(lcs, vec![(1, 0), (2, 2)]);

            let lcs = CommonRange::flatten(&kc_lcs(
                &tokens,
                &vec![a, b, a, c, a, d],
                &vec![b, a, c, a, d, a],
            ));
            assert_eq!(lcs, vec![(1, 0), (2, 1), (3, 2), (4, 3), (5, 4)]);
        }
    }
}
