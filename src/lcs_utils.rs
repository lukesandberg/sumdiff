use crate::token::{CommonRange, Token};
use std::mem::{transmute, MaybeUninit};

/// Peforms a simple 'counting' sort on the token array, but returns the indices of the sorted token locations
/// rather than sorting the tokens themselves.
/// So output[0] is the index of the first occurrence of the smallest token, and so on.
pub fn counting_sort(arr: &[Token], max_token: usize) -> Vec<usize> {
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

/// Returns the length of the common prefix of `a` and `b`
fn common_prefix<T: PartialEq + PartialOrd>(a: &[T], b: &[T]) -> usize {
    let mut i = 0;
    while i < a.len() && i < b.len() && a[i] == b[i] {
        i += 1;
    }
    i
}

// Returns the length of the common suffix of `a` and `b`
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

pub struct Trimmed<'a> {
    pub left: &'a [Token],
    pub right: &'a [Token],
    pub prefix: usize,
    pub suffix: usize,
    pub matches: Vec<CommonRange>,
}

/// Removes the common prefixes and suffixes
///
/// If the common prefixes and suffixes account for the full set of common ranges, just return them.
/// Otherwise return the remaining ranges and the length of the common prefix and suffix.
pub fn remove_suffixes_and_prefixes<'a>(
    left: &'a [Token],
    right: &'a [Token],
) -> Result<Vec<CommonRange>, Trimmed<'a>> {
    let mut n = left.len();
    let mut m = right.len();
    if n == 0 || m == 0 {
        return Ok(Vec::<CommonRange>::new());
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
            return Ok(matches);
        }
    }
    // wait to append the suffix matches until after the main loop
    let suffix = common_suffix(left, right);
    let left = &left[prefix..n - suffix];
    let right = &right[prefix..m - suffix];
    n = left.len();
    m = right.len();
    // If the suffix reduces both inputs to <=1 we are done since we have already compared prefixes and suffixes
    if n <= 1 && m <= 1 || n == 0 || m == 0 {
        if suffix > 0 {
            matches.push(CommonRange {
                left_start: n + prefix,
                right_start: m + prefix,
                length: suffix,
            });
        }
        return Ok(matches);
    }
    Err(Trimmed {
        left,
        right,
        matches,
        prefix,
        suffix,
    })
}

/// A branchfree implementation of binary search over a subrange
/// It is the callers responsibility to ensure that offset and end are in range
///
/// For our usecase this significantly outperforms the standard library binary search.
fn binary_search_range(arr: &[usize], offset: usize, end: usize, x: usize) -> usize {
    debug_assert!(offset < end && end <= arr.len() && offset < arr.len());
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
pub fn exponential_search_range(arr: &[usize], offset: usize, end: usize, x: usize) -> usize {
    debug_assert!(offset < end && end <= arr.len() && offset < arr.len());
    let mut bound = 1;
    let mut probe = offset + bound;
    while probe < end && unsafe { *arr.get_unchecked(probe) } < x {
        bound *= 2;
        probe = offset + bound;
    }
    let r = binary_search_range(
        arr,
        // start at the previous search location which we know is <=x
        offset + bound / 2,
        // end needs to be at most the size of the array or one past the last test which is >= x
        std::cmp::min(end, probe + 1),
        x,
    );
    debug_assert!(
        /* past the end */
        (r == end && arr[r - 1] < x)
        || /* before the beginning */ (r == offset && arr[r] >= x)
        || /* in the middle */ (arr[r] >= x && (r+1 == end || arr[r + 1] > x)),
        "offset={}, end={}, x={}, r={}, arr={:?}",
        offset,
        end,
        x,
        r,
        &arr[offset..end]
    );
    r
}

/// The trivial N^2 recurrence
fn naive_lcs_length(a: &[Token], b: &[Token]) -> usize {
    let n = a.len();
    let m = b.len();
    let mut cur = vec![0; m + 1];
    let mut prev = vec![0; m + 1];
    for i in 1..=n {
        for j in 1..=m {
            if a[i - 1] == b[j - 1] {
                cur[j] = prev[j - 1] + 1;
            } else {
                cur[j] = std::cmp::max(prev[j], cur[j - 1]);
            }
        }
        // swap rows
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[b.len()]
}

fn check_is_common_subsequence(
    lcs: &[(usize, usize)],
    left: &[Token],
    right: &[Token],
) -> Result<(), String> {
    let mut prev = None;
    for (l, r) in lcs {
        if let Some((pl, pr)) = prev {
            if pl >= *l || pr >= *r {
                return Err(format!(
                    "LCS sequence should be strictly increasing: ({}, {}) followed by ({}, {})",
                    pl, pr, l, r
                ));
            }
        }
        if left[*l] != right[*r] {
            return Err(format!(
                "LCS sequence should be a common subsequence: (left[{}] = {}, right[{}] = {})",
                l, left[*l], r, right[*r]
            ));
        }
        prev = Some((*l, *r));
    }
    Ok(())
}

/// Validate an LCS
pub fn check_is_lcs(lcs: &[(usize, usize)], left: &[Token], right: &[Token]) -> Result<(), String> {
    let expected_length = naive_lcs_length(left, right);
    if expected_length != lcs.len() {
        return Err(format!(
            "Expected length {} but got {}",
            expected_length,
            lcs.len()
        ));
    }
    check_is_common_subsequence(lcs, left, right)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::token::Tokens;

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

        // subrange tests
        assert_eq!(exponential_search_range(&arr, 1, 2, 0), 1);
        assert_eq!(exponential_search_range(&arr, 1, 2, 1), 1);
        assert_eq!(exponential_search_range(&arr, 1, 2, 2), 2);
        assert_eq!(exponential_search_range(&arr, 1, 2, 3), 2);
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
}
