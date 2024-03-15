use std::mem::{transmute, MaybeUninit};
// Inspired by https://en.algorithmica.org/hpc/data-structures/binary-search/

/// Returns an Eytzinger layout of the input vector as a copy.
pub fn eytzingerize(v: &[usize]) -> Vec<usize> {
    fn recur(
        output_index: usize,
        input_index: &mut usize,
        v: &[usize],
        output: &mut [MaybeUninit<usize>],
    ) {
        if output_index < output.len() {
            let next_output_index = output_index * 2;
            recur(next_output_index, input_index, v, output);
            // SAFETY: The compiler can eliminate the bounds check on output but not `v` :(
            // However this is safe because the recursion ensures that we never increment input_index past v.len()
            // and there is a debug check at the end to check that.
            output[output_index] = MaybeUninit::new(unsafe { *v.get_unchecked(*input_index) });
            *input_index += 1;
            // This will be turned into a tail call.
            recur(next_output_index + 1, input_index, v, output);
        }
    }
    // TODO: ensure that this buffer is cache line aligned
    // SAFETY: We use maybe uninitialized to avoid the overhead of zeroing the buffer since our recurrence will
    // assign every element.
    let mut output = vec![MaybeUninit::<usize>::uninit(); v.len() + 1];
    output[0] = MaybeUninit::new(0);
    let mut input_index = 0;
    recur(1, &mut input_index, v, &mut output);
    debug_assert!(input_index == v.len());
    unsafe { transmute(output) }
}

/// Returns the index of the first set bit in x (1 based indexing)
#[inline]
fn find_first_set(x: usize) -> usize {
    (x.trailing_zeros() + 1) as usize
}

/// Returns
pub fn eytzinger_search(v: &[usize], target: usize) -> usize {
    let mut k = 1;
    let n = v.len();
    // This has surprising results for n == 0 but we are guaranteed that n > 0
    // in that case this is exact the same as ilog2() but without the brances handling zero.
    let iters = (usize::BITS - 1) - n.leading_zeros();
    for _i in 0..iters {
        // TODO: add a prefetching instruction here too fetch the next cache line
        // SAFETY: We are guaranteed that k < n at this point because the loop is bound
        // by the floor(log2(n)) and k doubles in each iteration.
        let cmp = unsafe { *v.get_unchecked(k) } < target;

        // rust guarantees that booleans are represented as 0 and 1 for true and false.
        k = 2 * k + cmp as usize;
    }
    // We unroll the last iteration to make the loop branch fully predictable.
    // This branch is unpredictable however it will generally be compiled to a conditional move
    let last_index = if k < n { k } else { 0 };
    // SAFETY: The check above guarnatees that last_index is <n but we also need to ensure
    // that n>0 so that we don't panic.  This is guaranteed by the construction of the etz array.
    let cmp = unsafe { *v.get_unchecked(last_index) } < target;
    k = 2 * k + cmp as usize;
    k >>= find_first_set(!k);
    v[k]
}

#[cfg(test)]
mod tests {

    use super::*;

    fn check_etz(output: Vec<usize>, input: Vec<usize>) {
        let actual: Vec<usize> = eytzingerize(&input);
        assert_eq!(output, actual);
    }
    #[test]
    fn test_eytzingerize() {
        // 1
        check_etz(
            vec![0, 1], //
            vec![1],
        );
        // 2
        check_etz(
            vec![0, 2, 1], //
            vec![1, 2],
        );
        // 3
        check_etz(
            vec![0, 2, 1, 3], //
            vec![1, 2, 3],
        );
        // 4
        check_etz(
            vec![0, 3, 2, 4, 1], //
            vec![1, 2, 3, 4],
        );
        // 5
        check_etz(
            vec![0, 4, 2, 5, 1, 3], //
            vec![1, 2, 3, 4, 5],
        );
        // 6
        check_etz(
            vec![0, 4, 2, 6, 1, 3, 5], //
            vec![1, 2, 3, 4, 5, 6],
        );
        // 7
        check_etz(
            vec![0, 4, 2, 6, 1, 3, 5, 7], //
            vec![1, 2, 3, 4, 5, 6, 7],
        );
        // 8
        check_etz(
            vec![0, 5, 3, 7, 2, 4, 6, 8, 1], //
            vec![1, 2, 3, 4, 5, 6, 7, 8],
        );
        // 9
        check_etz(
            vec![0, 6, 4, 8, 2, 5, 7, 9, 1, 3], //
            vec![1, 2, 3, 4, 5, 6, 7, 8, 9],
        );
        // 10
        check_etz(
            vec![0, 7, 4, 9, 2, 6, 8, 10, 1, 3, 5], //
            vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        );
        // 14
        check_etz(
            vec![0, 8, 4, 12, 2, 6, 10, 14, 1, 3, 5, 7, 9, 11, 13], //
            vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
        );
        // 15
        check_etz(
            vec![0, 8, 4, 12, 2, 6, 10, 14, 1, 3, 5, 7, 9, 11, 13, 15], //
            vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        );
    }

    #[test]
    fn test_ffs() {
        assert_eq!(find_first_set(0b01), 1);
        assert_eq!(find_first_set(0b010), 2);
        assert_eq!(find_first_set(0b01000), 4);
        assert_eq!(find_first_set(0b01001), 1);
    }

    #[test]
    fn test_eytzinger_search() {
        let sorted = (0..100).map(|v| v * 2 + 1).collect::<Vec<usize>>();
        let etz = eytzingerize(&sorted);
        // exact matches
        for m in sorted.iter() {
            let bin_search = sorted[sorted.binary_search(m).unwrap()];
            assert_eq!(bin_search, *m);
            assert_eq!(eytzinger_search(&etz, *m), *m);
        }
        // exact misses
        for m in (0..100).map(|v| v * 2).collect::<Vec<usize>>().iter() {
            let bin_search = sorted[sorted.binary_search(m).unwrap_err()];
            assert_eq!(bin_search, *m + 1);
            assert_eq!(eytzinger_search(&etz, *m), *m + 1);
        }
        // matches past the end are different and basically depend on a sentinel value
        // for binary search it just points past the end of the array
        // for etz search it points at the zero value
        let bin_search = sorted.binary_search(&200).unwrap_err();
        assert_eq!(bin_search, 100); // points past the end of the array
        assert_eq!(eytzinger_search(&etz, 200), 0); // points at the zero value
    }

    #[test]
    fn test_eytzinger_search_small() {
        let sorted: Vec<usize> = vec![2];
        let etz = eytzingerize(&sorted);
        assert_eq!(eytzinger_search(&etz, 2), 2);
        assert_eq!(eytzinger_search(&etz, 1), 2);
        assert_eq!(eytzinger_search(&etz, 3), 0);
    }
}
