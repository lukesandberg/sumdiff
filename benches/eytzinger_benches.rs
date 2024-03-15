use std::collections::VecDeque;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion, Throughput};

use rand::{thread_rng, Rng};
use sumdifflib::eytzinger::{eytzinger_search, eytzingerize};
/// Transforms the input slice into EYT order using a BFS and a queue with a scratch buffer
/// This require a lot of scratch memory which is unfortunate, however it is a very simple algorithm
///
/// From benchmarking we can see this is the slowest of the three algorithms taking about 2.5x longer
fn bfs_eyt(v: &[usize]) -> Vec<usize> {
    let mut output = Vec::with_capacity(v.len() + 1);
    output.push(0);
    if v.len() >= 1 {
        let mut queue = VecDeque::with_capacity((v.len().next_power_of_two() / 2) as usize);
        queue.push_back((0 as usize, v.len()));
        while let Some((lower, upper)) = queue.pop_front() {
            let size = upper - lower;

            let mid = lower + size / 2;

            output.push(v[mid]);
            if mid != lower {
                queue.push_back((lower, mid));
            }
            if mid + 1 != upper {
                queue.push_back((mid + 1, upper));
            }
        }
    }
    output
}

/// Similar to the main implementation but avoids passing around mutable references
/// to `input_index`.
///
/// Shockingly this is also quite slow, taking  1.5x-2x longer than the main implementation
fn recursive_eyt_rv(v: &[usize]) -> Vec<usize> {
    fn recur(
        output_index: usize,
        mut input_index: usize,
        v: &[usize],
        output: &mut [usize],
    ) -> usize {
        let starting_input_index = input_index;
        if output_index < output.len() {
            let next_output_index = output_index * 2;
            input_index += recur(next_output_index, input_index, v, output);
            output[output_index] = v[input_index];
            input_index += 1;
            input_index += recur(next_output_index + 1, input_index, v, output);
        }
        input_index - starting_input_index
    }
    let mut output = vec![0; v.len() + 1];
    let final_index = recur(1, 0, v, &mut output);
    debug_assert!(final_index == v.len());
    output
}

const SIZE: &[usize] = &[
    1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072,
];
fn eyt_construction_bench(c: &mut Criterion) {
    {
        let mut group = c.benchmark_group("recursive_eyt_rv");
        for size in SIZE {
            group.throughput(Throughput::Elements(*size as u64));
            group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, size| {
                let input: Vec<usize> = (0..*size).collect();
                b.iter(|| recursive_eyt_rv(&input));
            });
        }
    }
    {
        let mut group = c.benchmark_group("eytzingerize");
        for size in SIZE {
            group.throughput(Throughput::Elements(*size as u64));
            group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, size| {
                let input: Vec<usize> = (0..*size).collect();
                b.iter(|| eytzingerize(&input));
            });
        }
    }
    {
        let mut group = c.benchmark_group("bfs_eyt");
        for size in SIZE {
            group.throughput(Throughput::Elements(*size as u64));
            group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, size| {
                let input: Vec<usize> = (0..*size).collect();
                b.iter(|| bfs_eyt(&input));
            });
        }
    }
}

/// Compares the performance of searching in an EYT array and a binary search
/// the eytzinger_search is consistently ~3x as fast as the binary search.
/// This comparision is not entirely fair since the algorithms do slightly different things
/// (one returns values, the other returns indexes), one returns out of bounds indexes the
/// other returns a sentinel. But this should be good enough to justify the usage.
fn search_eyt(c: &mut Criterion) {
    let vecs: Vec<(Vec<usize>, Vec<usize>)> = SIZE
        .iter()
        .map(|&size| {
            let sorted = (1..=size).map(|v| v * 2 + 1).collect::<Vec<usize>>();
            let eyt = eytzingerize(&sorted);
            (sorted, eyt)
        })
        .collect();
    {
        let mut group = c.benchmark_group("search_eyt");
        for (sorted, eyt) in &vecs {
            let size = sorted.len() as u64;
            let max = sorted.last().unwrap();
            group.throughput(Throughput::Elements(size));

            group.bench_with_input(BenchmarkId::from_parameter(size), &size, |b, _| {
                let mut rng = thread_rng();
                let cloned_eyt = eyt.clone(); // Clone the eyt variable
                b.iter_batched(
                    move || rng.gen_range(1..*max + 1),
                    move |i| eytzinger_search(&cloned_eyt, i), // Use the cloned_eyt variable
                    BatchSize::SmallInput,
                );
            });
        }
    }
    {
        let mut group = c.benchmark_group("search_binary");
        for (sorted, _) in &vecs {
            let size = sorted.len() as u64;
            let max = sorted.last().unwrap();
            group.throughput(Throughput::Elements(size));

            group.bench_with_input(BenchmarkId::from_parameter(size), &size, |b, _| {
                let mut rng = thread_rng();
                let cloned_sorted = sorted.clone(); // Clone the eyt variable
                b.iter_batched(
                    move || rng.gen_range(1..*max + 1),
                    move |i| match &cloned_sorted.binary_search(&i) {
                        Ok(i) => *i,
                        Err(i) => *i,
                    },
                    BatchSize::SmallInput,
                );
            });
        }
    }
}

criterion_group!(eytzinger_construction_benches, eyt_construction_bench,);
criterion_group!(search_benches, search_eyt);
criterion_main!(eytzinger_construction_benches, search_benches);
