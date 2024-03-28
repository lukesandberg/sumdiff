use criterion::{criterion_group, criterion_main, Criterion};
use rand::{thread_rng, Rng};
use sumdifflib::{
    dijkstra::dijkstra,
    kc::kc_lcs,
    lcs_utils::naive_lcs_length,
    token::{Token, Tokens},
};

const SIZE: &[usize] = &[1024, 65536, 1024 * 1024];

#[derive(Debug)]
enum EditSize {
    Small,
    Medium,
    Large,
}

impl EditSize {
    fn edit(&self, input: &[Token]) -> Vec<Token> {
        let mut right = input.to_vec();
        let size = input.len();
        let mut rnd = thread_rng();

        let edit_size = match self {
            EditSize::Small => {
                // edit 1% of the tokens
                size / 10
            }
            EditSize::Medium => {
                // edit 5% of the tokens
                size / 20
            }
            EditSize::Large => {
                // edit 10% of the tokens
                size / 10
            }
        };
        for _ in 0..std::cmp::max(1, edit_size) {
            let idx = rnd.gen_range(0..size);
            right[idx] = right[rnd.gen_range(0..size)];
        }
        return right.to_vec();
    }
}

#[derive(Debug)]
enum AlphabetSize {
    Small,
    Medium,
    Large,
}

impl AlphabetSize {
    fn size(&self) -> usize {
        match self {
            AlphabetSize::Small => 10,
            AlphabetSize::Medium => 100,
            AlphabetSize::Large => 1000,
        }
    }
}

fn lcs_benchmark(c: &mut Criterion) {
    for (size_index, &size) in SIZE.iter().enumerate() {
        for alphabet_size in [
            AlphabetSize::Small,
            AlphabetSize::Medium,
            AlphabetSize::Large,
        ]
        .iter()
        {
            let mut tokens = Tokens::new();
            let alpha_size = alphabet_size.size();
            let token_universe = (0..alpha_size)
                .map(|v: usize| tokens.get_token(v.to_string().into_bytes()))
                .collect::<Vec<Token>>();
            let mut rnd = thread_rng();
            let dist = rand::distributions::WeightedIndex::new(
                (1..=alpha_size)
                    .map(|i| std::cmp::max(1, alpha_size * (i / 10)))
                    .collect::<Vec<usize>>(),
            )
            .unwrap();
            let left = (0..size)
                .map(|_| token_universe[rnd.sample(&dist)])
                .collect::<Vec<Token>>();
            for edit_size in [EditSize::Small, EditSize::Medium, EditSize::Large].iter() {
                let right = edit_size.edit(&left);
                {
                    let mut group = c.benchmark_group(format!(
                        "lcs/alpha:{:?}/edit:{:?}/{:?}",
                        alphabet_size, edit_size, size
                    ));
                    group.throughput(criterion::Throughput::Elements(size as u64));
                    group.bench_function("kc", |b| {
                        b.iter(|| {
                            kc_lcs(&tokens, &left, &right);
                        });
                    });
                    if size_index == 0 {
                        // dijkstra is way too slow for large inputs
                        group.bench_function("dijkstra", |b| {
                            b.iter(|| {
                                dijkstra(&left, &right);
                            });
                        });
                    }
                    group.finish();
                }
                {
                    let mut group = c.benchmark_group(format!(
                        "lcs_len/alpha:{:?}/edit:{:?}/{:?}",
                        alphabet_size, edit_size, size
                    ));
                    group.throughput(criterion::Throughput::Elements(size as u64));
                    if size_index == 0 {
                        // The naive algorithm is too slow for large inputs
                        group.bench_function("naive", |b| {
                            b.iter(|| {
                                naive_lcs_length(&left, &right);
                            });
                        });
                    }
                    group.bench_function("meyers", |b| {
                        b.iter(|| {
                            dijkstra(&left, &right);
                        });
                    });
                    group.finish();
                }
            }
        }
    }
}

criterion_group!(benches, lcs_benchmark);
criterion_main!(benches);
