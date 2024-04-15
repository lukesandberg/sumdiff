use criterion::{criterion_group, criterion_main, Criterion};
use rand::{thread_rng, Rng};
use sumdifflib::{
    dijkstra::dijkstra,
    hyyro::hyyro_lcs_len,
    kc::kc_lcs,
    lcs_utils::naive_lcs_length,
    meyers::{meyers_lcs, meyers_lcs_length},
    token::{Token, Tokens},
};

const SIZE: &[usize] = &[1024, 2048, 4096, 8192, 16384];

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
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
                size / 100
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

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
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

#[derive(Debug, Clone, PartialEq, Eq)]
struct Scenario {
    left: Vec<Token>,
    right: Vec<Token>,
    edit_size: EditSize,
    alphabet_size: AlphabetSize,
    size: usize,
}

impl Scenario {
    fn sample_size(&self) -> usize {
        if self.size <= 1024 {
            100
        } else {
            10
        }
    }
}

fn lcs_benchmark(c: &mut Criterion) {
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
        for edit_size in [EditSize::Small, EditSize::Medium, EditSize::Large].iter() {
            let mut scenarios = Vec::new();
            for &size in SIZE {
                let left = (0..size)
                    .map(|_| token_universe[rnd.sample(&dist)])
                    .collect::<Vec<Token>>();
                let right = edit_size.edit(&left);
                scenarios.push(Scenario {
                    left: left.clone(),
                    right,
                    edit_size: *edit_size,
                    alphabet_size: *alphabet_size,
                    size,
                });
            }

            {
                let mut group = c.benchmark_group(format!(
                    "lcs/alpha:{:?}/edit:{:?}",
                    alphabet_size, edit_size
                ));
                for scenario in scenarios {
                    group.throughput(criterion::Throughput::Elements(scenario.size as u64));
                    group.sample_size(scenario.sample_size());
                    group.bench_function("kc", |b| {
                        b.iter(|| {
                            kc_lcs(&tokens, &scenario.left, &scenario.right);
                        });
                    });
                    if scenario.size == 1024 {
                        // dijkstra is way too slow for large inputs
                        group.bench_function("dijkstra", |b| {
                            b.iter(|| {
                                dijkstra(&scenario.left, &scenario.right);
                            });
                        });
                        // The naive algorithm is too slow for large inputs
                        group.bench_function("naive_len", |b| {
                            b.iter(|| {
                                naive_lcs_length(&scenario.left, &scenario.right);
                            });
                        });
                    }
                    group.bench_function("meyers_len", |b| {
                        b.iter(|| {
                            meyers_lcs_length(&scenario.left, &scenario.right);
                        });
                    });
                    group.bench_function("meyers", |b| {
                        b.iter(|| {
                            meyers_lcs(&scenario.left, &scenario.right);
                        });
                    });
                    group.bench_function("hyyro_len", |b| {
                        b.iter(|| {
                            hyyro_lcs_len(&tokens, &scenario.left, &scenario.right);
                        });
                    });
                }
                group.finish();
            }
        }
    }
}

criterion_group!(benches, lcs_benchmark);
criterion_main!(benches);
