use criterion::{criterion_group, criterion_main, Criterion};
use rand::{thread_rng, Rng};
use sumdifflib::{
    kc::kc_lcs,
    token::{Token, Tokens},
};

const SIZE: &[usize] = &[1024, 65536, 1024 * 1024];

#[derive(Debug)]
enum EditSize {
    Small,
    Medium,
    Large,
}
fn lcs_benchmark(c: &mut Criterion) {
    for edit_size in [EditSize::Small, EditSize::Medium, EditSize::Large].iter() {
        for &size in SIZE {
            let mut tokens = Tokens::new();
            let token_universe = (0..size / 2)
                .map(|v: usize| tokens.get_token(v.to_string().into_bytes()))
                .collect::<Vec<Token>>();

            let mut rnd = thread_rng();
            let dist = rand::distributions::WeightedIndex::new(
                (1..=size / 2)
                    .map(|i| std::cmp::max(1, (size / 2) * (i / 10)))
                    .collect::<Vec<usize>>(),
            )
            .unwrap();
            let left = (0..size)
                .map(|_| token_universe[rnd.sample(&dist)])
                .collect::<Vec<Token>>();
            let mut right = left.clone();
            match edit_size {
                EditSize::Small => {
                    // edit 1% of the tokens
                    for _ in 0..std::cmp::max(1, size / 100) {
                        let idx = rnd.gen_range(0..size);
                        right[idx] = token_universe[rnd.sample(&dist)];
                    }
                }
                EditSize::Medium => {
                    // edit 5% of the tokens
                    for _ in 0..std::cmp::max(1, size / 20) {
                        let idx = rnd.gen_range(0..size);
                        right[idx] = token_universe[rnd.sample(&dist)];
                    }
                }
                EditSize::Large => {
                    // edit 10% of the tokens
                    for _ in 0..std::cmp::max(1, size / 10) {
                        let idx = rnd.gen_range(0..size);
                        right[idx] = token_universe[rnd.sample(&dist)];
                    }
                }
            }
            let mut group = c.benchmark_group(format!("lcs/{:?}/{:?}", edit_size, size));
            group.throughput(criterion::Throughput::Elements(size as u64));
            group.bench_function("kc", |b| {
                b.iter(|| {
                    kc_lcs(&tokens, &left, &right);
                });
            });
        }
    }
}

criterion_group!(benches, lcs_benchmark);
criterion_main!(benches);
