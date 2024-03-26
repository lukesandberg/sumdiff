use rand::{thread_rng, Rng};
use sumdifflib::{
    dijkstra::dijkstra,
    kc::kc_lcs,
    lcs_utils::check_is_lcs,
    lex::lex_characters,
    token::{CommonRange, Token, Tokens},
};
#[cfg(test)]
extern crate quickcheck;
#[cfg(test)]
#[macro_use(quickcheck)]
extern crate quickcheck_macros;

struct LcsTest {
    left: &'static str,
    right: &'static str,
    length: Option<usize>,
    value: Option<Vec<(usize, usize)>>,
}

fn do_lcs_test(test: LcsTest) {
    let mut tokens = Tokens::new();
    let left_toks = lex_characters(&mut tokens, &test.left);
    let right_toks = lex_characters(&mut tokens, &test.right);
    let kc_lcs = CommonRange::flatten(&kc_lcs(&tokens, &left_toks, &right_toks));
    check_is_lcs(&kc_lcs, &left_toks, &right_toks).unwrap();
    let dijkstra_lcs = CommonRange::flatten(&dijkstra(&left_toks, &right_toks));
    check_is_lcs(&dijkstra_lcs, &left_toks, &right_toks).unwrap();
    match test.length {
        Some(len) => {
            assert_eq!(kc_lcs.len(), len);
            assert_eq!(dijkstra_lcs.len(), len);
        }
        None => {}
    }
    match test.value {
        Some(value) => {
            assert_eq!(kc_lcs, value);
            assert_eq!(dijkstra_lcs, value);
        }
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

// Shared test cases, these should only be used when there is a single unique LCS
lcs_test! {
    tiny: LcsTest{left: "a", right: "a", length: Some(1), value: Some(vec![(0, 0)])},
    left_empty: LcsTest{left: "", right: "a", length: Some(0), value: Some(vec![])},
    right_empty: LcsTest{left: "a", right: "", length: Some(0), value: Some(vec![])},
    non_trivial_zero_diff: LcsTest{left: "abcd", right: "efghi", length: Some(0), value: Some(vec![])},
    trivial: LcsTest{left: "aaaaaa", right: "aaaaaa", length: Some(6), value: None},
    example: LcsTest{
        left: "abacad",
        right: "bacada",
        length: Some(5),
        value: Some(vec![(1, 0), (2, 1), (3, 2), (4, 3), (5, 4)])},
    only_prefix: LcsTest{left: "ab", right: "ac", length: Some(1), value: Some(vec![(0, 0)])},
    middle_run: LcsTest{left: "abba", right: "cbbc", length: Some(2), value: Some(vec![(1, 1), (2, 2)])},
}
fn lex_nums(tokens: &mut Tokens, s: &[u8]) -> Vec<Token> {
    let mut toks = Vec::with_capacity(s.len());
    s.iter().for_each(|&c| {
        toks.push(tokens.get_token(vec![c]));
    });

    toks
}
#[quickcheck]
fn quickcheck_kc(left: Vec<u8>, right: Vec<u8>) -> bool {
    let mut tokens = Tokens::new();
    let left_toks = lex_nums(&mut tokens, &left);
    let right_toks = lex_nums(&mut tokens, &right);
    let kc_lcs = CommonRange::flatten(&kc_lcs(&tokens, &left_toks, &right_toks));
    check_is_lcs(&kc_lcs, &left_toks, &right_toks).is_ok()
}
#[quickcheck]
fn quickcheck_dijktra(left: Vec<u8>, right: Vec<u8>) -> bool {
    let mut tokens = Tokens::new();
    let left_toks = lex_nums(&mut tokens, &left);
    let right_toks = lex_nums(&mut tokens, &right);
    let dijkstra_lcs = CommonRange::flatten(&dijkstra(&left_toks, &right_toks));
    check_is_lcs(&dijkstra_lcs, &left_toks, &right_toks).is_ok()
}

const SIZE: &[usize] = &[10, 20, 64];

#[derive(Debug)]
enum EditSize {
    Small,
    Medium,
    Large,
}

#[test]
fn lcs_benchmark() {
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
            let kc_lcs = CommonRange::flatten(&kc_lcs(&tokens, &left, &right));
            check_is_lcs(&kc_lcs, &left, &right).unwrap();
            let dijkstra_lcs = CommonRange::flatten(&dijkstra(&left, &right));
            check_is_lcs(&dijkstra_lcs, &left, &right).unwrap();
        }
    }
}
