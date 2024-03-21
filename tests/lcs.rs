use sumdifflib::{
    kc::kc_lcs,
    token::{CommonRange, Token, Tokens},
};

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
