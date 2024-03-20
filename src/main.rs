use std::fs;
use std::io;
use std::{env, process};
use sumdifflib::{
    display::{display_diff, DisplayConfig},
    hunk::matches_to_hunks,
    kc, lex,
    token::{Parsed, Tokens},
};

struct FilePair {
    left: String,
    right: String,
}

fn main() {
    let file_pair = parse_args(env::args()).unwrap_or_else(|err| {
        eprintln!("Usage: sumdiff <left> <right>: {err}");
        process::exit(1)
    });
    let mut line_tokens = Tokens::new();
    let left_lines = tokenize_or_die(&file_pair.left, &mut line_tokens);
    let right_lines = tokenize_or_die(&file_pair.right, &mut line_tokens);
    let matches = kc::kc_lcs(&line_tokens, &left_lines.tokens, &right_lines.tokens);
    let hunks = matches_to_hunks(
        3,
        &matches,
        left_lines.tokens.len(),
        right_lines.tokens.len(),
    );
    if hunks.is_empty() {
        // By conventions empty diffs produce no output
        return;
    }
    println!("--- {}", file_pair.left);
    println!("+++ {}", file_pair.right);
    display_diff(
        &DisplayConfig {
            insert_prefix: b"+",
            delete_prefix: b"-",
            common_prefix: b" ",
            token_suffix: b"\n",
        },
        &line_tokens.printer(),
        &left_lines,
        &right_lines,
        &hunks,
        &mut std::io::stdout().lock(),
    )
    .unwrap_or_else(|err| {
        eprintln!("Error while displaying diff: {err}");
        process::exit(1)
    });
}

fn tokenize_or_die(file_name: &String, line_tokens: &mut Tokens) -> Parsed {
    let file = fs::OpenOptions::new()
        .read(true)
        .open(file_name)
        .unwrap_or_else(|err| {
            eprintln!("Unable to open {file_name}: {err}");
            process::exit(1)
        });
    let mut reader = io::BufReader::new(file);
    lex::lex_lines(line_tokens, &mut reader).unwrap_or_else(|err| {
        eprintln!("Unable to read {file_name}: {err}");
        process::exit(1)
    })
}

fn parse_args(mut args: impl Iterator<Item = String>) -> Result<FilePair, &'static str> {
    args.next(); // skip the program name
    let left = match args.next() {
        Some(arg) => arg,
        None => return Err("Expected 2 arguments, got 0"),
    };
    let right = match args.next() {
        Some(arg) => arg,
        None => return Err("Expected 2 arguments, got 1"),
    };
    Ok(FilePair { left, right })
}
