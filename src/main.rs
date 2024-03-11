mod display;
mod kc;
mod lex;
mod token;
use std::fs;
use std::io;
use std::{env, process};

use crate::token::Tokens;

struct FilePair {
    left: String,
    right: String,
}
fn main() {
    let file_pair = parse_args(env::args()).unwrap_or_else(|err| {
        eprintln!("Usage: sumdiff <left> <right>: {err}");
        process::exit(1)
    });
    fn buf_reader(file: &String) -> io::BufReader<fs::File> {
        io::BufReader::new(
            fs::OpenOptions::new()
                .read(true)
                .open(file)
                .unwrap_or_else(|err| {
                    eprintln!("Unable to open {file}: {err}");
                    process::exit(1)
                }),
        )
    }
    let mut line_tokens = Tokens::new();
    let mut left_reader = buf_reader(&file_pair.left);
    let left_lines = lex::lex_lines(&mut line_tokens, &mut left_reader).unwrap_or_else(|err| {
        eprintln!("Unable to read {0}: {err}", file_pair.left);
        process::exit(1)
    });
    let mut right_reader = buf_reader(&file_pair.right);
    let right_lines = lex::lex_lines(&mut line_tokens, &mut right_reader).unwrap_or_else(|err| {
        eprintln!("Unable to read {0}: {err}", file_pair.right);
        process::exit(1)
    });
    let matches = kc::kc_lcs(&line_tokens, &left_lines.tokens, &right_lines.tokens);
    display::display_diff(&line_tokens.printer(), &left_lines, &right_lines, &matches)
        .unwrap_or_else(|err| {
            eprintln!("Error while displaying diff: {err}");
            process::exit(1)
        });
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
