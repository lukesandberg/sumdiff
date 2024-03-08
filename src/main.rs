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
    // We are missing the implementation of the diff function, so just do some dummy lexing
    let mut left = buf_reader(&file_pair.left);
    let mut right = buf_reader(&file_pair.right);
    let mut line_tokens = Tokens::new();
    println!(
        "left-lines: {:?}",
        lex::lex_lines(&mut line_tokens, &mut left)
    );
    println!(
        "right-lines: {:?}",
        lex::lex_lines(&mut line_tokens, &mut right)
    );
    let mut fallback_tokens = Tokens::new();
    println!(
        "left-tokens: {:?}",
        lex::fallback_lexer(&mut fallback_tokens, &fs::read(&file_pair.left).unwrap())
    );
    println!(
        "right-tokens: {:?}",
        lex::fallback_lexer(&mut fallback_tokens, &fs::read(&file_pair.right).unwrap())
    );
}

fn parse_args(mut args: impl Iterator<Item = String>) -> Result<FilePair, &'static str> {
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
