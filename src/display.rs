use std::io::Error;
use std::io::StdoutLock;
use std::io::Write;

use crate::token::{Parsed, Token, TokenPrinter};

pub fn display_diff(
    printer: &TokenPrinter,
    left: &Parsed,
    right: &Parsed,
    matches: &[(usize, usize)],
) -> Result<(), Error> {
    let mut left_index = 0;
    let mut right_index = 0;
    let mut stdout = std::io::stdout().lock();
    for (left_match, right_match) in matches {
        flush_tokens(
            printer,
            &left.tokens,
            left_index,
            *left_match,
            b"-",
            &mut stdout,
        )?;
        flush_tokens(
            printer,
            &right.tokens,
            right_index,
            *right_match,
            b"+",
            &mut stdout,
        )?;

        // this means we matched some lines
        stdout.write_all(b" ")?;
        stdout.write_all(&printer.print(left.tokens[*left_match]))?;
        left_index = *left_match + 1;
        right_index = *right_match + 1;
    }
    flush_tokens(
        printer,
        &left.tokens,
        left_index,
        left.tokens.len(),
        b"-",
        &mut stdout,
    )?;
    flush_tokens(
        printer,
        &right.tokens,
        right_index,
        right.tokens.len(),
        b"+",
        &mut stdout,
    )?;
    Ok(())
}

fn flush_tokens(
    printer: &TokenPrinter,
    tokens: &[Token],
    start: usize,
    end: usize,
    prefix: &[u8],
    stdout: &mut StdoutLock,
) -> Result<(), Error> {
    for token in &tokens[start..end] {
        stdout.write_all(prefix)?;
        stdout.write_all(&printer.print(*token))?;
    }
    Ok(())
}
