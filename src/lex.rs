use crate::token::{Parsed, Token, Tokens};
/// This module contains main lexer entrypoints.
use std::io::{self, BufRead};

/// Parses the contents of a reader into line based tokens.
pub fn lex_lines(tokens: &mut Tokens, r: &mut impl BufRead) -> io::Result<Parsed> {
    let mut lines: Vec<Token> = Vec::new();
    let mut starts: Vec<usize> = Vec::new();
    let mut offset: usize = 0;
    let mut done: bool = false;
    let mut line: Vec<u8> = Vec::with_capacity(100);
    while !done {
        // Our lines include the newline character, so both \r\n and \n are supported
        let bytes_read = r.read_until(b'\n', &mut line)?;
        if bytes_read == 0 {
            done = true;
            // Always add an offset for the end of the file
            starts.push(offset);
        } else {
            let trailing_newline = line.last() == Some(&b'\n');
            if trailing_newline {
                line.pop(); // remove the trailiing newline
            } else {
                // To mark the lack of a newline we add a special trailing token with an embedded newline
                line.extend(b"\n\\ No newline at end of file");
            } // remove the newline
            let token = tokens.get_token_ref(&line);
            line.clear();
            lines.push(token);
            starts.push(offset);
            offset += bytes_read;
        }
    }
    Ok(Parsed {
        tokens: lines,
        starts,
    })
}

pub fn fallback_lexer(tokens: &mut Tokens, contents: &[u8]) -> io::Result<Parsed> {
    let mut toks: Vec<Token> = Vec::new();
    let mut starts: Vec<usize> = Vec::new();
    let mut offset: usize = 0;
    let mut current_class = CharClass::Word;

    for (index, next_class) in contents
        .iter()
        .map(|c: &u8| classify(*c))
        .enumerate()
        .chain(std::iter::once((contents.len(), CharClass::Punctuation)))
    {
        // Punctuation is always a token boundary
        if next_class != current_class || next_class == CharClass::Punctuation {
            static SPACE: &[u8; 1] = b" ";
            let token_image = match current_class {
                CharClass::Whitespace => SPACE,
                _ => &contents[offset..index],
            };
            let token = tokens.get_token_ref(token_image);
            toks.push(token);
            starts.push(offset);
            offset = index;
            current_class = next_class;
        }
    }
    starts.push(contents.len());
    Ok({
        Parsed {
            tokens: toks,
            starts,
        }
    })
}

/// A lexer that is mostly useful for tests
pub fn lex_characters(tokens: &mut Tokens, s: &str) -> Vec<Token> {
    let mut toks = Vec::with_capacity(s.len());
    s.chars().for_each(|c| {
        toks.push(tokens.get_token(vec![c as u8]));
    });
    toks
}

#[derive(Debug, Eq, PartialEq)]
enum CharClass {
    Word,
    Punctuation,
    Whitespace,
}

fn classify(c: u8) -> CharClass {
    match c {
        b'(' | b')' | b'{' | b'}' | b'[' | b']' | b'<' | b'>' | b'.' | b',' | b';' | b':'
        | b'?' | b'!' | b'-' | b'+' | b'*' | b'/' | b'%' | b'^' | b'&' | b'|' | b'=' | b'~'
        | b'@' | b'#' | b'$' | b'_' => CharClass::Punctuation,
        b' ' | b'\t' | b'\n' | b'\r' => CharClass::Whitespace,
        _ => CharClass::Word,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_lex_lines() {
        let mut tokens = Tokens::new();
        let mut r = io::BufReader::new("Hello, world!\nHow are you?\n".as_bytes());
        let parsed = lex_lines(&mut tokens, &mut r).unwrap();
        assert_eq!(parsed.tokens, vec![0, 1]);
        assert_eq!(
            tokens.get_token_images(&parsed.tokens),
            vec!["Hello, world!", "How are you?"]
        );
        assert_eq!(parsed.starts, vec![0, 14, 27]);
    }

    #[test]
    fn test_fallback_lexer() {
        let mut tokens = Tokens::new();
        let mut r = "Hello, world!\nHow are you?\n".as_bytes();
        let parsed = fallback_lexer(&mut tokens, &mut r).unwrap();
        assert_eq!(parsed.tokens, vec![0, 1, 2, 3, 4, 2, 5, 2, 6, 2, 7, 8, 2]);
        assert_eq!(
            tokens.get_token_images(&parsed.tokens),
            vec!["Hello", ",", " ", "world", "!", " ", "How", " ", "are", " ", "you", "?", " "]
        );

        assert_eq!(
            parsed.starts,
            vec![0, 5, 6, 7, 12, 13, 14, 17, 18, 21, 22, 25, 26, 27]
        );
    }
}
