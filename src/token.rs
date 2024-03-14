use std::collections::HashMap;

/// A Tokens is a collection of unique strings to unique token identifiers
#[derive(Debug)]
pub struct Tokens {
    // TODO(luke): Consider
    // 1. an arena allocator for the keys
    // 2. a trivial mechanism for looking up a key by value
    universe: HashMap<Vec<u8>, Token>,
}
impl Tokens {
    #[must_use]
    pub fn new() -> Self {
        Tokens {
            universe: HashMap::new(),
        }
    }
    // A token that is guaranteed to be > all known tokens
    pub fn token_upper_bound(&self) -> usize {
        self.universe.len()
    }
    /// Get the token for a string, creating a new one if necessary.
    pub fn get_token_ref(&mut self, s: &Vec<u8>) -> Token {
        match self.universe.get(s) {
            Some(token) => *token,
            None => {
                let id = self.get_next_token();
                // Since we only have an immutable reference we need to copy to insert
                self.universe.insert(s.clone(), id);
                id
            }
        }
    }
    /// Get the token for a string, creating a new one if necessary.
    pub fn get_token(&mut self, s: Vec<u8>) -> Token {
        let next_token = self.get_next_token();
        match self.universe.entry(s) {
            std::collections::hash_map::Entry::Occupied(entry) => *entry.get(),
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(next_token);
                next_token
            }
        }
    }

    fn get_next_token(&self) -> Token {
        (self.universe.len()) as Token
    }

    /// Get the string for a token, used only for testing.
    /// This is an O(n) operation in the size of the universe and the tokens
    #[cfg(test)]
    pub fn get_token_images(&self, tokens: &Vec<Token>) -> Vec<String> {
        let printer = self.printer();

        tokens
            .iter()
            .map(|&t| String::from_utf8_lossy(&printer.print(t)).to_string())
            .collect()
    }

    /// Returns a printer that can recover token images from tokens.
    pub fn printer<'a>(&'a self) -> TokenPrinter<'a> {
        let mut keys = vec![&[] as &[u8]; self.universe.len()];
        for (k, v) in self.universe.iter() {
            keys.as_mut_slice()[*v as usize] = k;
        }
        TokenPrinter { tokens: keys }
    }
}

pub struct TokenPrinter<'a> {
    tokens: Vec<&'a [u8]>,
}
impl<'a> TokenPrinter<'a> {
    pub fn print(&'a self, token: Token) -> Vec<u8> {
        self.tokens[token as usize].to_vec()
    }
    pub fn concat(&'a self, tokens: &Vec<Token>) -> Vec<u8> {
        let mut size: usize = 0;
        for token in tokens {
            size += self.tokens[*token as usize].len();
        }
        let mut result = Vec::with_capacity(size);
        for token in tokens {
            result.extend_from_slice(self.tokens[*token as usize]);
        }
        result
    }
}
/// A Token is a unique identifier for a string.
pub type Token = u32;

/// A Parsed file is a list of tokens and a list of start positions for each token.
#[derive(Debug)]
pub struct Parsed {
    /// The tokens in the file
    pub tokens: Vec<Token>,
    /// The start position of each token in the original file plus one additional entry marking the end of the last token
    pub starts: Vec<usize>,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_tokens() {
        let mut tokens = Tokens::new();
        let a = tokens.get_token(b"hello".to_vec());
        let b = tokens.get_token(b"world".to_vec());
        let c = tokens.get_token(b"hello".to_vec());
        assert_eq!(a, c);
        assert_ne!(a, b);
        assert_ne!(b, c);
        assert_eq!(tokens.token_upper_bound(), 2);
        assert_eq!(
            tokens.get_token_images(&vec![a, b, c]),
            vec!["hello", "world", "hello"]
        );
    }

    #[test]
    fn test_printer() {
        let mut tokens = Tokens::new();
        let a = tokens.get_token(b"hello".to_vec());
        let b = tokens.get_token(b"world".to_vec());
        let printer = tokens.printer();
        assert_eq!(printer.print(a), b"hello");
        assert_eq!(printer.print(b), b"world");
        assert_eq!(printer.concat(&vec![a, b]), b"helloworld");
    }
}
