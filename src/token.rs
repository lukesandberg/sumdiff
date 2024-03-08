use std::collections::HashMap;

/// A Tokens is a collection of unique strings to unique token identifiers
#[derive(Debug)]
pub struct Tokens {
    universe: HashMap<Vec<u8>, Token>,
}
impl<'a> Tokens {
    #[must_use]
    pub fn new() -> Self {
        Tokens {
            universe: HashMap::new(),
        }
    }
    /// Get the token for a string, creating a new one if necessary.
    pub fn get_token_ref(&mut self, s: &Vec<u8>) -> Token {
        match self.universe.get(s) {
            Some(token) => *token,
            None => {
                let id = self.universe.len() as Token;
                // Since we only have an immutable reference we need to copy to insert
                self.universe.insert(s.clone(), id);
                id
            }
        }
    }
    /// Get the token for a string, creating a new one if necessary.
    pub fn get_token(&mut self, s: Vec<u8>) -> Token {
        match self.universe.get(&s) {
            Some(token) => *token,
            None => {
                let id = self.universe.len() as Token;
                // Since we only have an immutable reference we need to copy to insert
                self.universe.insert(s, id);
                id
            }
        }
    }

    /// Get the string for a token, used only for testing
    #[cfg(test)]
    pub fn get_token_image(&self, token: Token) -> Option<String> {
        self.universe
            .iter()
            .find(|(_, &id)| id == token)
            .map(|(s, _)| String::from_utf8(s.to_vec()).unwrap())
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
