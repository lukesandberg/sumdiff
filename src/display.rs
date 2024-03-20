use std::io::Error;
use std::io::Write;
use std::ops::Range;

use crate::hunk::Hunk;
use crate::hunk::HunkPart;
use crate::token::{Parsed, Token, TokenPrinter};

pub struct DisplayConfig {
    pub insert_prefix: &'static [u8],
    pub delete_prefix: &'static [u8],
    pub common_prefix: &'static [u8],
    pub token_suffix: &'static [u8],
}

pub fn display_diff(
    config: &DisplayConfig,
    printer: &TokenPrinter,
    left: &Parsed,
    right: &Parsed,
    hunks: &[Hunk],
    mut writer: impl Write,
) -> Result<(), Error> {
    for hunk in hunks {
        writer.write_all(hunk.range().as_bytes())?;
        for part in hunk.parts() {
            match part {
                HunkPart::Delete(range) => {
                    flush_tokens(
                        printer,
                        &left.tokens,
                        range,
                        config.delete_prefix,
                        config.token_suffix,
                        &mut writer,
                    )?;
                }
                HunkPart::Insert(range) => {
                    flush_tokens(
                        printer,
                        &right.tokens,
                        range,
                        config.insert_prefix,
                        config.token_suffix,
                        &mut writer,
                    )?;
                }
                HunkPart::Common(range) => {
                    flush_tokens(
                        printer,
                        // could use left or right here, they are the same
                        &left.tokens,
                        &range.left_range(),
                        config.common_prefix,
                        config.token_suffix,
                        &mut writer,
                    )?;
                }
            }
        }
    }
    Ok(())
}

fn flush_tokens(
    printer: &TokenPrinter,
    tokens: &[Token],
    range: &Range<usize>,
    prefix: &[u8],
    suffix: &[u8],
    writer: &mut impl std::io::Write,
) -> Result<(), Error> {
    for token in &tokens[range.clone()] {
        writer.write_all(prefix)?;
        writer.write_all(&printer.print(*token))?;
        writer.write_all(suffix)?;
    }
    Ok(())
}
