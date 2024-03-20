/// This module contains the logic for converting matches into hunks
/// A hunk is a section of a diff that contains a common context and a set of insertions and deletions
use std::ops::Range;

use crate::token::CommonRange;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HunkPart {
    Common(CommonRange),
    Insert(Range<usize>),
    Delete(Range<usize>),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Hunk {
    // Storing a boxed slice would be slighlyt superior but we don't expect to have many hunks or parts
    // so the overhead of the extra allocation is not worth optimizing
    parts: Vec<HunkPart>,
}

impl Hunk {
    fn new(leading: CommonRange) -> Hunk {
        let mut hunk = Hunk {
            // Use a small default capacity since it is typically 'common', 'insert|delete' 'common'
            // Sometimes it will be smaller or longer but this is typical.
            parts: Vec::with_capacity(4),
        };
        hunk.parts.push(HunkPart::Common(leading));
        hunk
    }
    fn add_part(self: &mut Hunk, part: &HunkPart) {
        let (left, right) = self.hunk_ranges();
        match part {
            HunkPart::Common(c) => {
                debug_assert!(matches!(
                    self.parts.last(),
                    Some(HunkPart::Insert(_)) | Some(HunkPart::Delete(_))
                ));
                debug_assert_eq!(c.left_start, left.end);
                debug_assert_eq!(c.right_start, right.end);
            }
            HunkPart::Insert(r) => {
                debug_assert_eq!(r.start, right.end);
            }
            HunkPart::Delete(r) => {
                debug_assert_eq!(r.start, left.end);
            }
        }
        self.parts.push(part.clone());
    }

    /// Returns the ranges of the hunk
    fn hunk_ranges(self: &Hunk) -> (Range<usize>, Range<usize>) {
        let first = match &self.parts[0] {
            HunkPart::Common(c) => c,
            // We should never have a hunk that doesn't start with a common part
            _ => unreachable!(),
        };

        let left_start = first.left_start;
        let right_start = first.right_start;
        let mut left_length = first.length;
        let mut right_length = first.length;
        for part in &self.parts[1..] {
            match part {
                HunkPart::Common(CommonRange { length, .. }) => {
                    left_length += length;
                    right_length += length;
                }
                HunkPart::Insert(r) => {
                    right_length += r.len();
                }
                HunkPart::Delete(r) => {
                    left_length += r.len();
                }
            }
        }
        (
            Range {
                start: left_start,
                end: left_start + left_length,
            },
            Range {
                start: right_start,
                end: right_start + right_length,
            },
        )
    }

    /// Formats a unified diff range marker with a trailing newline
    pub fn range(self: &Hunk) -> String {
        let (left, right) = self.hunk_ranges();
        format!(
            "@@ -{},{} +{},{} @@\n",
            // By convention we use 1-based line numbers
            left.start + 1,
            left.len(),
            right.start + 1,
            right.len()
        )
    }

    pub fn parts(self: &Hunk) -> &[HunkPart] {
        &self.parts
    }
}

// TODO: implementing this as two generators (or custom iterators) would be
// interesting as it would allow us to interleave construction with printing.
// On the other hand these functions operate as a 'compression' on the matches
// so all returned vectors should be small. However it might be interesting to
// enable some kind of parallelisation in other parts of the tool.

fn matches_to_hunk_parts(
    matches: &[CommonRange],
    left_length: usize,
    right_length: usize,
) -> Vec<HunkPart> {
    let mut hunk_parts = Vec::new();
    // pretend there is a match at the very beginning of the files
    // This makes it simpler to understand inserts/deletes at the beginning
    let mut prev_left_end = 0;
    let mut prev_right_end = 0;

    // Append a dummy match to the end to detect inserts/deletes at the end
    for range in matches.iter().chain(std::iter::once(&CommonRange {
        left_start: left_length,
        right_start: right_length,
        length: 0,
    })) {
        if prev_left_end < range.left_start {
            hunk_parts.push(HunkPart::Delete(Range {
                start: prev_left_end,
                end: range.left_start,
            }));
        }
        if prev_right_end < range.right_start {
            hunk_parts.push(HunkPart::Insert(Range {
                start: prev_right_end,
                end: range.right_start,
            }));
        }
        if range.length > 0 {
            hunk_parts.push(HunkPart::Common(*range));
        }
        prev_left_end = range.left_start + range.length;
        prev_right_end = range.right_start + range.length;
    }

    hunk_parts
}

pub fn matches_to_hunks(
    max_context_tokens: usize,
    matches: &[CommonRange],
    left_length: usize,
    right_length: usize,
) -> Vec<Hunk> {
    let raw_parts = matches_to_hunk_parts(matches, left_length, right_length);
    let mut hunks: Vec<Hunk> = Vec::new();
    for (index, hunk) in raw_parts.iter().enumerate() {
        match hunk {
            HunkPart::Common(CommonRange {
                left_start,
                right_start,
                length,
            }) => {
                let more = index < raw_parts.len() - 1;
                if let Some(last_hunk) = hunks.last_mut() {
                    // If the common part added to the previous hunk would overlap with the next hunk
                    // we instead merge the hunks even thought the common context is longer.
                    if *length < 2 * max_context_tokens && more {
                        last_hunk.add_part(hunk);
                        continue; // We want to keep extending
                    } else {
                        // Common is either too long, or we are at the end
                        // divy it up and take the first max_context_tokens as a trailing common part
                        last_hunk.add_part(&HunkPart::Common(CommonRange {
                            left_start: *left_start,
                            right_start: *right_start,
                            length: std::cmp::min(max_context_tokens, *length),
                        }));
                        // fall through to maybe start a new hunk
                    }
                }
                // If there is an upcoming part then begin a new hunk
                if more {
                    // We only want to grab up to max_context_tokens tokens for the common prefix
                    let new_length = std::cmp::min(*length, max_context_tokens);
                    let delta = length - new_length;
                    hunks.push(Hunk::new(CommonRange {
                        left_start: *left_start + delta,
                        right_start: *right_start + delta,
                        length: new_length,
                    }));
                }
            }
            _ => {
                // we are either extending the previous hunk or starting a new one
                // By construction we know that any non-common part will be added onto the hunk since
                // we only split common ranges
                if let Some(last_hunk) = hunks.last_mut() {
                    last_hunk.add_part(hunk);
                } else {
                    // Otherwise we are starting a new hunk at the very beginning
                    let mut new_hunk = Hunk::new(CommonRange {
                        left_start: 0,
                        right_start: 0,
                        length: 0,
                    });
                    new_hunk.add_part(hunk);
                    hunks.push(new_hunk);
                }
            }
        }
    }
    hunks
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_matches_to_hunk_parts_empty() {
        // Empty
        assert_eq!(matches_to_hunk_parts(&[], 0, 0), vec![]);
    }
    #[test]
    fn test_matches_to_hunk_parts_change_beginning() {
        // Insert at the begnning
        assert_eq!(
            matches_to_hunk_parts(
                &[CommonRange {
                    left_start: 0,
                    right_start: 1,
                    length: 2
                }],
                2,
                3
            ),
            vec![
                HunkPart::Insert(0..1),
                HunkPart::Common(CommonRange {
                    left_start: 0,
                    right_start: 1,
                    length: 2
                })
            ]
        );

        // Delete at the begnning
        assert_eq!(
            matches_to_hunk_parts(
                &[CommonRange {
                    left_start: 1,
                    right_start: 0,
                    length: 2
                }],
                3,
                2
            ),
            vec![
                HunkPart::Delete(0..1),
                HunkPart::Common(CommonRange {
                    left_start: 1,
                    right_start: 0,
                    length: 2
                })
            ]
        );

        // Replace at the begnning
        assert_eq!(
            matches_to_hunk_parts(
                &[CommonRange {
                    left_start: 1,
                    right_start: 1,
                    length: 2
                }],
                3,
                3
            ),
            vec![
                HunkPart::Delete(0..1),
                HunkPart::Insert(0..1),
                HunkPart::Common(CommonRange {
                    left_start: 1,
                    right_start: 1,
                    length: 2
                })
            ]
        );
    }

    #[test]
    fn test_matches_to_hunk_parts_change_end() {
        // Insert at the end
        assert_eq!(
            matches_to_hunk_parts(
                &[CommonRange {
                    left_start: 0,
                    right_start: 0,
                    length: 2
                }],
                2,
                3
            ),
            vec![
                HunkPart::Common(CommonRange {
                    left_start: 0,
                    right_start: 0,
                    length: 2
                }),
                HunkPart::Insert(2..3),
            ]
        );

        // Delete at the begnning
        assert_eq!(
            matches_to_hunk_parts(
                &[CommonRange {
                    left_start: 0,
                    right_start: 0,
                    length: 2
                }],
                3,
                2
            ),
            vec![
                HunkPart::Common(CommonRange {
                    left_start: 0,
                    right_start: 0,
                    length: 2
                }),
                HunkPart::Delete(2..3),
            ]
        );

        // Replace at the begnning
        assert_eq!(
            matches_to_hunk_parts(
                &[CommonRange {
                    left_start: 0,
                    right_start: 0,
                    length: 2
                }],
                3,
                3
            ),
            vec![
                HunkPart::Common(CommonRange {
                    left_start: 0,
                    right_start: 0,
                    length: 2
                }),
                HunkPart::Delete(2..3),
                HunkPart::Insert(2..3),
            ]
        );
    }

    #[test]
    fn test_matches_to_hunk_parts_change_beginning_and_end() {
        // Replace at the beginning and end
        assert_eq!(
            matches_to_hunk_parts(
                &[CommonRange {
                    left_start: 1,
                    right_start: 1,
                    length: 2
                }],
                4,
                4
            ),
            vec![
                HunkPart::Delete(0..1),
                HunkPart::Insert(0..1),
                HunkPart::Common(CommonRange {
                    left_start: 1,
                    right_start: 1,
                    length: 2
                }),
                HunkPart::Delete(3..4),
                HunkPart::Insert(3..4),
            ]
        );
    }

    #[test]
    fn test_matches_to_hunks_change_beginning_and_end() {
        assert_eq!(
            matches_to_hunks(
                3,
                &[CommonRange {
                    left_start: 1,
                    right_start: 1,
                    length: 2
                }],
                4,
                4
            ),
            vec![Hunk {
                parts: vec![
                    HunkPart::Common(CommonRange {
                        left_start: 0,
                        right_start: 0,
                        length: 0
                    }),
                    HunkPart::Delete(0..1),
                    HunkPart::Insert(0..1),
                    HunkPart::Common(CommonRange {
                        left_start: 1,
                        right_start: 1,
                        length: 2
                    }),
                    HunkPart::Delete(3..4),
                    HunkPart::Insert(3..4),
                ]
            }]
        );
    }

    #[test]
    fn test_matches_to_hunks_2_hunks_change_beginning_and_end() {
        assert_eq!(
            matches_to_hunks(
                1,
                &[CommonRange {
                    left_start: 1,
                    right_start: 1,
                    length: 2
                }],
                4,
                4
            ),
            vec![
                Hunk {
                    parts: vec![
                        HunkPart::Common(CommonRange {
                            left_start: 0,
                            right_start: 0,
                            length: 0
                        }),
                        HunkPart::Delete(0..1),
                        HunkPart::Insert(0..1),
                        HunkPart::Common(CommonRange {
                            left_start: 1,
                            right_start: 1,
                            length: 1
                        })
                    ]
                },
                Hunk {
                    parts: vec![
                        HunkPart::Common(CommonRange {
                            left_start: 2,
                            right_start: 2,
                            length: 1
                        }),
                        HunkPart::Delete(3..4),
                        HunkPart::Insert(3..4),
                    ]
                },
            ]
        );
    }

    #[test]
    fn test_matches_to_hunks_change_middle() {
        assert_eq!(
            matches_to_hunks(
                3,
                &[
                    CommonRange {
                        left_start: 0,
                        right_start: 0,
                        length: 1
                    },
                    CommonRange {
                        left_start: 2,
                        right_start: 2,
                        length: 1
                    }
                ],
                3,
                3
            ),
            vec![Hunk {
                parts: vec![
                    HunkPart::Common(CommonRange {
                        left_start: 0,
                        right_start: 0,
                        length: 1
                    }),
                    HunkPart::Delete(1..2),
                    HunkPart::Insert(1..2),
                    HunkPart::Common(CommonRange {
                        left_start: 2,
                        right_start: 2,
                        length: 1
                    })
                ]
            }]
        );
    }
}
