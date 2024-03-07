# sumdiff
A tool to summarize large diffs

To aid in large refactorings it is necessary to understand the diff being applied while avoiding the burden of reviewing every line.

## Approach

The high level approach is to identify the crux of the changes and to ignore irrelevant _textual_ details.

1. Identify the changing parts of files (aka run `diff`)
1. Tokenize changing parts to understand lexical context
   - requires a language sensitive lexer
   - and a fallback lexer for ambiguous inputs (or unsupported languages)
1. Run `diff` again
1. Break the edit script into parts
1. Select and print the most popular edits that cover all changes.


## Huh?

What we are attempting to do is identify _logical_ and _prevalent_ changes.  For a set of scripted edits we can expect uninformity which should allow a small set of 'logical' edits to cover all the actual edits.

## Examples

TODO(luke)