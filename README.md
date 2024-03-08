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

TODO(luke): get things working first

## Implementation questions and decisions

### Utf8 Validation?

Should we validate our inputs as utf-8?  It seems reasonable however experience has shown that debugging a utf8 rejection error is difficult.  So on the one hand we could just make sure the error messages are useful, but on the other we could declare that this isn't _our_ problem.

For the purposes of a diff tool we should assume that the inputs are useful to someone (they bothered to ask us to compare them), so rejecting them isn't particularly useful.

### Struct of arrays vs array of structs

When 'parsing' files we want to track metadata bout the data such as file offset, token length, AST depth, etc.  It would be tempting to construct a `struct Token` to represent this but instead we have chosen a 'struct of arrays' apporach and store metadata in parallel datastructures.

Many algorithm runtimes will be dominated by the cost of scanning through token lists, while querying data about offset and depth is rarer, by moving that data to be stored elsewhere we slightly complicate some output algorithms in exchange for improving memory locality of our diff algorithms.