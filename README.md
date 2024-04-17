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

### Struct of arrays vs array of structs?

When 'parsing' files we want to track metadata bout the data such as file offset, token length, AST depth, etc.  It would be tempting to construct a `struct Token` to represent this but instead we have chosen a 'struct of arrays' apporach and store metadata in parallel datastructures.

Many algorithm runtimes will be dominated by the cost of scanning through token lists, while querying data about offset and depth is rarer, by moving that data to be stored elsewhere we slightly complicate some output algorithms in exchange for improving memory locality of our diff algorithms.


### Why compare Meyers, Dijkstra and the Kuo and Cross LCS algorithms?

The meyers and KC papers both have this fun detail where they criticize other algorithms, the criticisms are fair but tickled my curiousity. Furthermore, the Kuo and Cross algorithm facinated 
me.

* Kuo and Cross: This paper had some clear opportunities to improve performance (binary search, counting sort, caching offsets, compressing matches), and required some additional work over the paper to actually produce an LCS instead of a length.  By working on these I was able
to drop an `nlogn` term from the complexity, mitigate the `R` term but was unable to escape the `NL` term.  For standard edits (e.g. code diffs), the LCS is proportional to the length of the text so fundamentally we are still implementing a quadratic algorithm, though admittedly many of the constant factors are small.

* Dijkstra: the Meyers paper discusses that improving the performance of dijkstra on edit graphs requires 'complex linked list and bucket queues'.  Indeed the literature on these (Radix Queues, Fibonnaci Queues, Dials Algorithm) is complex.  But edit graphs are very simple with all edge weights as `0` or `1`, thus I was able to use a simple circular buffer (a `VecDeque`) as a priority queue, leading to a trivial `O(1)` solution.  The remaining problem is simply the datastructure management.  The size of our `predecessor` datastructures scale with the number of
nodes visited.  In principle this is `O(ND)` which is much worse than Meyers.
  - To do better I believe would require adopting more techniques from meyers, e.g. an explicit
    representation of the 'frontier diagonals'.


After working on all 3 it is clear that the criticisms of Meyers are completely accurate (except for the thing about 'complex priority queues').  The KC algorithm only becomes competitive when the alphabet size becomes large and the number of edits is also large (aka `D` is large), and for dijkstra to become competitive does require stealing representation techniques from Meyers.

The key insight of Meyers is to target `D` (the number of edits) instead of `L` (the LCS), which is useful for standard diff scenarios where our expecation is that `D` is small.

## TODO

### SIMD

The counting sort implementation requires a prefix-sum computation, this can be accomplished more quickly using SIMD a clear TODO, however since `std::simd` is unstable this is going to be complex. 
