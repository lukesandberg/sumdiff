/*!
An implementation of dijkstra's algorithm to find the shortest path through the edit graph.

After implementing it is pretty clear that the algorithm is not the best choice for this problem, but
it was a fun exercise to implement it.  While implementing it I had the following realizations:
1. The edit graph is a DAG and dijkstra is designed for general graphs.
2. For DAGs we can find the shortest path by traversing in topological order which is trivial through the
   edit graph (just traverse diagnonals). Then we can use dynamic programming to find the shortest path.
   The problem is that this is trivially an O(NM) algorithm, and this is more or less the insight of
   the meyers algorithm.
3. As implemented storage requirements are way too high.  This could be improved, but the improvements
   are tantamount to implementing meyers (since we only need to strore the previous diagonal).

The general problem with graph algorithms is that they scale with the size of the graph, and the edit graph
is O(NM) in size.  We need to take advantage of inherent structure in the edit graph to do better.

To do that in this case would involve things like:
 - aggressively following 'snakes' (diagonals) in the graph where the two substrings are aligned.
 - improving boundary conditions to prune unproductive paths (this is more or less A* search)
 - adopt a fibonnaci heap to decrease heap size and improve performance
 - use a more efficient encoding of the predecessor path (e.g. we only really need to store non matches)

The above are all more or less the insights of the meyers algorithm. Still dijkstra is fundamentally easier
to reason about (dijkstra wrote in on a napkin at a cafe while daydreaming about bridges!). The meyers algorithm
is the standard approach however for very good reasons.

Still finding a way to crank this to the max would be fun if only because Meyers said:
"""
While one could argue that further refinement [of dijkstra] leads to the simple algorithm of this paper,
the connection becomes so tenuous that the direct and easily motivated derivation used in this section
is preferable."
"""
*/

use crate::lcs_utils::{remove_suffixes_and_prefixes, Trimmed};
use crate::token::{CommonRange, Token};
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::collections::VecDeque;

type Vertex = (usize, usize);

#[derive(Debug, Clone, Copy)]
struct NodeMeta {
    pred: Vertex,
    distance: u32,
}

/// Returns the neighbors of this NodeMeta in the edit graph.
///
/// The neighbors of a NodeMeta (x,y) is either of sequences a and b of length n and m respectively
///
/// - if `x`` == `n-1`` then `(`x, m-1)`` with cost `m-y` meaning we insert the rest of b
/// - if `y`` == `m-1`` then `(`n-1, y)`` with cost `n-x` meaning we delete the rest of a
/// - else if a[x] == b[y] then `neighbors((`x+1, y+1))` with cost 0 meaning we match
/// - else `(`x+1, y)`` with cost 1 and `(`x, y+1)`` with cost 1 meaning we either delete one token from a or or insert one token into b
///
/// The special case is for the match, to reduce some overhead we follow the full diagonal.
fn neighbors((start_l, start_r): Vertex, left: &[Token], right: &[Token]) -> Vec<(Vertex, u32)> {
    let n = left.len();
    let m = right.len();
    if start_l == n || start_r == m {
        debug_assert!(
            start_l != n || start_r != m,
            "should never query neighbors of the end vertex"
        );
        // We are at the end of one of the sequences, just jump to the end
        // This is the one place that the distance is >1.  This might require a special case in the caller
        return vec![((n, m), std::cmp::max(n - start_l, m - start_r) as u32)];
    }
    // follow the whole diagonal
    // TODO: this operation can be optimized to O(1) using two suffix arrays
    let mut l = start_l;
    let mut r = start_r;
    if left[l] == right[r] {
        l += 1;
        r += 1;
        while l < n && r < m && left[l] == right[r] {
            l += 1;
            r += 1;
        }
        return vec![((l, r), 0)];
    }
    // only explore non diagonal if we can't advance on the diagonal
    // Visit deletions and then insertions, for equivalent cost paths this will prioritize deletions
    vec![((start_l + 1, start_r), 1), ((start_l, start_r + 1), 1)]
}

fn recover_path<'a>(
    mut vertex: Vertex,
    mut node: &'a NodeMeta,
    matches: &mut Vec<CommonRange>,
    prefix: usize,
    pool: &'a HashMap<Vertex, NodeMeta>,
) {
    while vertex != (0, 0) {
        let pred_vertex = node.pred;
        let left_delta = vertex.0 - pred_vertex.0;
        let right_delta = vertex.1 - pred_vertex.1;
        if left_delta == right_delta {
            debug_assert!(left_delta == right_delta);
            // This was a diagonal move
            matches.push(CommonRange {
                left_start: prefix + pred_vertex.0,
                right_start: prefix + pred_vertex.1,
                length: left_delta,
            });
        } else {
            // This was a horizontal or vertical move
            // TODO: find a way to avoid creating these NodeMetas
        }
        vertex = pred_vertex;
        node = &pool.get(&pred_vertex).unwrap();
    }
    if prefix > 0 {
        matches[1..].reverse();
    } else {
        matches.reverse();
    }
}

/// An implementation of dijkstra's algorithm to find the shortest path through the edit
/// graph between left and right.
///
/// The edit graph is an implicit graph where every vertex is a pair of indexes into the left and right slices
/// and edges only go 'forward', e.g. (x,y) has edges to (x+1,y) which represents a delete, (x,y+1) which
/// represents an insert, and maybe (x+1,y+1) if the tokens match.
///
/// The edge weights are 0 for match edges (diagonals) and 1 for inserts or deletes
pub fn dijkstra(left: &[Token], right: &[Token]) -> Vec<CommonRange> {
    let Trimmed {
        left,
        right,
        mut matches,
        prefix,
        suffix,
    } = match remove_suffixes_and_prefixes(left, right) {
        Ok(matches) => {
            return matches;
        }
        Err(t) => t,
    };

    let n = left.len();
    let m = right.len();
    debug_assert!((n > 1 && m > 0) || (n > 0 && m > 1), "n = {}, m = {}", n, m);
    // Limiting the size of this map is the most important optimization we can make.
    // Ideas:
    //   - We could store the keys as a single usize by packing both indices together. This
    //     is only possible if we require that n and m are less than 2^32.  This would be a
    //     constant size optimization
    //
    let mut pool: HashMap<Vertex, NodeMeta> = HashMap::with_capacity(n + m);

    let origin = (0, 0);
    pool.insert(
        origin,
        NodeMeta {
            distance: 0,
            pred: origin, // self link
        },
    );

    // This forms a simple 2 bucket bucket queue.
    // All the edges in our graph have a weight of 0 or 1, and dijkstra's algorithm is
    // monotonic meaning that having explored an NodeMeta at distance `d` we will never explore
    // a NodeMeta at distance `f < d `.  So we can just put discovered NodeMetas at distance zero
    // at the front and NodeMetas at distance 1 at the back.  This does change the order of
    // traversal, but it is still correct and gives us trivial O(1) insertions and deletions.
    // This is known as a bucket queue.
    let mut h: VecDeque<Vertex> = VecDeque::with_capacity(n + m);
    h.push_front(origin);
    // Tracks the offset where the current nodes at relative distance zero end which
    // tells us when to increment `front_distance`.
    let mut num_front = 1;
    let mut front_distance = 0;

    let end_vertex = (n, m);
    while !h.is_empty() {
        if num_front == 0 {
            // We have explored all the NodeMetas in the current frontier
            // We can now move to the next frontier
            num_front = h.len();
            front_distance += 1;
        }
        let vertex = h.pop_front().unwrap();

        if vertex == end_vertex {
            break;
        }
        let node = pool.get(&vertex).unwrap();
        num_front -= 1;
        let nd = node.distance;
        if nd < front_distance {
            // We found a shorter path to this NodeMeta.
            continue;
        }

        for (nv, edge_weight) in neighbors(vertex, left, right) {
            let neighbor_distance = nd + edge_weight;
            let updated = match pool.entry(nv) {
                Entry::Occupied(mut entry) => {
                    let cur = entry.get_mut();
                    if cur.distance > neighbor_distance {
                        cur.distance = neighbor_distance;
                        cur.pred = vertex;
                        // If we had an efficient mechanism it would be good to delete `nv` from the upcoming heap.
                        // This would require tracking offsets in the NodeMeta structure and then deleting an entry
                        // in our heap would be O(N).  Instead we detect the stale NodeMeta above when we compare the
                        // NodeMeta distance to the `front` distance.
                        true
                    } else {
                        false
                    }
                }
                Entry::Vacant(entry) => {
                    entry.insert(NodeMeta {
                        distance: neighbor_distance,
                        pred: vertex,
                    });
                    true
                }
            };
            if updated {
                match edge_weight {
                    0 => {
                        h.push_front(nv);
                        num_front += 1;
                    }
                    1 => h.push_back(nv),
                    _ => {
                        // This happens on a 'jump' to the end.
                        debug_assert!((n, m) == nv);
                        // In this case we don't update the heap since there are no neighbors to explore
                        // and we know it is larger than everything else.  The pool is already updated with the
                        // distance and predecessor.
                    }
                }
            }
        }
    }

    let end = match pool.get(&end_vertex) {
        Some(end) => end,
        None => {
            panic!(
                "Should have found a path to the end. left={:?} right={:?}",
                left, right
            );
        }
    };
    recover_path(end_vertex, &end, &mut matches, prefix, &pool);
    if suffix > 0 {
        matches.push(CommonRange {
            left_start: prefix + n,
            right_start: prefix + m,
            length: suffix,
        });
    }
    return matches;
}

#[cfg(test)]
mod tests {
    use crate::{lcs_utils::check_is_lcs, lex::lex_characters, token::Tokens};

    use super::*;

    // A reduced regression test for a bug where the algorithm would crash
    #[test]
    fn test_dijkstra_crash_trivial() {
        let mut tokens = Tokens::new();
        let left = lex_characters(&mut tokens, "ab");
        let right = lex_characters(&mut tokens, "cd");
        let lcs = CommonRange::flatten(&dijkstra(&left, &right));
        assert_eq!(lcs, vec![]);
    }

    #[test]
    fn test_dijkstra_crash() {
        let mut tokens = Tokens::new();
        let left = lex_characters(&mut tokens, "iic");
        let right = lex_characters(&mut tokens, "aicb");
        let lcs = CommonRange::flatten(&dijkstra(&left, &right));
        check_is_lcs(&lcs, &left, &right).unwrap();
        assert_eq!(lcs, vec![(1, 1), (2, 2)]);
    }

    #[test]
    fn test_dijkstra_crash_2() {
        let mut tokens = Tokens::new();
        let left = [15, 30, 14, 9, 22]
            .iter()
            .map(|s| tokens.get_token(s.to_string().as_bytes().to_vec()))
            .collect::<Vec<Token>>();
        let right = [29, 30, 26, 9, 14]
            .iter()
            .map(|s| tokens.get_token(s.to_string().as_bytes().to_vec()))
            .collect::<Vec<Token>>();
        let lcs = CommonRange::flatten(&dijkstra(&left, &right));
        check_is_lcs(&lcs, &left, &right).unwrap();
        assert_eq!(lcs, vec![(1, 1), (2, 4)]);
    }
}
