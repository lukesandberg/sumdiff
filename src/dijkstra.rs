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
use std::hash::{Hash, Hasher};

#[derive(Debug, Clone, Copy)]
enum Pred {
    CommonBegin(usize),
    CommonEnd(usize),
    None,
}

impl Pred {
    fn for_next(&self, i: usize) -> Pred {
        match self {
            // If node is the end of a common sequence, then we are branching so set our pred as the
            // common end of that sequence
            Pred::CommonBegin(_) => Pred::CommonEnd(i),
            // If node is pointing at the end of a common sequence point there as well.
            Pred::CommonEnd(pe) => Pred::CommonEnd(*pe),
            // If we are pointing at the origin, then there is just some sequence of inserts and deletes
            Pred::None => Pred::None,
        }
    }
}

struct Vertex {
    left: usize,
    right: usize,
    distance: u32,
    pred: Pred,
}

impl Hash for Vertex {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.left.hash(state);
        self.right.hash(state);
    }
}
impl PartialEq for Vertex {
    fn eq(&self, other: &Self) -> bool {
        self.left == other.left && self.right == other.right
    }
}
impl Eq for Vertex {}

struct NodePool {
    pool: Vec<Vertex>,
    index: HashMap<(usize, usize), usize>,
    max_distance: u32,
}

impl NodePool {
    fn new(n: usize, m: usize) -> Self {
        let mut pool = Vec::with_capacity(n + m);
        let mut index = HashMap::with_capacity(n + m);
        pool.push(Vertex {
            left: 0,
            right: 0,
            distance: 0,
            pred: Pred::None,
        });
        index.insert((0, 0), 0);
        let max_distance: u32 = (n + m + 1).try_into().unwrap();
        pool.push(Vertex {
            left: n,
            right: m,
            distance: max_distance,
            pred: Pred::None,
        });
        index.insert((n, m), 1);
        NodePool {
            pool,
            index,
            max_distance,
        }
    }

    fn dest(&self) -> &Vertex {
        self.get_by_index(1)
    }

    fn get_by_index(&self, index: usize) -> &Vertex {
        &self.pool[index]
    }

    fn get_or_insert(&mut self, vertex: (usize, usize)) -> (&mut Vertex, usize) {
        match self.index.entry(vertex) {
            Entry::Occupied(mut entry) => {
                let index = entry.get_mut();
                let vertex = &mut self.pool[*index];
                (vertex, *index)
            }
            Entry::Vacant(entry) => {
                let i = self.pool.len();
                self.pool.push(Vertex {
                    left: vertex.0,
                    right: vertex.1,
                    distance: self.max_distance,
                    pred: Pred::None,
                });
                entry.insert(i);
                (&mut self.pool[i], i)
            }
        }
    }
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
fn neighbors(start: &Vertex, left: &[Token], right: &[Token]) -> Vec<((usize, usize), bool)> {
    let n = left.len();
    let m = right.len();
    let start_l = start.left;
    let start_r = start.right;
    if start_l == n {
        debug_assert!(
            start_r != m,
            "should never query neighbors of the end vertex"
        );
        vec![((n, start_r + 1), false)]
    } else if start_r == m {
        vec![((start_l + 1, m), false)]
    } else {
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
            vec![((l, r), true)]
        } else {
            // only explore non diagonal if we can't advance on the diagonal
            // Visit deletions and then insertions, for equivalent cost paths this will prioritize deletions
            vec![
                ((start_l + 1, start_r), false),
                ((start_l, start_r + 1), false),
            ]
        }
    }
}

fn recover_path<'a>(
    mut vertex: &'a Vertex,
    matches: &mut Vec<CommonRange>,
    prefix: usize,
    pool: &'a NodePool,
) {
    vertex = match vertex.pred {
        Pred::None => {
            return; // There were no matches, just stop
        }
        Pred::CommonBegin(_) => unreachable!(), // The end node ends on a common sequence, this is impossible due to our suffix removal
        Pred::CommonEnd(i) => pool.get_by_index(i),
    };
    while let Pred::CommonBegin(b) = vertex.pred {
        let pred = pool.get_by_index(b);
        let left_delta = vertex.left - pred.left;
        let right_delta = vertex.right - pred.right;
        debug_assert!(left_delta > 0 && right_delta == left_delta);
        // This was a diagonal move, produce a range.
        matches.push(CommonRange {
            left_start: prefix + pred.left,
            right_start: prefix + pred.right,
            length: left_delta,
        });

        vertex = match pred.pred {
            Pred::None => {
                break; // There were no more matches, just stop
            }
            Pred::CommonBegin(_) => unreachable!(), // A begin can point at a begin
            Pred::CommonEnd(i) => pool.get_by_index(i),
        };
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
    let mut pool = NodePool::new(n, m);
    // This forms a simple 2 bucket bucket queue.
    // All the edges in our graph have a weight of 0 or 1, and dijkstra's algorithm is
    // monotonic meaning that having explored an NodeMeta at distance `d` we will never explore
    // a NodeMeta at distance `f < d `.  So we can just put discovered NodeMetas at distance zero
    // at the front and NodeMetas at distance 1 at the back.  This does change the order of
    // traversal, but it is still correct and gives us trivial O(1) insertions and deletions.
    // This is known as a bucket queue.
    let mut h: VecDeque<usize> = VecDeque::with_capacity(n + m);
    h.push_front(0);

    // Tracks the offset where the current nodes at relative distance zero end which
    // tells us when to increment `front_distance`.
    let mut num_front = 1;
    let mut front_distance: u32 = 0;

    let mut explored_points: usize = 0;
    while let Some(i) = h.pop_front() {
        if num_front == 0 {
            // We have explored all the NodeMetas in the current frontier
            // We can now move to the next frontier
            num_front = h.len();
            front_distance += 1;
        } else {
            num_front -= 1;
        }
        debug_assert!(
            front_distance as usize <= n + m,
            "front_distance = {}, n = {}, m = {}",
            front_distance,
            n,
            m
        );
        let node = pool.get_by_index(i);
        let nd = node.distance;
        if nd < front_distance {
            // We already found a shorter path to this NodeMeta.
            continue;
        }
        let nl = node.left;
        let nr = node.right;
        explored_points += 1;
        if nl == n && nr == m {
            break;
        }
        let np = node.pred;
        for (nv, is_match) in neighbors(node, left, right) {
            let neighbor_distance = nd + if is_match { 0 } else { 1 };
            let (neighbor, index) = pool.get_or_insert(nv);
            if neighbor.distance > neighbor_distance {
                neighbor.distance = neighbor_distance;
                match is_match {
                    true => {
                        neighbor.pred = Pred::CommonBegin(i);
                        h.push_front(index);
                        num_front += 1;
                    }
                    false => {
                        neighbor.pred = np.for_next(i);
                        h.push_back(index);
                    }
                }
            }
        }
    }

    let end = pool.dest();
    if end.distance as usize > n + m {
        panic!(
            "Should have found a path to the end. left={:?} right={:?}",
            left, right
        );
    };
    debug_assert!(explored_points <= (n + m) * end.distance as usize);
    recover_path(end, &mut matches, prefix, &pool);
    if suffix > 0 {
        matches.push(CommonRange {
            left_start: prefix + n,
            right_start: prefix + m,
            length: suffix,
        });
    }
    matches
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
        assert_eq!(lcs, vec![(1, 1), (3, 3)]);
    }

    #[test]
    fn test_diagonal_crash() {
        let mut tokens = Tokens::new();
        let left = lex_characters(&mut tokens, "abcbd");
        let right = lex_characters(&mut tokens, "dabcddb");
        let lcs = CommonRange::flatten(&dijkstra(&left, &right));
        check_is_lcs(&lcs, &left, &right).unwrap();
        assert_eq!(lcs, vec![(0, 1), (1, 2), (2, 3), (4, 4)]);
    }
}
