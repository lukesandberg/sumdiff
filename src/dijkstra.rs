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

use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::{BinaryHeap, HashMap};

use crate::lcs_utils::{remove_suffixes_and_prefixes, Trimmed};
use crate::token::{CommonRange, Token};
use std::rc::Rc;

#[derive(Debug, Clone)]
struct Node {
    left_index: usize,
    right_index: usize,
    distance: u32,
    pred: Option<Rc<Node>>,
}

impl Node {
    /// Returns the neighbors of this node in the edit graph.
    ///
    /// The neighbors of a node (x,y) is either of sequences a and b of length n and m respectively
    ///
    /// - if `x`` == `n-1`` then `(`x, m-1)`` with cost `m-y` meaning we insert the rest of b
    /// - if `y`` == `m-1`` then `(`n-1, y)`` with cost `n-x` meaning we delete the rest of a
    /// - else if a[x] == b[y] then `neighbors((`x+1, y+1))` with cost 0 meaning we match
    /// - else `(`x+1, y)`` with cost 1 and `(`x, y+1)`` with cost 1 meaning we either delete one token from a or or insert one token into b
    ///
    /// The special case is for the match, to reduce some overhead we follow the full diagonal.
    fn neighbors(&self, left: &[Token], right: &[Token]) -> Vec<(usize, usize, u32)> {
        let start_l = self.left_index;
        let start_r = self.right_index;
        let n = left.len();
        let m = right.len();
        if start_l == n || start_r == m {
            // We are at the end of one of the sequences, just jump to the end
            // This is the one place that the distance is >1.  This might require a special case in the caller
            return vec![(n, m, std::cmp::max(n - start_l, m - start_r) as u32)];
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
            return vec![(l, r, 0)];
        }
        // only explore non diagonal if we can't advance on the diagonal
        vec![(start_l + 1, start_r, 1), (start_l, start_r + 1, 1)]
    }
}

impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.distance.eq(&other.distance)
    }
}

impl Eq for Node {}

impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Node {
    /// Nodes are ordered by decreasing distance
    fn cmp(&self, other: &Self) -> Ordering {
        other.distance.cmp(&self.distance)
    }
}

fn recover_path(mut node: &Node, matches: &mut Vec<CommonRange>, prefix: usize) {
    while let Some(ref pred) = &node.pred {
        let left_delta = node.left_index - pred.left_index;
        let right_delta = node.right_index - pred.right_index;
        if left_delta == right_delta {
            debug_assert!(left_delta == right_delta);
            // This was a diagonal move
            matches.push(CommonRange {
                left_start: prefix + pred.left_index,
                right_start: prefix + pred.right_index,
                length: left_delta,
            });
        } else {
            // This was a horizontal or vertical move
            // TODO: find a way to avoid creating these nodes
        }
        node = pred;
    }
    // We should have made it back to the beginning.
    debug_assert!(node.left_index == 0 && node.right_index == 0);
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
    let mut distances: HashMap<(usize, usize), u32> = HashMap::new();
    // Because our distances are integers comparing them is incredibly fast
    // So we should use a d-ary heap where d is probably at least 16.
    // Additionally this would allow us to not store the distance in the node.
    // Or use a bucket queue, which allows for O(1) insertions and deletions.
    // In the case of a bucket queue I believe we would only need two buckets
    // So perhaps we could just use a vecdeque and an integer to track the transition?
    // neighbors with distance zero would go to the front and neighbors with distance one would go to the back.
    // Or we could use two hashsets (or btreesets) to allow for O(1) deletions.
    let mut heap: BinaryHeap<Rc<Node>> = BinaryHeap::with_capacity(std::cmp::max(n, m));

    distances.insert((0, 0), 0);
    heap.push(Rc::new(Node {
        left_index: 0,
        right_index: 0,
        distance: 0,
        pred: None,
    }));
    while let Some(node) = heap.pop() {
        let sl = node.left_index;
        let sr = node.right_index;
        if distances.get(&(sl, sr)) != Some(&node.distance) {
            // This means we found a shorter path to this node already.
            continue;
        }
        if sl == n && sr == m {
            recover_path(&node, &mut matches, prefix);
            if suffix > 0 {
                matches.push(CommonRange {
                    left_start: prefix + n,
                    right_start: prefix + m,
                    length: suffix,
                });
            }
            return matches;
        }
        for (nleft, nright, edge_weight) in node.neighbors(left, right) {
            let neighbor_distance = node.distance + edge_weight;
            let updated = match distances.entry((nleft, nright)) {
                Entry::Occupied(entry) => {
                    let cur = entry.into_mut();
                    if *cur > neighbor_distance {
                        *cur = neighbor_distance;
                        // TODO: delete the old node from the heap
                        true
                    } else {
                        false
                    }
                }
                Entry::Vacant(entry) => {
                    entry.insert(neighbor_distance);
                    true
                }
            };
            if updated {
                // This node may already be in the heap, but with a larger distance.
                // our if-statement above will ensure it is ignored in that case.
                heap.push(Rc::new(Node {
                    left_index: nleft,
                    right_index: nright,
                    distance: neighbor_distance,
                    pred: Some(node.clone()),
                }));
            }
        }
    }
    panic!("Should have found a path to the end");
}
