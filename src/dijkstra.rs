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

use crate::lcs_utils::{common_prefix, common_suffix};
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
    fn neighbors(self: &Self, left: &[Token], right: &[Token]) -> Vec<(usize, usize, u32)> {
        let can_advance_left = self.left_index + 1 < left.len();
        let can_advance_right = self.right_index + 1 < left.len();
        if can_advance_left || can_advance_right {
            if can_advance_left && can_advance_right {
                if left[self.left_index] == right[self.right_index] {
                    return vec![
                        (self.left_index + 1, self.right_index + 1, 0),
                        (self.left_index + 1, self.right_index, 1),
                        (self.left_index, self.right_index + 1, 1),
                    ];
                }
                return vec![
                    (self.left_index + 1, self.right_index, 1),
                    (self.left_index, self.right_index + 1, 1),
                ];
            } else {
                if can_advance_left {
                    return vec![(self.left_index + 1, self.right_index, 1)];
                } else {
                    return vec![(self.left_index, self.right_index + 1, 1)];
                }
            }
        }
        vec![]
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

/// An implementation of dijkstra's algorithm to find the shortest path through the edit
/// graph between left and right.
///
/// The edit graph is an implicit graph where every vertex is a pair of indexes into the left and right slices
/// and edges only go 'forward', e.g. (x,y) has edges to (x+1,y) which represents a delete, (x,y+1) which
/// represents an insert, and maybe (x+1,y+1) if the tokens match.
///
/// The edge weights are 0 for match edges (diagonals) and 1 for inserts or deletes
pub fn dijkstra(left: &[Token], right: &[Token]) -> Vec<CommonRange> {
    let mut n = left.len();
    let mut m = right.len();
    if n == 0 || m == 0 {
        return vec![];
    }
    let prefix = common_prefix(left, right);
    let mut matches = Vec::new();
    if prefix > 0 {
        matches.push(CommonRange {
            left_start: 0,
            right_start: 0,
            length: prefix,
        });
        // Exact match or one sequence is a prefix of the other
        if prefix == n || prefix == m {
            return matches;
        }
    }
    // wait to append the suffix matches until after the main loop
    let suffix = common_suffix(left, right);
    n -= suffix;
    m -= suffix;
    let left = &left[prefix..n];
    let right = &right[prefix..m];
    n -= prefix;
    m -= prefix;
    let mut distances: HashMap<(usize, usize), u32> = HashMap::new();
    let mut heap: BinaryHeap<Rc<Node>> = BinaryHeap::with_capacity(std::cmp::max(n, m));

    distances.insert((0, 0), 0);
    let source = Rc::new(Node {
        left_index: 0,
        right_index: 0,
        distance: 0,
        pred: None,
    });
    heap.push(source.clone());
    let mut dest = None;
    while let Some(node) = heap.pop() {
        if distances.get(&(node.left_index, node.right_index)) != Some(&node.distance) {
            // This means we found a shorter path to this node already.
            continue;
        }
        if node.left_index == n - 1 && node.right_index == m - 1 {
            dest = Some(node);
            break;
        }
        for (left, right, edge_weight) in node.neighbors(left, right) {
            let neighbor_distance = node.distance + edge_weight;
            let updated = match distances.entry((left, right)) {
                Entry::Occupied(entry) => {
                    let cur = entry.into_mut();
                    if *cur > neighbor_distance {
                        *cur = neighbor_distance;
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
                    left_index: left,
                    right_index: right,
                    distance: neighbor_distance,
                    pred: Some(node.clone()),
                }));
            }
        }
    }
    drop(heap);
    drop(distances);
    let mut node = &dest.unwrap();
    let mut c: CommonRange = CommonRange {
        left_start: prefix + node.left_index,
        right_start: prefix + node.right_index,
        length: if left[node.left_index] == right[node.right_index] {
            1
        } else {
            0
        },
    };
    loop {
        match &node.pred {
            Some(pred) => {
                let is_match = left[pred.left_index] == right[pred.right_index];
                if is_match
                    && prefix + pred.left_index + 1 == c.left_start + c.length
                    && prefix + pred.right_index + 1 == c.right_start + c.length
                {
                    c.length += 1;
                    c.left_start -= 1;
                    c.right_start -= 1;
                } else {
                    if c.length > 0 {
                        matches.push(c);
                    }
                    c = CommonRange {
                        left_start: prefix + pred.left_index,
                        right_start: prefix + pred.right_index,
                        length: if is_match { 1 } else { 0 },
                    };
                }
                node = pred;
            }
            None => {
                if c.length > 0 {
                    matches.push(c);
                }
                break;
            }
        }
    }

    debug_assert!(*node == source);
    if prefix > 0 {
        matches[1..].reverse();
    } else {
        matches.reverse();
    }
    if suffix > 0 {
        matches.push(CommonRange {
            left_start: prefix + n,
            right_start: prefix + m,
            length: suffix,
        });
    }
    matches
}
