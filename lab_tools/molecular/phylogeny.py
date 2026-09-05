#!/usr/bin/env python3
"""Phylogenetic Tree CLI – Neighbor-Joining tree builder."""

import argparse
import sys
from typing import Dict, List, Tuple, Optional


def parse_phylip(text: str) -> Tuple[List[str], List[List[float]]]:
    """Parse a PHYLIP distance matrix (square or lower-triangular).

    Returns:
        labels: ordered list of taxon names.
        matrix: full symmetric distance matrix as list of lists.

    Raises:
        ValueError: on malformed input.
    """
    lines = [line.strip() for line in text.strip().split('\n') if line.strip() and not line.strip().startswith('#')]
    if not lines:
        raise ValueError("Empty input.")
    try:
        n = int(lines[0].split()[0])
    except (IndexError, ValueError):
        raise ValueError("First line must be an integer (number of taxa).")
    if n < 2:
        raise ValueError("Need at least 2 taxa.")

    labels = []
    rows = []

    for i, line in enumerate(lines[1:1+n], start=1):
        parts = line.split()
        if len(parts) < 2:
            raise ValueError(f"Line {i}: expected label and at least one distance.")
        label = parts[0]
        if label in labels:
            raise ValueError(f"Duplicate label '{label}' on line {i}.")
        labels.append(label)
        distances = []
        for j, val in enumerate(parts[1:]):
            try:
                distances.append(float(val))
            except ValueError:
                raise ValueError(f"Line {i}: invalid distance '{val}'.")
        rows.append(distances)

    if len(rows) != n:
        raise ValueError(f"Expected {n} data lines, got {len(rows)}.")

    # Determine if matrix is lower-triangular (row length increases) or square
    is_lower_tri = (len(rows[0]) == 1 and len(rows[-1]) == n)
    is_square = all(len(r) == n for r in rows)

    if is_square:
        matrix = [r[:] for r in rows]
    elif is_lower_tri:
        # Build full symmetric matrix
        matrix = [[0.0]*n for _ in range(n)]
        for i in range(n):
            if len(rows[i]) != i+1:
                raise ValueError(f"Line {i+1}: expected {i+1} distances for lower-triangular matrix, got {len(rows[i])}.")
            for j in range(i+1):
                matrix[i][j] = rows[i][j]
                matrix[j][i] = rows[i][j]
    else:
        # Possibly ragged; try to infer as lower-tri if last row has n elements
        if len(rows[-1]) == n and rows[0][0] == 0.0:
            # Rebuild as lower triangular
            matrix = [[0.0]*n for _ in range(n)]
            for i in range(n):
                for j, val in enumerate(rows[i]):
                    if j > i:
                        raise ValueError(f"Line {i+1}: too many distances for lower-triangular format.")
                    matrix[i][j] = val
                    matrix[j][i] = val
        else:
            raise ValueError("Matrix format not recognised; must be square or lower-triangular.")

    # Validate symmetry
    for i in range(n):
        for j in range(i+1, n):
            if abs(matrix[i][j] - matrix[j][i]) > 1e-6:
                raise ValueError(f"Matrix not symmetric: D[{i},{j}]={matrix[i][j]} vs D[{j},{i}]={matrix[j][i]}.")
        if abs(matrix[i][i]) > 1e-6:
            raise ValueError(f"Diagonal element {i} must be zero.")

    return labels, matrix


class NJNode:
    """Node in the tree."""
    __slots__ = ('label', 'children', 'length')

    def __init__(self, label: Optional[str] = None, length: float = 0.0):
        self.label = label
        self.children: List['NJNode'] = []
        self.length = length

    def is_leaf(self) -> bool:
        return not self.children


def neighbor_joining(labels: List[str], D: List[List[float]]) -> NJNode:
    """Build a tree using the Neighbor-Joining algorithm.

    Args:
        labels: taxon names.
        D: symmetric distance matrix (n x n).

    Returns:
        Root node of the unrooted tree (placed at the last join).
    """
    n = len(labels)
    if n == 2:
        # Just create a node linking the two leaves
        root = NJNode()
        a = NJNode(labels[0], D[0][1] / 2.0)
        b = NJNode(labels[1], D[0][1] / 2.0)
        root.children = [a, b]
        return root

    # Work with mutable copies
    d = [row[:] for row in D]
    active = list(range(n))          # indices still in the matrix
    nodes: Dict[int, NJNode] = {i: NJNode(labels[i]) for i in range(n)}

    while len(active) > 2:
        m = len(active)
        # Compute S_i = sum of distances to all other active taxa / (m - 2)
        S = [0.0] * m
        for i_idx, i in enumerate(active):
            total = sum(d[i_idx][j_idx] for j_idx in range(m) if j_idx != i_idx)
            S[i_idx] = total / (m - 2) if m > 2 else 0.0

        # Find pair (i,j) with minimum Q = (m-2)*d[i][j] - S[i] - S[j]
        min_q = float('inf')
        pair = (0, 1)
        for i_idx in range(m):
            for j_idx in range(i_idx + 1, m):
                q = (m - 2) * d[i_idx][j_idx] - S[i_idx] - S[j_idx]
                if q < min_q:
                    min_q = q
                    pair = (i_idx, j_idx)

        i_idx, j_idx = pair
        i_orig = active[i_idx]
        j_orig = active[j_idx]

        # Branch lengths
        dist_ij = d[i_idx][j_idx]
        if m > 2:
            len_i = (dist_ij + (S[i_idx] - S[j_idx]) / (m - 2)) / 2.0
        else:
            len_i = dist_ij / 2.0
        len_j = dist_ij - len_i

        # Create internal node, attaching the existing nodes so subtrees are kept
        new_node = NJNode()
        child_i = nodes[i_orig]
        child_j = nodes[j_orig]
        child_i.length = max(len_i, 0.0)
        child_j.length = max(len_j, 0.0)
        new_node.children = [child_i, child_j]

        # Update distance matrix: new row for (i,j) combined
        new_dist_row = []
        for k_idx, k in enumerate(active):
            if k_idx == i_idx or k_idx == j_idx:
                continue
            dik = d[i_idx][k_idx]
            djk = d[j_idx][k_idx]
            d_new = (dik + djk - dist_ij) / 2.0
            new_dist_row.append(max(d_new, 0.0))

        # Rebuild active list and distance matrix
        new_active = [k for idx, k in enumerate(active) if idx not in (i_idx, j_idx)]
        new_idx = len(new_active)  # index of the new internal node
        new_active.append(max(active) + 1)  # temporary unique index
        new_d = [[0.0] * (m - 1) for _ in range(m - 1)]

        # Copy old distances for remaining taxa
        old_to_new = {}
        new_pos = 0
        for old_idx, k in enumerate(active):
            if old_idx in (i_idx, j_idx):
                continue
            old_to_new[old_idx] = new_pos
            new_pos += 1
        # fill the new matrix
        for a_old_idx, a_new_idx in old_to_new.items():
            for b_old_idx, b_new_idx in old_to_new.items():
                new_d[a_new_idx][b_new_idx] = d[a_old_idx][b_old_idx]
            # distance to new node
            new_d[a_new_idx][new_idx] = new_dist_row[a_new_idx]
            new_d[new_idx][a_new_idx] = new_dist_row[a_new_idx]

        # Register new node
        nodes[len(nodes)] = new_node
        # Replace the two joined indices
        active = [active[k] for k in range(m) if k not in (i_idx, j_idx)]
        active.append(len(nodes) - 1)
        d = new_d

    # Final two: join them
    if len(active) == 2:
        i_idx, j_idx = 0, 1
        i_orig = active[i_idx]
        j_orig = active[j_idx]
        root = NJNode()
        nodes[i_orig].length = d[i_idx][j_idx] / 2.0
        nodes[j_orig].length = d[i_idx][j_idx] / 2.0
        root.children = [nodes[i_orig], nodes[j_orig]]
        return root

    # Should never reach
    raise RuntimeError("Unexpected state in NJ algorithm.")


def to_newick(node: NJNode) -> str:
    """Convert an NJNode tree to Newick format with branch lengths."""
    if node.is_leaf():
        return f"{node.label}:{node.length:.6f}"
    children_str = ','.join(to_newick(child) for child in node.children)
    return f"({children_str})"


def ascii_tree(node: NJNode, prefix: str = "", is_last: bool = True) -> str:
    """Generate a simple ASCII representation of the tree."""
    if node.is_leaf():
        connector = "└── " if is_last else "├── "
        return prefix + connector + f"{node.label} ({node.length:.4f})\n"
    result = ""
    if node.label is None and node.length == 0.0:
        # internal node, could print as [+] but hide for root
        pass
    child_count = len(node.children)
    for i, child in enumerate(node.children):
        child_is_last = (i == child_count - 1)
        connector = "└── " if child_is_last else "├── "
        if node.length != 0.0 or node.label is not None:
            # this internal node has length or label (only in final tree root is not printed)
            result += prefix + connector + f"({(node.label or '+')}:{node.length:.4f})\n"
            new_pref = prefix + ("    " if child_is_last else "│   ")
        else:
            new_pref = prefix + ("    " if child_is_last else "│   ")
        result += ascii_tree(child, new_pref, child_is_last)
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Build a phylogenetic tree from a distance matrix using Neighbor-Joining."
    )
    parser.add_argument(
        "input",
        help="Path to distance matrix file (PHYLIP format)."
    )
    parser.add_argument(
        "-o", "--output",
        help="Write Newick tree to file instead of stdout."
    )
    parser.add_argument(
        "--ascii",
        action="store_true",
        help="Print an ASCII representation of the tree."
    )
    args = parser.parse_args(argv)

    try:
        with open(args.input) as f:
            text = f.read()
    except IOError as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        return 1

    try:
        labels, matrix = parse_phylip(text)
    except ValueError as e:
        print(f"Error parsing matrix: {e}", file=sys.stderr)
        return 1

    tree = neighbor_joining(labels, matrix)
    newick_str = to_newick(tree) + ";"

    if args.output:
        try:
            with open(args.output, 'w') as out:
                out.write(newick_str + "\n")
        except IOError as e:
            print(f"Error writing output: {e}", file=sys.stderr)
            return 1
    else:
        print(newick_str)

    if args.ascii:
        print(ascii_tree(tree))

    return 0


if __name__ == "__main__":
    sys.exit(main())
