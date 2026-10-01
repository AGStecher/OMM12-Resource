#!/usr/bin/env python3
"""Validate neighbor_joining() against a hand-constructed, exactly-additive
distance matrix derived from a known tree with known branch lengths. NJ is
mathematically guaranteed to exactly recover the topology and branch
lengths for any perfectly additive distance matrix, so this is a strong,
self-contained correctness check (no reliance on recalled literature
numbers).

Ground-truth (unrooted) tree, as nested Newick for reference:
  ((A:1,B:2):3,(C:1.5,D:2.5):2,((E:1,F:1):1,G:4):2);
Pairwise leaf distances below were hand-computed by summing branch lengths
along the unique path between each pair (shown in the comments).
"""
import math
from nj import neighbor_joining, bipartitions, to_newick

labels = ["A", "B", "C", "D", "E", "F", "G"]

# hand-computed exact patristic distances (see module docstring / derivation)
raw = {
    "A": {"A": 0, "B": 3, "C": 7.5, "D": 8.5, "E": 8, "F": 8, "G": 10},
    "B": {"A": 3, "B": 0, "C": 8.5, "D": 9.5, "E": 9, "F": 9, "G": 11},
    "C": {"A": 7.5, "B": 8.5, "C": 0, "D": 4, "E": 7.5, "F": 7.5, "G": 9.5},
    "D": {"A": 8.5, "B": 9.5, "C": 4, "D": 0, "E": 8.5, "F": 8.5, "G": 10.5},
    "E": {"A": 8, "B": 9, "C": 7.5, "D": 8.5, "E": 0, "F": 2, "G": 6},
    "F": {"A": 8, "B": 9, "C": 7.5, "D": 8.5, "E": 2, "F": 0, "G": 6},
    "G": {"A": 10, "B": 11, "C": 9.5, "D": 10.5, "E": 6, "F": 6, "G": 0},
}

tree = neighbor_joining(labels, raw)
print("Reconstructed Newick:", to_newick(tree))

bp = bipartitions(tree, labels)
print("\nReconstructed internal-edge bipartitions (canonicalized, side w/o 'A'):")
for side, bl in bp.items():
    print(f"  {sorted(side)}: {bl:.6f}")

expected = {
    frozenset(["C", "D"]): 2.0,
    frozenset(["E", "F", "G"]): 2.0,
    frozenset(["E", "F"]): 1.0,
    frozenset(["C", "D", "E", "F", "G"]): 3.0,
}

failures = []
if set(bp.keys()) != set(expected.keys()):
    failures.append(f"Bipartition SET mismatch.\n  got:      {[sorted(s) for s in bp]}\n  expected: {[sorted(s) for s in expected]}")
else:
    for side, exp_bl in expected.items():
        got_bl = bp[side]
        if not math.isclose(got_bl, exp_bl, abs_tol=1e-6):
            failures.append(f"Branch length mismatch for {sorted(side)}: got {got_bl}, expected {exp_bl}")

# also check pendant (leaf) branch lengths by inspecting the tree directly
def leaf_branch_lengths(node, out):
    for child, bl in node.children:
        if child.is_leaf():
            out[child.name] = bl
        else:
            leaf_branch_lengths(child, out)

leaf_bls = {}
leaf_branch_lengths(tree, leaf_bls)
expected_leaf_bls = {"A": 1.0, "B": 2.0, "C": 1.5, "D": 2.5, "E": 1.0, "F": 1.0, "G": 4.0}
print("\nPendant (leaf) branch lengths:")
for name in labels:
    print(f"  {name}: got={leaf_bls.get(name)}, expected={expected_leaf_bls[name]}")
    if not math.isclose(leaf_bls.get(name, float('nan')), expected_leaf_bls[name], abs_tol=1e-6):
        failures.append(f"Leaf branch length mismatch for {name}: got {leaf_bls.get(name)}, expected {expected_leaf_bls[name]}")

print()
if failures:
    print("FAILURES:")
    for f in failures:
        print(" -", f)
    raise SystemExit(1)
else:
    print("PASS: neighbor_joining() exactly reconstructed the known additive tree (topology + all branch lengths).")
