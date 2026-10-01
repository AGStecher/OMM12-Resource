#!/usr/bin/env python3
"""Neighbor-Joining tree construction (Saitou & Nei 1987), pure Python.

Builds an unrooted tree (final trifurcation at the base, the standard NJ
output), with a Newick serializer, a midpoint-rooting step for nicer
display, and a bipartition extractor used both for self-validation and for
bootstrap support tallying.
"""
import copy
import itertools


class Node:
    __slots__ = ("name", "children", "id")
    _counter = itertools.count()

    def __init__(self, name=None):
        self.name = name
        self.children = []  # list of [Node, branch_length]
        self.id = next(Node._counter)

    def is_leaf(self):
        return len(self.children) == 0

    def leaves(self):
        if self.is_leaf():
            return [self.name]
        out = []
        for child, _bl in self.children:
            out.extend(child.leaves())
        return out


def neighbor_joining(labels, D):
    """labels: list of names. D: dict-of-dicts distance matrix (symmetric,
    D[a][b] for a!=b). Returns the root Node of an unrooted NJ tree (final
    node has 3 children -- the standard unrooted-tree trifurcation)."""
    active = {name: Node(name) for name in labels}
    dist = {a: {b: D[a][b] for b in labels if b != a} for a in labels}

    while len(active) > 2:
        names = list(active.keys())
        r = len(names)
        total = {a: sum(dist[a].values()) for a in names}

        # find pair minimizing the Q criterion
        best_q = None
        best_pair = None
        for a, b in itertools.combinations(names, 2):
            q = (r - 2) * dist[a][b] - total[a] - total[b]
            if best_q is None or q < best_q:
                best_q = q
                best_pair = (a, b)
        a, b = best_pair

        dab = dist[a][b]
        delta_a = 0.5 * dab + (total[a] - total[b]) / (2 * (r - 2))
        delta_b = dab - delta_a
        # guard against tiny negative branch lengths from floating point noise
        delta_a = max(delta_a, 0.0) if abs(delta_a) < 1e-9 else delta_a
        delta_b = max(delta_b, 0.0) if abs(delta_b) < 1e-9 else delta_b

        u = Node(None)
        u.children.append([active[a], delta_a])
        u.children.append([active[b], delta_b])
        uname = f"__internal_{u.id}__"

        new_dist_u = {}
        for k in names:
            if k in (a, b):
                continue
            dku = 0.5 * (dist[a][k] + dist[b][k] - dab)
            new_dist_u[k] = dku

        # rebuild active/dist without a,b, with u added
        del active[a]; del active[b]
        for k in list(dist.keys()):
            if k in (a, b):
                del dist[k]
                continue
            del dist[k][a]
            del dist[k][b]
            dist[k][uname] = new_dist_u[k]
        dist[uname] = new_dist_u
        active[uname] = u

    # exactly 2 left -> join with a single edge, split evenly is arbitrary;
    # standard convention: the remaining distance is the branch between them.
    names = list(active.keys())
    a, b = names[0], names[1]
    dab = dist[a][b]
    root = Node(None)
    root.children.append([active[a], dab])
    # NOTE: with 2 nodes left this is really just one edge; we attach both
    # under a synthetic root with the full length on one side and 0 on the
    # other only if there isn't already a 3-way trifurcation path. To keep
    # the standard unrooted trifurcating-base convention (matches most NJ
    # implementations, e.g. ape::nj), we instead special-case: if root's
    # single child 'a' is already internal, graft b's single child onto it.
    if not active[a].is_leaf() and len(active[a].children) == 2:
        active[a].children.append([active[b], dab])
        root = active[a]
    else:
        root.children.append([active[b], 0.0])
    return root


def to_newick(node, support=None):
    """support: optional dict node.id -> bootstrap percentage (int), applied
    to internal nodes only."""
    def rec(n):
        if n.is_leaf():
            return n.name
        parts = [rec(c) + f":{bl:.6f}" for c, bl in n.children]
        s = "(" + ",".join(parts) + ")"
        if support is not None and n.id in support:
            s += str(support[n.id])
        return s
    return rec(node) + ";"


def bipartitions(node, all_leaves):
    """Return dict: frozenset(one side of the split, canonicalized to the
    side NOT containing all_leaves[0]) -> branch_length, for every internal
    edge (edges leading to a leaf are not included)."""
    anchor = all_leaves[0]
    out = {}

    def rec(n):
        if n.is_leaf():
            return {n.name}
        collected = set()
        for child, bl in n.children:
            child_leaves = rec(child)
            collected |= child_leaves
            if not child.is_leaf():
                side = frozenset(child_leaves)
                if anchor in side:
                    side = frozenset(set(all_leaves) - side)
                out[side] = bl
        return collected

    rec(node)
    return out


def patristic_distance(node, la, lb):
    """Path length between two leaves (for midpoint rooting / validation)."""
    def dfs(n, target, acc):
        if n.name == target:
            return acc
        for child, bl in n.children:
            r = dfs(child, target, acc + bl)
            if r is not None:
                return r
        return None

    # find path from root to la and to lb via parent tracking
    parent = {}
    bl_to_parent = {}

    def build_parent(n, p=None, bl=0.0):
        parent[n.id] = p
        bl_to_parent[n.id] = bl
        for child, cbl in n.children:
            build_parent(child, n, cbl)

    build_parent(node)

    def find(n, target):
        if n.name == target:
            return n
        for child, _ in n.children:
            r = find(child, target)
            if r is not None:
                return r
        return None

    na = find(node, la)
    nb = find(node, lb)

    def path_to_root(n):
        p = []
        cur = n
        while cur is not None:
            p.append(cur)
            cur = parent.get(cur.id)
        return p

    pa = path_to_root(na)
    pb = path_to_root(nb)
    ids_a = {n.id: i for i, n in enumerate(pa)}
    lca = None
    for n in pb:
        if n.id in ids_a:
            lca = n
            break
    d = 0.0
    for n in pa:
        if n.id == lca.id:
            break
        d += bl_to_parent[n.id]
    for n in pb:
        if n.id == lca.id:
            break
        d += bl_to_parent[n.id]
    return d


def bipartitions_with_nodes(node, all_leaves):
    """Like bipartitions(), but returns frozenset -> (Node, branch_length)
    so callers can map computed support values back onto specific nodes."""
    anchor = all_leaves[0]
    out = {}

    def rec(n):
        if n.is_leaf():
            return {n.name}
        collected = set()
        for child, bl in n.children:
            child_leaves = rec(child)
            collected |= child_leaves
            if not child.is_leaf():
                side = frozenset(child_leaves)
                if anchor in side:
                    side = frozenset(set(all_leaves) - side)
                out[side] = (child, bl)
        return collected

    rec(node)
    return out
