#!/usr/bin/env python3
"""Jukes-Cantor (1969) corrected pairwise distance from an MSA."""
import math


def p_distance(seq_a, seq_b):
    """Pairwise-deletion p-distance: fraction mismatched among columns
    where NEITHER sequence has a gap. Returns (p, n_sites_compared)."""
    n_compared = 0
    n_diff = 0
    for a, b in zip(seq_a, seq_b):
        if a == "-" or b == "-":
            continue
        n_compared += 1
        if a != b:
            n_diff += 1
    if n_compared == 0:
        return None, 0
    return n_diff / n_compared, n_compared


def jc69_correct(p):
    """Jukes-Cantor 1969 distance correction. Undefined (returns None) for
    p >= 0.75 (saturation -- infinite distance under the JC model)."""
    if p is None:
        return None
    if p >= 0.75:
        return None
    return -0.75 * math.log(1 - (4.0 / 3.0) * p)


def distance_matrix(msa, labels, correct=True):
    """msa: dict name->aligned seq. Returns dict-of-dicts distance matrix,
    plus a report of any pair that hit JC saturation (p>=0.75, distance
    undefined) so the caller can decide how to handle it."""
    D = {a: {} for a in labels}
    saturated = []
    raw_p = {a: {} for a in labels}
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            p, ncomp = p_distance(msa[a], msa[b])
            raw_p[a][b] = raw_p[b][a] = p
            d = jc69_correct(p) if correct else p
            if d is None:
                saturated.append((a, b, p))
                d = 3.0  # cap: treat as "very far" rather than crash; flagged separately
            D[a][b] = d
            D[b][a] = d
        D[a][a] = 0.0
    return D, saturated, raw_p


if __name__ == "__main__":
    # sanity checks against hand-computed JC values
    import math as _m
    tests = [
        (0.0, 0.0),
        (0.1, -0.75 * _m.log(1 - 4/3*0.1)),
        (0.5, -0.75 * _m.log(1 - 4/3*0.5)),
    ]
    for p, expected in tests:
        got = jc69_correct(p)
        assert math.isclose(got, expected, rel_tol=1e-9), (p, got, expected)
        print(f"p={p} -> JC d={got:.6f} (matches direct formula)")
    assert jc69_correct(0.75) is None
    assert jc69_correct(0.9) is None
    print("JC correction sanity checks passed.")
