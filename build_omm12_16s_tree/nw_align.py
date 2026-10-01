#!/usr/bin/env python3
"""Global pairwise nucleotide alignment with affine gap penalties (Gotoh's
algorithm, three-state formulation as in Durbin et al., Biological Sequence
Analysis, ch. 2). Pure numpy + Python stdlib, no external alignment tools.
"""
import math
import numpy as np

NEG_INF = -1e18


def align_affine(x, y, match=2.0, mismatch=-1.0, gap_open=10.0, gap_extend=0.5):
    """Global alignment of x, y. Gap of length L costs gap_open + L*gap_extend
    (i.e. gap_open is charged once when a gap starts, on top of the per-
    residue gap_extend that every gap position -- including the first --
    pays). Returns (aligned_x, aligned_y, score).
    """
    n, m = len(x), len(y)
    xa = np.frombuffer(x.encode("ascii"), dtype=np.uint8)
    ya = np.frombuffer(y.encode("ascii"), dtype=np.uint8)

    M = np.full((n + 1, m + 1), NEG_INF)
    Ix = np.full((n + 1, m + 1), NEG_INF)  # gap in y (x consumed, y not)
    Iy = np.full((n + 1, m + 1), NEG_INF)  # gap in x (y consumed, x not)
    M[0, 0] = 0.0
    for i in range(1, n + 1):
        Ix[i, 0] = -gap_open - (i - 1) * gap_extend
    for j in range(1, m + 1):
        Iy[0, j] = -gap_open - (j - 1) * gap_extend

    for i in range(1, n + 1):
        row_match = np.where(xa[i - 1] == ya, match, mismatch).astype(np.float64)
        prev_best = np.maximum(np.maximum(M[i - 1, 0:m], Ix[i - 1, 0:m]), Iy[i - 1, 0:m])
        M[i, 1:m + 1] = row_match + prev_best

        Ix[i, 1:m + 1] = np.maximum(M[i - 1, 1:m + 1] - gap_open, Ix[i - 1, 1:m + 1] - gap_extend)

        # Iy has a sequential (within-row) dependency -> plain loop
        Iy_row = Iy[i]
        M_row = M[i]
        prev_iy = NEG_INF
        for j in range(1, m + 1):
            cand = max(M_row[j - 1] - gap_open, prev_iy - gap_extend)
            Iy_row[j] = cand
            prev_iy = cand

    # pick best final state
    finals = {"M": M[n, m], "Ix": Ix[n, m], "Iy": Iy[n, m]}
    state = max(finals, key=finals.get)
    score = finals[state]

    ax, ay = [], []
    i, j = n, m
    tol = 1e-6
    while i > 0 or j > 0:
        if state == "M":
            s = match if x[i - 1] == y[j - 1] else mismatch
            target = M[i, j] - s
            cands = {"M": M[i - 1, j - 1], "Ix": Ix[i - 1, j - 1], "Iy": Iy[i - 1, j - 1]}
            ax.append(x[i - 1]); ay.append(y[j - 1])
            i -= 1; j -= 1
        elif state == "Ix":
            cands = {"M": M[i - 1, j] - gap_open, "Ix": Ix[i - 1, j] - gap_extend}
            target = Ix[i, j]
            ax.append(x[i - 1]); ay.append("-")
            i -= 1
        else:  # Iy
            cands = {"M": M[i, j - 1] - gap_open, "Iy": Iy[i, j - 1] - gap_extend}
            target = Iy[i, j]
            ax.append("-"); ay.append(y[j - 1])
            j -= 1
        # pick predecessor whose value matches target (within tolerance)
        state = min(cands, key=lambda k: abs(cands[k] - target))
        if not math.isclose(cands[state], target, abs_tol=tol, rel_tol=1e-6):
            # shouldn't happen with consistent float64 arithmetic, but don't
            # silently produce a wrong alignment if it does
            raise RuntimeError(f"Traceback mismatch at i={i},j={j}: {cands} vs {target}")

    return "".join(reversed(ax)), "".join(reversed(ay)), float(score)


if __name__ == "__main__":
    # --- toy validation tests ---
    failures = []

    # 1. Identical sequences -> perfect match, no gaps
    ax, ay, sc = align_affine("ACGTACGT", "ACGTACGT")
    print("Test 1 (identical):", ax, "|", ay, "score=", sc)
    if ax != "ACGTACGT" or ay != "ACGTACGT" or "-" in ax or "-" in ay:
        failures.append("Test 1 failed: expected gap-free identical alignment")

    # 2. Single mismatch, equal length, no gaps expected
    ax, ay, sc = align_affine("AAAA", "AAAT")
    print("Test 2 (single mismatch):", ax, "|", ay, "score=", sc)
    if "-" in ax or "-" in ay:
        failures.append("Test 2 failed: expected no gaps for a single internal mismatch")
    if sum(1 for a, b in zip(ax, ay) if a != b) != 1:
        failures.append("Test 2 failed: expected exactly one mismatch column")

    # 3. Contiguous 8bp deletion should align as ONE gap, not scattered
    x3 = "AAAACCCCTTTTGGGG"
    y3 = "AAAAGGGG"
    ax, ay, sc = align_affine(x3, y3)
    print("Test 3 (contiguous deletion):")
    print("  x:", ax)
    print("  y:", ay)
    # count gap *runs* in ay (should be 1 run of length 8)
    runs = 0
    prev = None
    for ch in ay:
        if ch == "-" and prev != "-":
            runs += 1
        prev = ch
    n_gaps = ay.count("-")
    print(f"  gap runs={runs}, total gap chars={n_gaps}")
    if runs != 1 or n_gaps != 8:
        failures.append(f"Test 3 failed: expected 1 contiguous gap run of length 8, got runs={runs} len={n_gaps}")
    # and the ungapped columns should exactly reconstruct "AAAA"+"GGGG" matched to y
    reconstructed = "".join(a for a, b in zip(ax, ay) if b != "-")
    if reconstructed != y3:
        failures.append(f"Test 3 failed: reconstructed non-gap x-columns {reconstructed!r} != y {y3!r}")

    # 4. Symmetry: aligning (x,y) vs (y,x) should give the same score
    _, _, sc_a = align_affine("GATTACAGATTACA", "GATCACAGATCACA")
    _, _, sc_b = align_affine("GATCACAGATCACA", "GATTACAGATTACA")
    print("Test 4 (symmetry): score_xy=", sc_a, "score_yx=", sc_b)
    if not math.isclose(sc_a, sc_b, rel_tol=1e-9):
        failures.append("Test 4 failed: alignment score not symmetric under swapping x,y")

    print()
    if failures:
        print("FAILURES:")
        for f in failures:
            print(" -", f)
        raise SystemExit(1)
    else:
        print("All alignment validation tests passed.")
