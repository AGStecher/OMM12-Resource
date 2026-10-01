#!/usr/bin/env python3
"""Center-star multiple sequence alignment of the 12 OMM12 16S sequences.

1. Compute all C(12,2)=66 pairwise global alignments (Gotoh affine gaps).
2. Pick the medoid sequence (minimizes total p-distance to the other 11) as
   the alignment center.
3. Merge the 11 pairwise alignments (center vs each other sequence) into one
   MSA by the standard center-star merge: at each "slot" between consecutive
   center residues, take the longest insertion seen in any pairwise
   alignment at that slot, and pad every sequence's insertion at that slot
   with trailing gaps up to that length.
"""
import itertools
import sys
import time

from nw_align import align_affine


def read_fasta(path):
    seqs = {}
    order = []
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].strip()
                order.append(name)
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs, order


def p_distance(ax, ay):
    """Fraction of mismatched sites among columns where NEITHER sequence has
    a gap (pairwise deletion)."""
    n_compared = 0
    n_diff = 0
    for a, b in zip(ax, ay):
        if a == "-" or b == "-":
            continue
        n_compared += 1
        if a != b:
            n_diff += 1
    if n_compared == 0:
        return None
    return n_diff / n_compared, n_compared


def build_star_msa(center_seq, others):
    """others: dict of name -> (aligned_center_form, aligned_other_form) all
    against `center_seq`. Returns dict name -> final aligned string
    (including an entry "__center__" for the center itself) all of equal
    length.
    """
    L = len(center_seq)
    # For each other sequence, compute per-slot insertion strings (slot 0..L)
    per_seq_slot_ins = {}   # name -> list of length L+1, each an insertion string
    per_seq_slot_core = {}  # name -> list of length L, char aligned to that center residue
    for name, (ac, ao) in others.items():
        slot_ins = [""] * (L + 1)
        slot_core = [None] * L
        p = 0  # index of next center residue (0-based)
        col = 0
        buf = []
        while col < len(ac):
            if ac[col] != "-":
                # flush any pending insertion into slot p
                slot_ins[p] += "".join(buf)
                buf = []
                slot_core[p] = ao[col]
                p += 1
            else:
                buf.append(ao[col])
            col += 1
        # trailing insertion after the last residue -> slot L
        slot_ins[L] += "".join(buf)
        assert p == L, f"{name}: consumed {p} of {L} center residues"
        per_seq_slot_ins[name] = slot_ins
        per_seq_slot_core[name] = slot_core

    max_insert = [0] * (L + 1)
    for name, slot_ins in per_seq_slot_ins.items():
        for s in range(L + 1):
            if len(slot_ins[s]) > max_insert[s]:
                max_insert[s] = len(slot_ins[s])

    # build final center row
    out = {}
    parts = []
    for s in range(L + 1):
        parts.append("-" * max_insert[s])
        if s < L:
            parts.append(center_seq[s])
    out["__center__"] = "".join(parts)

    for name in others:
        slot_ins = per_seq_slot_ins[name]
        slot_core = per_seq_slot_core[name]
        parts = []
        for s in range(L + 1):
            ins = slot_ins[s]
            pad = max_insert[s] - len(ins)
            parts.append(ins + "-" * pad)
            if s < L:
                parts.append(slot_core[s])
        out[name] = "".join(parts)

    return out


def main():
    fasta_path = "/sessions/upbeat-wizardly-ramanujan/mnt/outputs/omm12_tree/omm12_16s_raw.fasta"
    seqs, order = read_fasta(fasta_path)
    names = order
    n = len(names)
    print(f"Loaded {n} sequences")

    t0 = time.time()
    pairwise = {}   # (i,j) i<j -> (ax, ay, score)
    dist = {n_: {} for n_ in names}
    for i, j in itertools.combinations(range(n), 2):
        a, b = names[i], names[j]
        ax, ay, sc = align_affine(seqs[a], seqs[b])
        pairwise[(a, b)] = (ax, ay, sc)
        pd, ncomp = p_distance(ax, ay)
        dist[a][b] = pd
        dist[b][a] = pd
    print(f"All {len(pairwise)} pairwise alignments done in {time.time()-t0:.1f}s")

    # pick medoid = center
    totals = {a: sum(dist[a][b] for b in names if b != a) for a in names}
    center = min(totals, key=totals.get)
    print(f"Center (medoid) sequence: {center}  (sum p-distance to others = {totals[center]:.4f})")
    print("Per-bacterium total p-distance to all others (sanity check, sorted):")
    for a, tot in sorted(totals.items(), key=lambda kv: kv[1]):
        print(f"  {a}: {tot:.4f}")

    others = {}
    for name in names:
        if name == center:
            continue
        if (center, name) in pairwise:
            ax, ay, _ = pairwise[(center, name)]
        else:
            ay, ax, _ = pairwise[(name, center)]
        others[name] = (ax, ay)

    msa = build_star_msa(seqs[center], others)
    msa[center] = msa.pop("__center__")

    lengths = {k: len(v) for k, v in msa.items()}
    assert len(set(lengths.values())) == 1, f"MSA columns not uniform: {lengths}"
    L = next(iter(lengths.values()))
    print(f"\nFinal MSA length: {L} columns, {len(msa)} sequences")

    # --- correctness checks: degapping each aligned row must reconstruct the original sequence exactly
    fail = []
    for name in names:
        degapped = msa[name].replace("-", "")
        if degapped != seqs[name]:
            fail.append(name)
    if fail:
        print("FAIL: degapped MSA rows don't match originals for:", fail)
        sys.exit(1)
    print("All 12 rows verified: degapping the MSA exactly reconstructs the original extracted sequence.")

    out_path = "/sessions/upbeat-wizardly-ramanujan/mnt/outputs/omm12_tree/omm12_16s_msa.fasta"
    with open(out_path, "w") as fh:
        for name in names:
            fh.write(f">{name}\n")
            s = msa[name]
            for i in range(0, len(s), 70):
                fh.write(s[i:i + 70] + "\n")
    print(f"MSA written -> {out_path}")

    # quick composition sanity: fraction of fully-gap columns (shouldn't be huge)
    all_gap_cols = sum(1 for c in range(L) if all(msa[name][c] == "-" for name in names))
    print(f"Fully-gap columns (should be 0): {all_gap_cols}")


if __name__ == "__main__":
    main()
