#!/usr/bin/env python3
"""Bootstrap support for the OMM12 16S NJ tree: resample MSA columns with
replacement, rebuild JC-distance + NJ tree per replicate, and tally, for
each internal split in the best tree, the fraction of replicate trees that
contain the identical bipartition (standard Felsenstein bootstrap).
"""
import math
import time

import numpy as np

from jc_distance import jc69_correct
from nj import neighbor_joining, bipartitions, bipartitions_with_nodes, to_newick

SEED = 20260909  # fixed for reproducibility; documented in the app caveat too
N_REPLICATES = 200


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


def main():
    msa, labels = read_fasta("omm12_16s_msa.fasta")
    n = len(labels)
    L = len(msa[labels[0]])

    M = np.array([[c for c in msa[name]] for name in labels])  # n x L char array
    gap = (M == "-")

    rng = np.random.default_rng(SEED)

    # best tree (recompute here from the un-resampled MSA for a self-contained script)
    def dist_from_matrix(Mx, gapx):
        D = {a: {} for a in labels}
        for i in range(n):
            for j in range(i + 1, n):
                valid = ~gapx[i] & ~gapx[j]
                ncomp = int(valid.sum())
                if ncomp == 0:
                    p = None
                else:
                    ndiff = int(((Mx[i] != Mx[j]) & valid).sum())
                    p = ndiff / ncomp
                d = jc69_correct(p) if p is not None else 3.0
                D[labels[i]][labels[j]] = d
                D[labels[j]][labels[i]] = d
            D[labels[i]][labels[i]] = 0.0
        return D

    D_best = dist_from_matrix(M, gap)
    best_tree = neighbor_joining(labels, D_best)
    best_bp = bipartitions_with_nodes(best_tree, labels)
    print(f"Best tree has {len(best_bp)} internal splits to support-test.")

    counts = {side: 0 for side in best_bp}

    t0 = time.time()
    for rep in range(N_REPLICATES):
        idx = rng.integers(0, L, size=L)
        Mr = M[:, idx]
        gapr = gap[:, idx]
        Dr = dist_from_matrix(Mr, gapr)
        tree_r = neighbor_joining(labels, Dr)
        bp_r = bipartitions(tree_r, labels)
        rep_sides = set(bp_r.keys())
        for side in counts:
            if side in rep_sides:
                counts[side] += 1
        if (rep + 1) % 50 == 0:
            print(f"  {rep + 1}/{N_REPLICATES} replicates done ({time.time()-t0:.1f}s elapsed)")

    print(f"\nAll {N_REPLICATES} bootstrap replicates done in {time.time()-t0:.1f}s")

    support = {}
    print("\nBootstrap support per internal split:")
    for side, (node, bl) in best_bp.items():
        pct = round(100 * counts[side] / N_REPLICATES)
        support[node.id] = pct
        print(f"  {sorted(side)}: {pct}%  (branch length {bl:.6f})")

    newick_with_support = to_newick(best_tree, support=support)
    print("\nFinal Newick with bootstrap support:")
    print(newick_with_support)

    with open("omm12_16s_tree_final.nwk", "w") as fh:
        fh.write(newick_with_support + "\n")
    print("\nWritten -> omm12_16s_tree_final.nwk")

    with open("omm12_16s_pipeline_report.txt", "w") as fh:
        fh.write("OMM12 16S rRNA gene phylogeny -- pipeline report\n")
        fh.write("=" * 60 + "\n\n")
        fh.write("Method: one representative 16S rRNA gene copy per genome (longest\n")
        fh.write("annotated copy), aligned with a custom Gotoh affine-gap-penalty\n")
        fh.write("center-star multiple alignment (medoid-sequence-anchored), pairwise\n")
        fh.write("distances corrected with the Jukes-Cantor (1969) substitution model,\n")
        fh.write("tree built by Neighbor-Joining (Saitou & Nei 1987), support from\n")
        fh.write(f"{N_REPLICATES} bootstrap replicates (column resampling with replacement,\n")
        fh.write(f"RNG seed {SEED} for reproducibility).\n\n")
        fh.write(f"Alignment: {n} taxa x {L} columns.\n")
        fh.write("No external alignment/tree tools were available in this environment\n")
        fh.write("(no internet access to install MAFFT/IQ-TREE/etc.) -- this is a from-\n")
        fh.write("scratch pure-Python/numpy implementation. Each component (pairwise\n")
        fh.write("aligner, JC correction, Neighbor-Joining) was unit-tested against\n")
        fh.write("known-correct reference cases before being run on real data; see\n")
        fh.write("nw_align.py, jc_distance.py and test_nj.py in this same folder.\n\n")
        fh.write("This is a real single-locus (16S) molecular phylogeny with alignment,\n")
        fh.write("a substitution model, and bootstrap support -- a genuine improvement\n")
        fh.write("over a k-mer/UPGMA similarity sketch. It is NOT a whole-genome GBDP\n")
        fh.write("tree (which would need many more loci / whole-genome distances) and,\n")
        fh.write("like any single-gene tree, has limited power to resolve the deepest,\n")
        fh.write("most ancient splits between very different phyla -- treat weakly\n")
        fh.write("supported deep branches (~<70% bootstrap) with appropriate caution.\n\n")
        fh.write("Final Newick (with bootstrap support as internal node labels):\n")
        fh.write(newick_with_support + "\n")


if __name__ == "__main__":
    main()
