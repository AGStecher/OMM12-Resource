#!/usr/bin/env python3
import sys
from jc_distance import distance_matrix, p_distance
from nj import neighbor_joining, to_newick, bipartitions


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


msa, labels = read_fasta("omm12_16s_msa.fasta")
L = len(next(iter(msa.values())))
print(f"MSA: {len(labels)} sequences x {L} columns")

# alignment stats
n_var = 0
n_gapfree = 0
for c in range(L):
    col = [msa[name][c] for name in labels]
    if "-" not in col:
        n_gapfree += 1
    if len(set(col)) > 1:
        n_var += 1
print(f"Variable columns: {n_var}/{L}  |  gap-free columns: {n_gapfree}/{L}")

D, saturated, raw_p = distance_matrix(msa, labels, correct=True)
print(f"\nSaturated pairs (p>=0.75, JC undefined): {saturated}")

print("\nJC-corrected distance matrix:")
print("".ljust(24), *[l[:10].ljust(11) for l in labels])
for a in labels:
    print(a[:22].ljust(24), *[f"{D[a][b]:.4f}".ljust(11) for b in labels])

tree = neighbor_joining(labels, D)
newick = to_newick(tree)
print("\nBest-tree Newick:")
print(newick)

with open("omm12_16s_besttree.nwk", "w") as fh:
    fh.write(newick + "\n")
print("\nWritten -> omm12_16s_besttree.nwk")
