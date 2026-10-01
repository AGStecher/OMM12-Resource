#!/usr/bin/env python3
"""Extract one representative 16S rRNA gene sequence per OMM12 bacterium.

Picks the longest annotated 16S rRNA feature per genome (Prokka/barrnap
usually annotates several near-identical rRNA operon copies; the longest
copy is least likely to be a partial/edge-of-contig fragment). Handles the
known acutalibacter_muris_kb18 GFF quirk where the COG-enriched
"_prokka_1.gff" has its seqid column corrupted (set to each feature's own
locus_tag instead of the real contig accession) by falling back to the
plain "_prokka.gff" whenever the enriched file's seqid doesn't match any
contig name actually present in the genome FASTA.
"""
import os
import re
import sys

DATA_DIR = "/sessions/upbeat-wizardly-ramanujan/mnt/OMM12_website/data"

bacteria_filenames = [
    "acutalibacter_muris_kb18",
    "akkermansia_muciniphila_YL44",
    "bacteroides_caecimuris_I48",
    "bifidobacterium_animalis_YL2",
    "blautia_coccoides_YL58",
    "clostridium_innocuum_I46",
    "enterocloster_clostridioformis_YL32",
    "enterococcus_faecalis_KB1",
    "flavonifractor_plautii_YL31",
    "limosilactobacillus_reuteri_I49",
    "muribaculum_intestinales_YL27",
    "turicimonas_muris_YL45",
]

bacteria_names = {
    "acutalibacter_muris_kb18": "Acutalibacter muris KB18",
    "akkermansia_muciniphila_YL44": "Akkermansia muciniphila YL44",
    "bacteroides_caecimuris_I48": "Bacteroides caecimuris I48",
    "bifidobacterium_animalis_YL2": "Bifidobacterium animalis YL2",
    "blautia_coccoides_YL58": "Blautia coccoides YL58",
    "clostridium_innocuum_I46": "Clostridium innocuum I46",
    "enterocloster_clostridioformis_YL32": "Enterocloster clostridioformis YL32",
    "enterococcus_faecalis_KB1": "Enterococcus faecalis KB1",
    "flavonifractor_plautii_YL31": "Flavonifractor plautii YL31",
    "limosilactobacillus_reuteri_I49": "Limosilactobacillus reuteri I49",
    "muribaculum_intestinales_YL27": "Muribaculum intestinales YL27",
    "turicimonas_muris_YL45": "Turicimonas muris YL45",
}

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def parse_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = "".join(chunks).upper()
    return seqs


def parse_gff_rrna(path):
    """Return list of dicts: seqid, start, end, strand, locus_tag, product for rRNA rows."""
    rows = []
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                break  # hit an embedded ##FASTA section's header lines
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts[:9]
            if ftype != "rRNA":
                continue
            m_product = re.search(r"product=([^;]+)", attrs)
            m_locus = re.search(r"locus_tag=([^;]+)", attrs)
            product = m_product.group(1) if m_product else ""
            locus_tag = m_locus.group(1) if m_locus else ""
            if "16S" not in product:
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            rows.append({
                "seqid": seqid, "start": start_i, "end": end_i,
                "strand": strand, "locus_tag": locus_tag, "product": product,
            })
    return rows


def best_gff_for(bacterium):
    """Prefer the COG-enriched GFF, but fall back to the plain one if its
    seqid values don't actually match any contig in the genome FASTA
    (the acutalibacter_muris_kb18 corruption case)."""
    enriched = os.path.join(DATA_DIR, f"{bacterium}_prokka_1.gff")
    plain = os.path.join(DATA_DIR, f"{bacterium}_prokka.gff")
    fasta_path = os.path.join(DATA_DIR, f"{bacterium}_genome.fasta")
    contigs = parse_fasta(fasta_path)
    contig_names = set(contigs.keys())

    candidate = enriched if os.path.exists(enriched) else plain
    if not os.path.exists(candidate):
        return None, contigs

    rows = parse_gff_rrna(candidate)
    if rows and not any(r["seqid"] in contig_names for r in rows):
        # seqid column looks corrupted (never matches a real contig) -- fall back
        if os.path.exists(plain) and plain != candidate:
            rows_plain = parse_gff_rrna(plain)
            if rows_plain and any(r["seqid"] in contig_names for r in rows_plain):
                return rows_plain, contigs
    return rows, contigs


def main():
    out_fasta = []
    report_lines = []
    for bf in bacteria_filenames:
        rows, contigs = best_gff_for(bf)
        if not rows:
            report_lines.append(f"{bf}: NO 16S rRNA rows found")
            continue
        # keep only rows whose seqid is an actual contig (defensive)
        rows = [r for r in rows if r["seqid"] in contigs]
        if not rows:
            report_lines.append(f"{bf}: 16S rows found but none matched a real contig")
            continue
        # pick the longest copy
        for r in rows:
            r["length"] = r["end"] - r["start"] + 1
        rows.sort(key=lambda r: r["length"], reverse=True)
        best = rows[0]
        contig_seq = contigs[best["seqid"]]
        s, e = best["start"], best["end"]
        if s < 1 or e > len(contig_seq) or s > e:
            report_lines.append(f"{bf}: chosen 16S row has out-of-bounds coordinates ({s}-{e}, contig len {len(contig_seq)})")
            continue
        seq = contig_seq[s - 1:e]
        if best["strand"] == "-":
            seq = revcomp(seq)
        out_fasta.append((bf, seq, best))
        report_lines.append(
            f"{bf}: {len(rows)} copies found, using locus_tag={best['locus_tag']} "
            f"len={best['length']}bp contig={best['seqid']} strand={best['strand']}"
        )

    fasta_path = "/sessions/upbeat-wizardly-ramanujan/mnt/outputs/omm12_tree/omm12_16s_raw.fasta"
    with open(fasta_path, "w") as fh:
        for bf, seq, meta in out_fasta:
            fh.write(f">{bf}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")

    print(f"Extracted {len(out_fasta)}/{len(bacteria_filenames)} sequences -> {fasta_path}\n")
    print("\n".join(report_lines))

    lengths = [len(seq) for _, seq, _ in out_fasta]
    if lengths:
        print(f"\nLength range: {min(lengths)}-{max(lengths)} bp (n={len(lengths)})")


if __name__ == "__main__":
    main()
