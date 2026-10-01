#!/usr/bin/env python3
"""Extract antiSMASH zip output into flat TSVs the Shiny app can load
directly:
  1. data/{bacterium}_antismash_regions.tsv -- one row per predicted BGC
     region, with its predicted product type and (if antiSMASH found one)
     the most similar known cluster from the MIBiG database.
  2. data/{bacterium}_antismash_genes.tsv -- one row per gene INSIDE a
     detected region, with antiSMASH's own "gene_kind" classification
     (biosynthetic / biosynthetic-additional / transport / regulatory /
     other / resistance) -- the same categories antiSMASH's own region
     viewer color-codes gene arrows by. Used to draw a matching gggenes
     arrow diagram for a selected region in the app.

For each data/{bacterium}_antismash.zip (matched case-insensitively against
the app's known bacteria_filenames):
  1. Unzip to a scratch dir, locate the main antiSMASH results JSON
     (named "{something}.json", NOT one of the per-region .gbk files).
  2. For every contig record's "areas" (= BGC regions antiSMASH called),
     pull region number, coordinates, predicted product(s)/category, and
     the top knownclusterblast hit if any -> regions.tsv.
  3. For every CDS feature whose coordinates fall inside one of those
     areas, pull its qualifiers (locus_tag, gene name, product, gene_kind,
     strand) -> genes.tsv.
"""
import csv
import glob
import json
import os
import re
import shutil
import tempfile
import zipfile

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

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
lower_to_canonical = {b.lower(): b for b in bacteria_filenames}


def find_zips():
    out = {}
    # match "_antismash.zip" in any capitalisation (uploads arrive as e.g.
    # "_antiSmash.zip"), and allow an optional "_genome" before it
    for path in glob.glob(os.path.join(DATA_DIR, "*.zip")):
        base = os.path.basename(path)
        if not base.lower().endswith("_antismash.zip"):
            continue
        stem = base[:-len("_antismash.zip")]
        if stem.lower().endswith("_genome"):
            stem = stem[:-len("_genome")]
        canon = lower_to_canonical.get(stem.lower())
        if canon in out:
            print(f"  NOTE: two antiSMASH zips for {canon}; using {base}")
        if canon is None:
            print(f"  WARNING: {base} doesn't match any known bacterium filename (stem={stem!r}) -- skipped")
            continue
        out[canon] = path
    return out


def load_main_json(zip_path, bacterium):
    tmpdir = tempfile.mkdtemp(prefix="antismash_")
    try:
        with zipfile.ZipFile(zip_path) as zf:
            names = zf.namelist()
            # main results json: top-level *.json that is NOT under a subfolder
            # and isn't a per-region file (those are .gbk, not .json)
            candidates = [n for n in names if n.endswith(".json") and "/" not in n]
            if not candidates:
                candidates = [n for n in names if n.endswith(".json")]
            if not candidates:
                return None
            # prefer the one whose name looks like the bacterium/input filename
            chosen = candidates[0]
            for c in candidates:
                if bacterium.lower() in c.lower() or bacterium.lower().replace("_kb18", "_kb18") in c.lower():
                    chosen = c
                    break
            zf.extract(chosen, tmpdir)
            with open(os.path.join(tmpdir, chosen)) as fh:
                return json.load(fh)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


LOCATION_RE = re.compile(r"\[?<?(\d+):>?(\d+)\]?\(([+-])\)")


def parse_location(loc_str):
    """antiSMASH/Biopython location strings look like '[3913:4321](-)' (and
    sometimes with '<'/'>' fuzzy-boundary markers). Returns (start, end,
    strand) with start/end already 1-based-friendly (Biopython's `start` is
    0-based internally but antiSMASH's own JSON keeps that raw, so this is
    used only for containment checks against the equally-raw area
    start/end, not displayed to the user directly)."""
    m = LOCATION_RE.search(loc_str or "")
    if not m:
        return None, None, None
    return int(m.group(1)), int(m.group(2)), m.group(3)


def process(bacterium, zip_path):
    data = load_main_json(zip_path, bacterium)
    if data is None:
        print(f"  {bacterium}: no top-level results JSON found inside {zip_path} -- skipped")
        return

    region_rows = []
    gene_rows = []
    for rec in data.get("records", []):
        contig = rec.get("id", "")
        areas = rec.get("areas", []) or []
        cb = (rec.get("modules", {}) or {}).get("antismash.modules.clusterblast", {}) or {}
        known_results = ((cb.get("knowncluster") or {}).get("results")) or []
        known_by_region = {r["region_number"]: r for r in known_results}

        # pre-parse every CDS feature's coordinates once per record
        cds_feats = []
        for feat in rec.get("features", []) or []:
            if feat.get("type") != "CDS":
                continue
            fstart, fend, fstrand = parse_location(feat.get("location", ""))
            if fstart is None:
                continue
            cds_feats.append((fstart, fend, fstrand, feat.get("qualifiers", {}) or {}))

        for i, area in enumerate(areas, start=1):
            region_id = f"region{i:03d}"
            a_start, a_end = area.get("start", 0) or 0, area.get("end", 0) or 0
            products = area.get("products", []) or []
            categories = sorted({
                pc.get("category", "") for pc in (area.get("protoclusters", {}) or {}).values()
            })
            known = known_by_region.get(i)
            best_hit_accession = ""
            best_hit_desc = ""
            best_hit_similarity = ""
            if known and known.get("total_hits", 0) > 0 and known.get("ranking"):
                top_hit, top_score = known["ranking"][0]
                best_hit_accession = top_hit.get("accession", "")
                best_hit_desc = top_hit.get("description", "")
                best_hit_similarity = top_score.get("similarity", "")

            region_rows.append({
                "Region": region_id,
                "Contig": contig,
                "Start": a_start,
                "End": a_end,
                "Length_bp": a_end - a_start,
                "Predicted_Product": ";".join(products),
                "Category": ";".join(c for c in categories if c),
                "Known_Cluster_Accession": best_hit_accession,
                "Known_Cluster_Description": best_hit_desc,
                "Known_Cluster_Similarity_pct": best_hit_similarity,
            })

            # genes whose CDS falls inside this region's coordinates
            for fstart, fend, fstrand, q in cds_feats:
                if fstart >= a_start and fend <= a_end:
                    locus_tag = (q.get("locus_tag") or q.get("ID") or [""])[0]
                    gene_name = (q.get("gene") or [""])[0]
                    product = (q.get("product") or [""])[0]
                    gene_kind = (q.get("gene_kind") or ["other"])[0] or "other"
                    gene_functions = "; ".join(q.get("gene_functions", []) or [])
                    gene_rows.append({
                        "Region": region_id,
                        "Contig": contig,
                        "Locus_Tag": locus_tag,
                        "Gene_Name": gene_name,
                        "Product": product,
                        "Start": fstart,
                        "End": fend,
                        "Strand": fstrand,
                        "Gene_Kind": gene_kind,
                        "Gene_Functions": gene_functions,
                    })

    regions_path = os.path.join(DATA_DIR, f"{bacterium}_antismash_regions.tsv")
    with open(regions_path, "w", newline="") as fh:
        fieldnames = ["Region", "Contig", "Start", "End", "Length_bp", "Predicted_Product",
                      "Category", "Known_Cluster_Accession", "Known_Cluster_Description",
                      "Known_Cluster_Similarity_pct"]
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in region_rows:
            w.writerow(r)

    genes_path = os.path.join(DATA_DIR, f"{bacterium}_antismash_genes.tsv")
    with open(genes_path, "w", newline="") as fh:
        fieldnames = ["Region", "Contig", "Locus_Tag", "Gene_Name", "Product", "Start", "End",
                      "Strand", "Gene_Kind", "Gene_Functions"]
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in gene_rows:
            w.writerow(r)

    print(f"  {bacterium}: {len(region_rows)} BGC region(s) -> {os.path.basename(regions_path)}, "
          f"{len(gene_rows)} gene(s) in regions -> {os.path.basename(genes_path)}")


def main():
    zips = find_zips()
    if not zips:
        print("No *_antismash.zip files found in data/.")
        return
    print(f"Found {len(zips)} antiSMASH zip(s):")
    for bacterium, path in sorted(zips.items()):
        print(f"Processing {bacterium} <- {os.path.basename(path)}")
        process(bacterium, path)


if __name__ == "__main__":
    main()
