#!/usr/bin/env python3
"""Turn eggNOG-mapper output into the two files the Shiny app reads:

  data/{bacterium}_eggnog.gff  -> eggNOG section (COG category, description,
                                  orthologous groups, Pfam, preferred name)
  data/{bacterium}_kegg.tsv    -> KEGG Pathways section (one row per gene x
                                  KEGG module, grouped by KEGG's own module
                                  categories)

Input: put the eggNOG-mapper result file for a strain here
         data/{bacterium}_eggnog/{anything}.emapper.annotations
(the plain-text .emapper.annotations file, not the Excel version).

KEGG module names/categories are downloaded once from rest.kegg.jp
(br:ko00002) and cached as data/kegg_modules.tsv -- needs internet the
first time only.

Usage (from the OMM12_website folder):
    python build_crispr_bgc_pipeline/eggnog_to_app.py           # skip strains that already have output
    python build_crispr_bgc_pipeline/eggnog_to_app.py --force   # overwrite existing output
"""
import csv
import glob
import json
import os
import sys
import urllib.request

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
MODULE_CACHE = os.path.join(DATA_DIR, "kegg_modules.tsv")

bacteria_filenames = [
    "acutalibacter_muris_kb18", "akkermansia_muciniphila_YL44",
    "bacteroides_caecimuris_I48", "bifidobacterium_animalis_YL2",
    "blautia_coccoides_YL58", "clostridium_innocuum_I46",
    "enterocloster_clostridioformis_YL32", "enterococcus_faecalis_KB1",
    "flavonifractor_plautii_YL31", "limosilactobacillus_reuteri_I49",
    "muribaculum_intestinales_YL27", "turicimonas_muris_YL45",
]

KEGG_COLS = ["locus_tag", "Gene_Symbol", "Description", "KO", "KO_module",
             "Pathway_1", "Pathway", "Main_Pathway", "Sub_pathway"]


def load_modules():
    """module id -> (name, level1, level2, level3)"""
    if not os.path.exists(MODULE_CACHE):
        print("Downloading KEGG module hierarchy (br:ko00002) from rest.kegg.jp ...")
        with urllib.request.urlopen("https://rest.kegg.jp/get/br:ko00002/json", timeout=60) as r:
            tree = json.load(r)
        rows = []

        def walk(node, path):
            if "children" not in node:
                name = node["name"]
                if name.startswith("M") and name[1:6].isdigit():
                    mid, rest = name[:6], name[6:].strip()
                    rest = rest.split(" [")[0].strip()  # drop trailing [PATH:..] / [RN:..]
                    rows.append([mid, rest] + (path[:3] + ["", "", ""])[:3])
                return
            for c in node["children"]:
                walk(c, path + [c["name"]])

        for c in tree["children"]:
            walk(c, [c["name"]])
        with open(MODULE_CACHE, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["module", "name", "level1", "level2", "level3"])
            w.writerows(rows)
        print(f"  cached {len(rows)} modules -> {MODULE_CACHE}")
    mods = {}
    with open(MODULE_CACHE, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mods[r["module"]] = (r["name"], r["level1"], r["level2"], r["level3"])
    return mods


def read_annotations(path):
    header, rows = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("#query"):
                header = line[1:].split("\t")
                continue
            if not line or line.startswith("#"):
                continue
            if header is None:
                sys.exit(f"{path}: no '#query' header line -- is this an .emapper.annotations file?")
            rows.append(dict(zip(header, line.split("\t"))))
    return rows


def clean(v):
    v = (v or "").strip()
    return "" if v == "-" else v


def gff_escape(v):
    return v.replace(";", ",").replace("=", ":").replace("\t", " ")


def write_gff(bacterium, rows):
    out = os.path.join(DATA_DIR, f"{bacterium}_eggnog.gff")
    with open(out, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("##gff-version 3\n## created from eggNOG-mapper annotations by eggnog_to_app.py\n")
        for r in rows:
            attrs = [f"ID={r['query']}"]
            for key, col in [("em_OGs", "eggNOG_OGs"), ("em_COG_cat", "COG_category"),
                             ("em_desc", "Description"), ("em_Preferred_name", "Preferred_name"),
                             ("em_PFAMs", "PFAMs"), ("em_score", "score"), ("em_evalue", "evalue")]:
                attrs.append(f"{key}={gff_escape(clean(r.get(col)))}")
            fh.write("\t".join([r["query"], "eggNOG-mapper", "CDS", "1", "1",
                                clean(r.get("score")) or ".", "+", ".", ";".join(attrs)]) + "\n")
    return out


def write_kegg(bacterium, rows, modules):
    out = os.path.join(DATA_DIR, f"{bacterium}_kegg.tsv")
    n = 0
    # The app reads this file as Latin-1, so write it in Latin-1.
    with open(out, "w", newline="", encoding="latin-1", errors="replace") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(KEGG_COLS)
        for r in rows:
            mods = [m for m in clean(r.get("KEGG_Module")).split(",") if m]
            if not mods:
                continue
            kos = [k.replace("ko:", "") for k in clean(r.get("KEGG_ko")).split(",") if k]
            for m in mods:
                name, l1, l2, l3 = modules.get(m, (m, "", "", ""))
                w.writerow([r["query"], clean(r.get("Preferred_name")), clean(r.get("Description")),
                            ",".join(kos), m, "; ".join(x for x in (l1, l2, l3) if x) or "Unclassified",
                            name, l2 or "Unclassified", l3 or "Unclassified"])
                n += 1
    return out, n


def main():
    force = "--force" in sys.argv
    modules = None
    found = False
    for b in bacteria_filenames:
        files = sorted(glob.glob(os.path.join(DATA_DIR, f"{b}_eggnog", "*.emapper.annotations")))
        if not files:
            continue
        found = True
        gff_out = os.path.join(DATA_DIR, f"{b}_eggnog.gff")
        kegg_out = os.path.join(DATA_DIR, f"{b}_kegg.tsv")
        if not force and (os.path.exists(gff_out) or os.path.exists(kegg_out)):
            print(f"  {b}: eggnog.gff or kegg.tsv already exists -- skipped (use --force to overwrite)")
            continue
        if modules is None:
            modules = load_modules()
        rows = read_annotations(files[-1])
        write_gff(b, rows)
        _, n = write_kegg(b, rows, modules)
        print(f"  {b}: {len(rows)} annotated proteins -> _eggnog.gff, {n} gene-module rows -> _kegg.tsv")
    if not found:
        print("No data/{bacterium}_eggnog/*.emapper.annotations files found.")


if __name__ == "__main__":
    main()
