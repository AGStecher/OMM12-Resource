#!/usr/bin/env python3
"""Extract CRISPRCasFinder zip output into flat TSVs the Shiny app can load
directly, mirroring the existing AMRFinderPlus/DefenseFinder pattern.

For each data/{bacterium}_crispercasfinder.zip (bacterium name matched
case-insensitively against the app's known bacteria_filenames, since
CRISPRCasFinder zips have been uploaded with inconsistent casing, e.g.
"acutalibacter_muris_KB18" vs the app's "acutalibacter_muris_kb18"):

  1. Unzip to a scratch dir, locate result.json (CRISPRCasFinder nests it
     under a session-hash subfolder).
  2. Flatten every CRISPR array across all sequences/contigs into
     data/{bacterium}_crispr_arrays.tsv
  3. Flatten every Cas gene (grouped by cluster) into
     data/{bacterium}_cas_systems.tsv

Safe to re-run: only (re)writes the two TSVs for bacteria whose zip is
present, never touches the zip itself or any other bacterium's files.
"""
import csv
import json
import os
import zipfile
import tempfile
import shutil
import glob

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
    # Uploaded CRISPRCasFinder zips have appeared under several misspelled
    # suffixes ("_crispercasfinder.zip", "_crsipercasfinder.zip",
    # "_crsipercastfinder.zip", ...). Rather than hard-coding every typo,
    # match any *.zip whose name (lowercased) contains both "cas" and
    # "find" and isn't an antiSMASH zip, then resolve the bacterium by
    # case-insensitive prefix match against the known filenames.
    out = {}
    for path in glob.glob(os.path.join(DATA_DIR, "*.zip")):
        base = os.path.basename(path)
        lower = base.lower()
        if "antismash" in lower:
            continue
        if "cas" not in lower or "find" not in lower:
            continue
        canon = None
        for b in bacteria_filenames:
            if lower.startswith(b.lower() + "_"):
                canon = b
                break
        if canon is None:
            print(f"  WARNING: {base} doesn't match any known bacterium filename -- skipped")
            continue
        out[canon] = path
    return out


def load_result_json(zip_path):
    tmpdir = tempfile.mkdtemp(prefix="ccf_")
    try:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tmpdir)
        matches = glob.glob(os.path.join(tmpdir, "**", "result.json"), recursive=True)
        if not matches:
            return None
        with open(matches[0]) as fh:
            return json.load(fh)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def process(bacterium, zip_path):
    data = load_result_json(zip_path)
    if data is None:
        print(f"  {bacterium}: no result.json found inside {zip_path} -- skipped")
        return

    array_rows = []
    cas_rows = []
    for seq in data.get("Sequences", []):
        seq_id = seq.get("Id", "")
        for arr in seq.get("Crisprs", []) or []:
            array_rows.append({
                "Sequence": seq_id,
                "CRISPR_Id": arr.get("Name", ""),
                "Start": arr.get("Start", ""),
                "End": arr.get("End", ""),
                "Length": (arr.get("End", 0) or 0) - (arr.get("Start", 0) or 0) + 1,
                "Orientation": arr.get("Potential_Orientation", ""),
                "DR_Consensus": arr.get("DR_Consensus", ""),
                "DR_Length": arr.get("DR_Length", ""),
                "Spacers_Nb": arr.get("Spacers", ""),
                "Evidence_Level": arr.get("Evidence_Level", ""),
                "Conservation_DRs_pct": arr.get("Conservation_DRs", ""),
                "Conservation_Spacers_pct": arr.get("Conservation_Spacers", ""),
            })
        for cas_cluster in seq.get("Cas", []) or []:
            cluster_type = cas_cluster.get("Type", "")
            cluster_start = cas_cluster.get("Start", "")
            cluster_end = cas_cluster.get("End", "")
            genes = cas_cluster.get("Genes", []) or []
            if not genes:
                cas_rows.append({
                    "Sequence": seq_id, "Cas_Cluster_Type": cluster_type,
                    "Cluster_Start": cluster_start, "Cluster_End": cluster_end,
                    "Gene_Subtype": "", "Gene_Start": "", "Gene_End": "", "Gene_Orientation": "",
                })
            for g in genes:
                cas_rows.append({
                    "Sequence": seq_id,
                    "Cas_Cluster_Type": cluster_type,
                    "Cluster_Start": cluster_start,
                    "Cluster_End": cluster_end,
                    "Gene_Subtype": g.get("Sub_type", ""),
                    "Gene_Start": g.get("Start", ""),
                    "Gene_End": g.get("End", ""),
                    "Gene_Orientation": g.get("Orientation", ""),
                })

    arr_path = os.path.join(DATA_DIR, f"{bacterium}_crispr_arrays.tsv")
    with open(arr_path, "w", newline="") as fh:
        fieldnames = ["Sequence", "CRISPR_Id", "Start", "End", "Length", "Orientation",
                      "DR_Consensus", "DR_Length", "Spacers_Nb", "Evidence_Level",
                      "Conservation_DRs_pct", "Conservation_Spacers_pct"]
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in array_rows:
            w.writerow(r)

    cas_path = os.path.join(DATA_DIR, f"{bacterium}_cas_systems.tsv")
    with open(cas_path, "w", newline="") as fh:
        fieldnames = ["Sequence", "Cas_Cluster_Type", "Cluster_Start", "Cluster_End",
                      "Gene_Subtype", "Gene_Start", "Gene_End", "Gene_Orientation"]
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in cas_rows:
            w.writerow(r)

    print(f"  {bacterium}: {len(array_rows)} CRISPR array row(s) -> {os.path.basename(arr_path)}, "
          f"{len(cas_rows)} Cas gene row(s) -> {os.path.basename(cas_path)}")


def main():
    zips = find_zips()
    if not zips:
        print("No *_crispercasfinder.zip files found in data/.")
        return
    print(f"Found {len(zips)} CRISPRCasFinder zip(s):")
    for bacterium, path in sorted(zips.items()):
        print(f"Processing {bacterium} <- {os.path.basename(path)}")
        process(bacterium, path)


if __name__ == "__main__":
    main()
