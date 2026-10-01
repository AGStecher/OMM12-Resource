#!/usr/bin/env python3
"""Flatten raw DeepLoc 2.0 output (data/{bacterium}_deeploc/results_*.csv)
into the data/{bacterium}_localization.csv the app actually reads.

Raw DeepLoc CSV columns: unnamed index, ACC, Localization, "Cell wall &
surface", Extracellular, Cytoplasmic, "Cytoplasmic Membrane",
"Outer Membrane", Periplasmic.

App-expected columns (see load_localization_csv() / the 8 existing
*_localization.csv files for the pattern): Protein_name, Localization,
Cell_wall_surface, Extracellular, Cytoplasmic, Cytoplasmic_Membrane,
Outer_Membrane, Periplasmic.

Safe to re-run: only (re)writes data/{bacterium}_localization.csv for
bacteria whose raw data/{bacterium}_deeploc/results_*.csv is present.
"""
import csv
import glob
import os
import sys

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

COL_MAP = {
    "ACC": "Protein_name",
    "Localization": "Localization",
    "Cell wall & surface": "Cell_wall_surface",
    "Extracellular": "Extracellular",
    "Cytoplasmic": "Cytoplasmic",
    "Cytoplasmic Membrane": "Cytoplasmic_Membrane",
    "Outer Membrane": "Outer_Membrane",
    "Periplasmic": "Periplasmic",
}
OUT_FIELDS = ["Protein_name", "Localization", "Cell_wall_surface", "Extracellular",
              "Cytoplasmic", "Cytoplasmic_Membrane", "Outer_Membrane", "Periplasmic"]


def find_raw(bacterium):
    folder = os.path.join(DATA_DIR, f"{bacterium}_deeploc")
    if not os.path.isdir(folder):
        return None
    candidates = sorted(glob.glob(os.path.join(folder, "results_*.csv")))
    # The DeepLocPro web server takes max 500 proteins per job, so a genome is
    # usually submitted in chunks -> merge ALL results_*.csv files in the folder.
    return candidates if candidates else None


def prokka_ids(bacterium):
    """Protein IDs (Prokka locus tags) of this strain, from data/{bacterium}_protein.faa."""
    faa = os.path.join(DATA_DIR, f"{bacterium}_protein.faa")
    if not os.path.exists(faa):
        return None
    with open(faa, encoding="utf-8") as fh:
        return {l[1:].split()[0] for l in fh if l.startswith(">")}


def convert(bacterium, raw_paths):
    out_path = os.path.join(DATA_DIR, f"{bacterium}_localization.csv")
    rows, seen = [], set()
    valid = prokka_ids(bacterium)
    for raw_path in raw_paths:
        # Skip result files from runs on other protein sets (e.g. NCBI WP_ IDs):
        # only proteins with this strain's Prokka IDs can be linked to the rest of the app.
        if valid is not None:
            with open(raw_path, newline="", encoding="utf-8") as fh:
                ids = [r.get("ACC") for r in csv.DictReader(fh)]
            if not any(i in valid for i in ids):
                print(f"  {bacterium}: {os.path.basename(raw_path)} has no Prokka protein IDs "
                      f"(older run on another protein set) -- ignored")
                continue
        with open(raw_path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            missing = [c for c in COL_MAP if c not in reader.fieldnames]
            if missing:
                print(f"  WARNING: {bacterium}: {os.path.basename(raw_path)} missing expected "
                      f"column(s) {missing} (found: {reader.fieldnames}) -- skipped")
                return 0
            for row in reader:
                if row["ACC"] in seen:  # same protein in two chunks -> keep first
                    continue
                if valid is not None and row["ACC"] not in valid:
                    continue
                seen.add(row["ACC"])
                rows.append({OUT_FIELDS[i]: row[src] for i, src in enumerate(COL_MAP)})
    if not rows:
        print(f"  {bacterium}: no usable proteins in the result files -- existing table left unchanged")
        return 0
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_FIELDS)
        w.writeheader()
        w.writerows(rows)
    return len(rows)


def main():
    force = "--force" in sys.argv
    for bacterium in bacteria_filenames:
        existing_csv = os.path.join(DATA_DIR, f"{bacterium}_localization.csv")
        raw = find_raw(bacterium)
        if raw is None:
            continue
        # rebuild when the raw DeepLocPro results are newer than the existing table
        if not force and os.path.exists(existing_csv) and \
                max(os.path.getmtime(f) for f in raw) <= os.path.getmtime(existing_csv):
            print(f"  {bacterium}: _localization.csv is up to date -- skipped")
            continue
        n = convert(bacterium, raw)
        if n:
            print(f"  {bacterium}: {n} proteins <- {len(raw)} result file(s) -> {bacterium}_localization.csv")


if __name__ == "__main__":
    main()
