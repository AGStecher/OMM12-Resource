#!/usr/bin/env bash
# Run AMRFinderPlus for every OMM12 bacterium's Prokka output, one call per genome.
#
# Expects a sibling "prokka" folder laid out as:
#   ../prokka/{bacterium}/{bacterium}.faa
#   ../prokka/{bacterium}/{bacterium}.fna
#   ../prokka/{bacterium}/{bacterium}.gff
# (i.e. run this script from inside the folder that "prokka" is a sibling of --
#  adjust PROKKA_DIR below if your layout differs.)
#
# Writes output to amrfinder_out/{bacterium}_amrfinder.tsv -- copy those files into
# this app's data/ folder (as data/{bacterium}_amrfinder.tsv) and the AMR & Specialty
# Genes page will pick them up automatically, no other changes needed.

set -uo pipefail

PROKKA_DIR="../prokka"
OUT_DIR="amrfinder_out"

bacteria=(
  acutalibacter_muris_kb18
  akkermansia_muciniphila_YL44
  bacteroides_caecimuris_I48
  bifidobacterium_animalis_YL2
  blautia_coccoides_YL58
  clostridium_innocuum_I46
  enterocloster_clostridioformis_YL32
  enterococcus_faecalis_KB1
  flavonifractor_plautii_YL31
  limosilactobacillus_reuteri_I49
  muribaculum_intestinales_YL27
  turicimonas_muris_YL45
)

mkdir -p "$OUT_DIR"

n_ok=0
n_skip=0
n_fail=0

for b in "${bacteria[@]}"; do
  faa="${PROKKA_DIR}/${b}/${b}.faa"
  fna="${PROKKA_DIR}/${b}/${b}.fna"
  gff="${PROKKA_DIR}/${b}/${b}.gff"
  out="${OUT_DIR}/${b}_amrfinder.tsv"
  prot_out="${OUT_DIR}/${b}_protein_amrfinder.faa"

  if [[ ! -f "$faa" || ! -f "$fna" || ! -f "$gff" ]]; then
    echo "[SKIP] ${b} -- missing one of: ${faa} / ${fna} / ${gff}"
    n_skip=$((n_skip + 1))
    continue
  fi

  echo "[RUN ] ${b}"
  if amrfinder \
       -p "$faa" \
       -n "$fna" \
       -g "$gff" \
       --annotation_format prokka \
       --plus \
       -o "$out" \
       --protein_output "$prot_out"; then
    echo "[ OK ] ${b} -> ${out}"
    n_ok=$((n_ok + 1))
  else
    echo "[FAIL] ${b}"
    n_fail=$((n_fail + 1))
  fi
done

echo ""
echo "Done: ${n_ok} succeeded, ${n_skip} skipped (missing files), ${n_fail} failed."
echo "Copy ${OUT_DIR}/*_amrfinder.tsv into the app's data/ folder to load them in the resource."
