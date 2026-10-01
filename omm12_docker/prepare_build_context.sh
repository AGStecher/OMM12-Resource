#!/usr/bin/env bash
# Run this from inside omm12_docker/ (where this script lives), with
# OMM12_SRC pointing at your actual project folder (the one containing
# app.R, data/, www/) -- e.g. the OMM12_website folder itself, one level up.
#
# It copies app.R and www/ as-is, and copies data/ EXCLUDING the raw,
# never-read-by-the-app files: per-strain raw Prokka output folders
# (*_prokka/), raw COGclassifier output folders (*_cog/ -- the app only
# reads the flat *_cog.tsv file that sits alongside, not this folder),
# raw tool zips (antiSMASH/CRISPRCasFinder), raw DeepLoc folders
# (*_deeploc/), and the stray non-OMM12 E. coli files that ended up in
# data/ but aren't one of the 12 community members.
#
# This matters because data/ is currently ~854 MB, but the app itself
# only ever reads the flattened top-level files (e.g. *_prokka.gff,
# *_prokka.tsv, *_antismash_regions.tsv, *_localization.csv, etc.) --
# trimming this down keeps the Docker image small and avoids shipping
# raw intermediate files publicly for no reason. Verified: trims the
# data folder from ~854 MB down to ~250 MB, with all files app.R
# actually reads (genome_info.tsv, provenance.tsv, ortholog_best_hits.tsv,
# per-strain *_prokka.gff/.tsv/.gbk, *_cog.tsv, etc.) confirmed present.

set -euo pipefail

OMM12_SRC="${1:-..}"   # defaults to the parent folder (OMM12_website itself)

if [ ! -f "$OMM12_SRC/app.R" ]; then
  echo "Could not find app.R under: $OMM12_SRC"
  echo "Usage: ./prepare_build_context.sh /path/to/OMM12_website"
  exit 1
fi

# Wipe any previous trimmed copy first -- rsync's --delete-excluded does not
# reliably remove non-empty excluded directories left over from a prior run
# on all rsync versions, so start clean every time instead.
rm -rf app/data app/www
mkdir -p app/data app/www

cp "$OMM12_SRC/app.R" app/app.R
rsync -a "$OMM12_SRC/www/" app/www/

rsync -a \
  --exclude '*_prokka/' \
  --exclude '*_cog/' \
  --exclude '*.zip' \
  --exclude '*_deeploc/' \
  --exclude 'escherichia_coli_strain_Mt1B1*' \
  "$OMM12_SRC/data/" app/data/

echo "Done. Trimmed data size:"
du -sh app/data
echo "(compare to original:)"
du -sh "$OMM12_SRC/data"
