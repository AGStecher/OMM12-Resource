#!/usr/bin/env python3
"""Generate the taxonomy diagrams shown on each OMM12 member's page
(www/images/{bacterium_filename}.png).

All 12 images share one fixed canvas size and layout ("single-colour
ladder"): Phylum -> Class -> Order -> Family -> Genus -> Species -> Strain as
rounded boxes in one navy colour scale, light (broad rank) to dark (strain),
each box a little narrower than the one above; rank names on the left, genus
and species names in italics, and a source line at the bottom.

The lineage is read from data/genome_info.tsv (column taxonomic_classification,
taken from NCBI Taxonomy for the strain's GenBank genome record); species
names are NCBI's current names (see SPECIES below). Re-run after changing
genome_info.tsv:

    python build_crispr_bgc_pipeline/generate_taxonomy_diagrams.py

Needs Pillow (pip install pillow) and the Liberation Sans fonts; set FONT_DIR
if they live elsewhere on your system.
"""
import csv
import os
import sys

from PIL import Image, ImageDraw, ImageFont

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
INFO = os.path.join(ROOT, "data", "genome_info.tsv")
OUT_DIR = os.path.join(ROOT, "www", "images")
ACCESSED = "30 Sep 2026"   # date the lineage was checked on NCBI Taxonomy

FONT_DIR = os.environ.get("FONT_DIR", "/usr/share/fonts/truetype/liberation")
F_BOLD = os.path.join(FONT_DIR, "LiberationSans-Bold.ttf")
F_BOLD_IT = os.path.join(FONT_DIR, "LiberationSans-BoldItalic.ttf")
F_REG = os.path.join(FONT_DIR, "LiberationSans-Regular.ttf")

# Species names as they appear in NCBI Taxonomy (the strain's genome record)
SPECIES = {
    "acutalibacter_muris_kb18": ("Acutalibacter muris", "KB18"),
    "akkermansia_muciniphila_YL44": ("Akkermansia muciniphila", "YL44"),
    "bacteroides_caecimuris_I48": ("Bacteroides caecimuris", "I48"),
    "bifidobacterium_animalis_YL2": ("Bifidobacterium animalis", "YL2"),
    "blautia_coccoides_YL58": ("Blautia pseudococcoides", "YL58"),
    "clostridium_innocuum_I46": ("[Clostridium] innocuum", "I46"),
    "enterocloster_clostridioformis_YL32": ("Enterocloster clostridioformis", "YL32"),
    "enterococcus_faecalis_KB1": ("Enterococcus faecalis", "KB1"),
    "flavonifractor_plautii_YL31": ("Flavonifractor plautii", "YL31"),
    "limosilactobacillus_reuteri_I49": ("Limosilactobacillus reuteri", "I49"),
    "muribaculum_intestinales_YL27": ("Muribaculum intestinale", "YL27"),
    "turicimonas_muris_YL45": ("Turicimonas muris", "YL45"),
}

RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]
# light -> dark navy scale (same navy as the app header); text colour per box
COLORS = ["#dfe7f1", "#c6d4e6", "#a9bfd9", "#86a3c7", "#5f82b0", "#3d6196", "#1b263b"]
TEXT_COLORS = ["#1b263b"] * 4 + ["#ffffff"] * 3
ITALIC = {"Genus", "Species"}

# Fixed layout (pixels). Drawn at 2x for sharp display, same for every image.
W, H = 1000, 1420
LABEL_X = 40            # rank labels
BOX_X, BOX_W = 230, 700  # widest (phylum) box
INSET = 18               # each rank's box is 2*INSET px narrower than the one above
BOX_H, GAP = 150, 40
TOP = 40
RADIUS = 18
TEXT_MAX = 44           # max font size inside boxes
TEXT_MIN = 26
ONE_LINE_MIN = 38      # below this size a name is split over two lines instead
SPECIES_SIZE = 38      # species box: fixed size, two lines


def font(path, size):
    return ImageFont.truetype(path, size)


def fit_lines(draw, text, path, max_w):
    """Largest font size (TEXT_MAX..ONE_LINE_MIN) at which text fits in one line,
    else split into two lines (at the space nearest the middle)."""
    for size in range(TEXT_MAX, ONE_LINE_MIN - 1, -2):
        f = font(path, size)
        if draw.textlength(text, font=f) <= max_w:
            return [text], f
    words = text.split(" ")
    if len(words) > 1:
        best = min(range(1, len(words)),
                   key=lambda i: abs(len(" ".join(words[:i])) - len(" ".join(words[i:]))))
        lines = [" ".join(words[:best]), " ".join(words[best:])]
        for size in range(TEXT_MAX, TEXT_MIN - 1, -2):
            f = font(path, size)
            if all(draw.textlength(l, font=f) <= max_w for l in lines):
                return lines, f
    return [text], font(path, TEXT_MIN)


def arrow(draw, x, y0, y1):
    draw.line([(x, y0), (x, y1 - 18)], fill="#333333", width=8)
    draw.polygon([(x - 18, y1 - 22), (x + 18, y1 - 22), (x, y1)], fill="#333333")


def make(bacterium, lineage, taxid):
    species, strain = SPECIES[bacterium]
    ranks = [r.strip() for r in lineage.rstrip(".").split(";") if r.strip()]
    # lineage = Bacteria; <kingdom>; phylum; class; order; family; genus
    if len(ranks) < 7:
        sys.exit(f"{bacterium}: lineage has fewer ranks than expected: {lineage}")
    phylum, cls, order, family, genus = ranks[-5:]
    labels = [phylum, cls, order, family, genus, species, strain]

    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    f_rank = font(F_REG, 30)
    for i, (rank, text, col, tcol) in enumerate(zip(RANKS, labels, COLORS, TEXT_COLORS)):
        y = TOP + i * (BOX_H + GAP)
        x0, x1 = BOX_X + i * INSET, BOX_X + BOX_W - i * INSET
        d.rounded_rectangle([x0, y, x1, y + BOX_H], radius=RADIUS, fill=col)
        d.text((LABEL_X, y + BOX_H / 2), rank, font=f_rank, fill="#6b6b7b", anchor="lm")
        if rank == "Species":
            # always genus / epithet on two lines at one fixed size, so every strain looks the same
            g, _, ep = text.partition(" ")
            lines, f = [g, ep], font(F_BOLD_IT, SPECIES_SIZE)
        else:
            lines, f = fit_lines(d, text, F_BOLD_IT if rank in ITALIC else F_BOLD, (x1 - x0) - 50)
        lh = f.size * 1.2
        y0 = y + BOX_H / 2 - lh * (len(lines) - 1) / 2
        for k, line in enumerate(lines):
            d.text(((x0 + x1) / 2, y0 + k * lh), line, font=f, fill=tcol, anchor="mm")
    src = f"Lineage: NCBI Taxonomy (taxid {taxid}), accessed {ACCESSED}"
    d.text((W / 2, H - 40), src, font=font(F_REG, 26), fill="#6b6b7b", anchor="mm")
    img.save(os.path.join(OUT_DIR, f"{bacterium}.png"), optimize=True)


def main():
    with open(INFO, encoding="utf-8") as fh:
        rows = {r["bacterium_name_clean"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for b in SPECIES:
        r = rows[b]
        make(b, r["taxonomic_classification"], r["ncbi_taxonomy_id"])
        print("wrote", os.path.join(OUT_DIR, f"{b}.png"))


if __name__ == "__main__":
    main()
