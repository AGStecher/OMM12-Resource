# OMM12 Resource — Genomic & Functional Browser

An R Shiny web application for exploring genomic, functional, and phylogenetic data for the 12 bacterial strains of the **Oligo-Mouse-Microbiota-12 (OMM12)** synthetic gut community. This document describes every page and plot in the app, the tools/versions behind each dataset, and walks through one worked example (*Acutalibacter muris* KB18) so a new user can see exactly what each output means.

Last updated: 2026-09-14. The app is a single file, `app.R`, reading pre-computed data from `data/`. Where a dataset hasn't been generated for a given strain yet, the affected section says so explicitly rather than showing an empty chart.

---

## 1. What is OMM12?

OMM12 is a synthetic, defined 12-member bacterial community representing the major phyla of the mouse gut microbiota, designed by Brugiroux et al. (*Nature Microbiology*, 2016; [10.1038/nmicrobiol.2016.215](https://doi.org/10.1038/nmicrobiol.2016.215)). It's stable across mouse generations, reproducible across facilities, and confers colonization resistance against *Salmonella* Typhimurium — which has made it a widely used gnotobiotic model.

The 12 members (Prokka annotation strain codes in parentheses):

| # | Species | Strain |
|---|---------|--------|
| 1 | *Acutalibacter muris* | KB18 |
| 2 | *Akkermansia muciniphila* | YL44 |
| 3 | *Bacteroides caecimuris* | I48 |
| 4 | *Bifidobacterium animalis* | YL2 |
| 5 | *Blautia pseudococcoides* (type strain; called *B. coccoides* in Brugiroux et al. 2016) | YL58 |
| 6 | *[Clostridium] innocuum* (NCBI name; called *C. innocuum* in Brugiroux et al. 2016) | I46 |
| 7 | *Enterocloster clostridioformis* | YL32 |
| 8 | *Enterococcus faecalis* | KB1 |
| 9 | *Flavonifractor plautii* | YL31 |
| 10 | *Limosilactobacillus reuteri* | I49 |
| 11 | *Muribaculum intestinale* | YL27 |
| 12 | *Turicimonas muris* | YL45 |

---

## 2. How the app is organized

The sidebar has two kinds of pages:

- **Community-wide pages** — one page, all 12 members shown together (search, compare, phylogeny, defense/CRISPR/BGC/AMR overviews, pan-genome, primer design).
- **Bacteria Details** — pick one member from the OMM12 grid (or any "view details" link) to get a deep-dive page with ~13 sub-sections specific to that strain (genome viewer, annotation, COG, KEGG, eggNOG, MobileOG, localization, metabolic model, defense, CRISPR, BGC, AMR).

Every page also has a **Tools Guide** (sidebar, right after About) — a single reference listing what every tool identifies and what it's used for, and every section links back to it via a short one-line blurb at the top.

---

## 3. Tools & versions behind the data

All analyses were run **offline** by the lab and the flat output files dropped into `data/`; the app itself never re-runs any of these tools, it only reads and visualizes their output. Versions below are what's actually embedded in the output files themselves (not assumed) — where a tool doesn't record its version in its output, that's noted rather than guessed.

| Tool | Used for | Version (as recorded in the data) | Reference (paper / website) |
|---|---|---|---|
| **Prokka** | Genome annotation (CDS/tRNA/rRNA calls, the base gene catalog) | **1.13** (from KB18's `.log`; internally calls aragorn 1.2, barrnap 0.9, blastp 2.12) | Seemann 2014, *Bioinformatics* 30(14):2068–2069 — [doi.org/10.1093/bioinformatics/btu153](https://doi.org/10.1093/bioinformatics/btu153) · [github.com/tseemann/prokka](https://github.com/tseemann/prokka) |
| **COGclassifier** | COG functional category assignment | **2.0.0** (from each strain's `cogclassifier.log`) | [github.com/moshi4/COGclassifier](https://github.com/moshi4/COGclassifier) (no formal paper) |
| **eggNOG-mapper** | Deeper orthology annotation (KB18 only) | **2.1.12** (`emapper-2.1.12`, from the GFF header comment) | Cantalapiedra et al. 2021, *Mol Biol Evol* 38(12):5825–5829 — [doi.org/10.1093/molbev/msab293](https://doi.org/10.1093/molbev/msab293) · [eggnog-mapper.embl.de](http://eggnog-mapper.embl.de) |
| **KEGG** (Orthology/pathway reference) | KEGG Orthology (KO) + pathway/module mapping | Reference database, not versioned per-run | Kanehisa et al. 2023, *Nucleic Acids Res* 51(D1):D587–D592 · [genome.jp/kegg](https://www.genome.jp/kegg/) |
| **DeepLocPro** | Subcellular localization prediction | Not recorded in output | Moreno, Nielsen, Winther & Teufel 2024, *Bioinformatics* 40(12):btae677 — [doi.org/10.1093/bioinformatics/btae677](https://doi.org/10.1093/bioinformatics/btae677) · [services.healthtech.dtu.dk/services/DeepLocPro-1.0](https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/) |
| **mobileOG-db** (DIAMOND homology + keyword search) | Mobile genetic element detection | DB version not recorded | Brown et al. 2022, *Appl Environ Microbiol* — [doi.org/10.1128/aem.00991-22](https://doi.org/10.1128/aem.00991-22) · [mobileogdb.flsi.cloud.vt.edu](https://mobileogdb.flsi.cloud.vt.edu/) |
| **gapseq** | Genome-scale metabolic model reconstruction (all 12 strains; models from Zenodo [10.5281/zenodo.17358311](https://doi.org/10.5281/zenodo.17358311)) | **1.4.0** (from the SBML model notes) | Zimmermann, Kaleta & Waschina 2021, *Genome Biology* 22:81 — [doi.org/10.1186/s13059-021-02295-1](https://doi.org/10.1186/s13059-021-02295-1) · [github.com/jotech/gapseq](https://github.com/jotech/gapseq) |
| **AMRFinderPlus** | Curated antimicrobial resistance / virulence gene calls | Not recorded in output (reference DB accessions like `NF033117.2` are per-hit, not a tool version) | Feldgarden et al. 2021, *Sci Rep* 11:12728 — [doi.org/10.1038/s41598-021-91456-0](https://doi.org/10.1038/s41598-021-91456-0) · [github.com/ncbi/amr](https://github.com/ncbi/amr) |
| **DefenseFinder** | Anti-phage defense system detection | Not recorded in output | Tesson et al. 2022, *Nat Commun* 13:2561 — [doi.org/10.1038/s41467-022-30269-9](https://doi.org/10.1038/s41467-022-30269-9) · [defensefinder.mdmlab.fr](https://defensefinder.mdmlab.fr/) |
| **CRISPRCasFinder** | CRISPR array + Cas gene detection | **4.2.30** (from `result.json`) | Couvin et al. 2018, *Nucleic Acids Res* 46(W1):W246–W251 — [doi.org/10.1093/nar/gky425](https://doi.org/10.1093/nar/gky425) · [crisprcas.i2bc.paris-saclay.fr](https://crisprcas.i2bc.paris-saclay.fr) |
| **antiSMASH** | Biosynthetic gene cluster (BGC) detection | **8.0.4** (from the results JSON; schema version 4) | Blin et al. 2025, *Nucleic Acids Res* 53(W1):W32–W38 — [doi.org/10.1093/nar/gkaf334](https://doi.org/10.1093/nar/gkaf334) · [antismash.secondarymetabolites.org](https://antismash.secondarymetabolites.org) |
| **TYGS** | Whole-genome + 16S GBDP phylogenetic tree (KB18 only, per-bacterium "Phylogeny" tab) | Tree search via **FastME 2.1.4** (BioNJ starting tree + SPR postprocessing) | Meier-Kolthoff & Göker 2019, *Nat Commun* 10:2182 — [doi.org/10.1038/s41467-019-10210-3](https://doi.org/10.1038/s41467-019-10210-3) · [tygs.dsmz.de](https://tygs.dsmz.de) |
| **FastME** | Distance-based tree search algorithm used by TYGS | 2.1.4 | Lefort, Desper & Gascuel 2015, *Mol Biol Evol* 32(10):2798–2800 — [doi.org/10.1093/molbev/msv150](https://doi.org/10.1093/molbev/msv150) · [atgc-montpellier.fr/fastme](http://www.atgc-montpellier.fr/fastme/) |
| **Primer3** (methodology reference) | Primer Design tool is a **native R reimplementation** inspired by Primer3's approach — the real Primer3 binary is not called | Native reimplementation, not the original tool | Untergasser et al. 2012, *Nucleic Acids Res* 40(15):e115 — [doi.org/10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596) · [primer3.org](https://primer3.org) |
| *(this app's own pipeline)* | OMM12 Relatedness Tree (all 12 members, community-wide page) | From-scratch Python: Gotoh affine-gap pairwise alignment → center-star MSA → Jukes-Cantor (1969) distance → Neighbor-Joining → 200-replicate bootstrap. No external phylogenetics tools were available in the build environment, so this pipeline was written and unit-tested from scratch (see `build_omm12_16s_tree/`). | NJ method: Saitou & Nei 1987, *Mol Biol Evol* 4(4):406–425 — [doi.org/10.1093/oxfordjournals.molbev.a040454](https://doi.org/10.1093/oxfordjournals.molbev.a040454) |
| *(this app's own features)* | Cross-Genome Search, Sequence Search, Compare Bacteria, Gene Neighborhood, Pan-genome/Orthogroups | Native to this app — no external tool | — |

### Data completeness (as of 2026-09-14)

Not every tool has been run for every strain yet. The app tracks this explicitly per dataset (green "ran clean" vs red "not run yet" messaging) rather than silently showing zeros:

| Dataset | Strains with data |
|---|---|
| Prokka annotation | 12 / 12 |
| AMRFinderPlus | 12 / 12 |
| MobileOG | 12 / 12 |
| DefenseFinder | 12 / 12 |
| COG | 12 / 12 |
| CRISPR arrays | 12 / 12 (CRISPRCasFinder 4.2.30) |
| Localization (DeepLocPro) | 12 / 12 (KB18 and YL27 still use an older run that doesn't cover all proteins) |
| antiSMASH (BGC) | 12 / 12 (antiSMASH 8.0.4) |
| KEGG pathways | 2 / 12 (KB18, *B. caecimuris* I48) |
| eggNOG annotation | 1 / 12 (KB18 only) |
| Genome-scale metabolic model (GEM) | 12 / 12 (gapseq 1.4.0; Zimmermann & Burrichter 2025, Zenodo [10.5281/zenodo.17358311](https://doi.org/10.5281/zenodo.17358311)) |
| TYGS whole-genome / 16S GBDP tree | 1 / 12 (KB18 only) |

---

## 4. Community-wide pages

### About
Project description, OMM12 background, and a publication dashboard built from the screened OMM12 literature list: 334 research articles, reviews and book chapters (2015–2026) from a Dimensions.ai full-text search, each checked by hand for whether it actually uses OMM12 (133 do). Shows summary figures (publications, countries, Scopus citations, h-index, open access), publications and citations per year, a world map by author country (all / OMM12 users only) and the top 10 journals. The app reads only aggregated tables (`data/omm12_pubstats_*.tsv`); the per-publication list (`data/omm12_publications.tsv`) comes from the Dimensions free version and Scopus, whose terms don't allow redistributing it, so it is kept out of GitHub and the Docker image.

**To update it:** put new Dimensions exports, the updated `omm12_usage_check.csv` and `omm12_citations.csv` in `data/statistics_OMM12/figures_v7/OMM12_worldmap/`, then run `python build_publication_stats/build_publication_stats.py` (same screening rules and country matching as `OMM12_Publications_2015-2026.ipynb`; it also writes the aggregated `data/omm12_pubstats_*.tsv` via `make_publication_summary.py`) and rebuild the app. New publications that are not yet in `omm12_usage_check.csv` should be screened there first.

### Tools Guide
Reference page listing every tool/section, what it identifies, and what it's used for — grouped into Genome & Functional Annotation, Defense & Resistance, Biosynthesis, Comparative & Evolutionary, and Practical Tools.

### OMM12 Resource (home)
A 12-tile grid of all community members (taxonomy diagram picture + name); click any tile to open that strain's Bacteria Details page.

### Cross-Genome Search
Keyword/gene-symbol search (e.g. "bile salt hydrolase", "flagellin") across all 12 genomes' Prokka product descriptions at once. Returns a table of every matching gene, which bacterium it's in, and its product description.

### Compare Bacteria
Pick two members; get side-by-side genome-overview stats plus four charts: COG functional profile, subcellular localization profile, mobile genetic element categories, and a shared-genes (orthologs) panel.

**Example (KB18 vs. any other member):** the COG chart shows each bacterium's gene counts per COG letter category side by side as grouped bars — e.g. how many genes KB18 has annotated as category "J" (translation) vs. category "E" (amino acid metabolism), compared to the second strain.

### Gene Neighborhood / Synteny View
Pick a gene (by ortholog group or locus tag); see a gggenes arrow diagram of that gene and its flanking genes in its own genome, plus the same neighborhood in every other member that has an ortholog there, stacked for visual comparison of whether gene order (synteny) is conserved.

### Pan-genome / Ortholog Matrix
All 3,800+ genes across all 12 genomes clustered into orthologous groups (via reciprocal-best-hit sequence search, not just annotation-name matching), classified as:
- **Core** — present in all 12 members
- **Soft-core** — present in most (a configurable threshold below 12)
- **Shell** — present in a few
- **Unique** — present in only one member

A stacked bar chart shows each bacterium's genome composition by these four tiers. A browsable table lets you pick any orthogroup and see a gene tree (ggtree, built from a fast alignment-free distance) for the sequences in that group.

**Example:** KB18 contributes 3,818 CDS to the pan-genome; the orthogroup table shows, for any one of them, exactly which of the other 11 members share an ortholog and how many total members carry it.

### OMM12 Relatedness Tree
Two distinct visualizations, clearly separated because they measure different things:
1. **16S rRNA gene phylogeny** — a real molecular phylogeny built from this app's own from-scratch alignment/NJ/bootstrap pipeline (see §3 above), with bootstrap support shown on internal nodes. This is the one to treat as an actual evolutionary tree.
2. **Pairwise proteome similarity heatmap** — mean k-mer containment score between reciprocal-best-hit orthologs for every pair of the 12 genomes. This is a whole-proteome similarity signal, *not* a phylogeny (no alignment, no substitution model) — useful as a quick sanity check but not for inferring relationships.

### Sequence Search
Paste a raw or FASTA protein sequence; searches it against one or all 12 proteomes (BLAST-style) and returns matches with alignment stats. Use this when you have an actual sequence rather than a gene name.

### Primer Design
Pick a gene from any member; a native R Primer3-lite implementation designs candidate PCR primer pairs (forward/reverse, product size, melting temp, GC content) for that target sequence — a practical wet-lab tool for strain-specific detection or qPCR.

### Defense Systems Overview
Community-wide DefenseFinder results: bar charts of defense system types per bacterium and defense-vs-antidefense counts, plus a browsable table of every detected system (type, subtype, activity, genes involved) across all 12 members.

**Example:** KB18 has 14 defense systems across 23 genes detected by DefenseFinder, including anti-CRISPR (`acriia21`) hits.

### CRISPR Arrays & Cas Systems
Community-wide CRISPRCasFinder results: summary stats, and two browsable tables — CRISPR arrays (spacer count, repeat consensus, orientation, evidence level) and Cas gene clusters (subtype, gene positions), each filterable by bacterium.

**Example:** KB18 has exactly 1 CRISPR array (76 bp, 1 spacer, evidence level 1) and 0 Cas gene clusters detected — i.e. it carries CRISPR spacer memory but no complete adjacent Cas machinery was found by CRISPRCasFinder.

### Biosynthetic Gene Clusters (BGC)
Community-wide antiSMASH results, shown as antiSMASH's own native visualization style rather than a generic table: pick a bacterium and one of its predicted regions to see a gggenes arrow diagram of that region's genes, colored by antiSMASH's own gene-kind classification (core biosynthetic / additional biosynthetic / transport / regulatory / resistance / other), plus the region's predicted product/category and any MIBiG known-cluster match.

**Example:** KB18 has 7 predicted BGC regions covering 151 genes total; region 1 can be selected to see its exact gene-by-gene arrow diagram and, if antiSMASH found a match, the closest known reference cluster from the MIBiG database with a similarity percentage.

### AMR & Specialty Genes
Two complementary views, stacked so the more reliable one is checked first:
1. **Curated AMRFinderPlus calls** — real % identity / % coverage / matched NCBI reference accession per hit, run offline per strain. This takes priority wherever available.
2. **Keyword screen** — a fallback scan of Prokka product descriptions for resistance/virulence-related keywords, used community-wide as an approximate cross-check, clearly labeled as less reliable than the curated calls.

**Example:** KB18's curated AMRFinderPlus table shows one hit — `vanR` (VanR-ABDEGLN family response regulator), core scope, AMR type, glycopeptide class / vancomycin subclass, 100% coverage, 73.19% identity to the NCBI reference `WP_063856721.1`.

---

## 5. Bacteria Details (per-strain deep dive)

Selecting any member opens a page with a General Information box (taxonomy, genome size, GC%, shape, lifestyle, NCBI links) and a Data Provenance box (which tool version produced which dataset, and when), followed by these sub-sections (navigable via the left content index):

| Sub-section | What it shows |
|---|---|
| **Genome Viewer** | Linear (Prokka GFF-based) and circular (CGView) genome maps. |
| **Phylogenetic Tree** | The TYGS whole-genome/16S GBDP tree — **KB18 only**; other strains show a note explaining none has been built yet and how to add one. |
| **Annotation (Prokka)** | Full Prokka feature table (locus tag, gene name, product, coordinates) plus the raw Prokka summary stats (contig count, bases, tRNA/tmRNA/rRNA/CDS/gene counts). |
| **COG** | Bar chart of gene counts per COG functional category, with the category legend. |
| **eggNOG** | Supplemental annotation (ortholog groups, Pfam domains, preferred gene names) — **KB18 only** so far. |
| **KEGG Pathways** | Collapsible accordion (Main pathway → Sub-pathway → gene table) of KEGG Orthology assignments. Data present for 2/12 strains. |
| **MobileOG** | Bar chart of mobile genetic element categories detected (DIAMOND homology + keyword search against mobileOG-db). |
| **Localization** | Bar chart of predicted subcellular localization (DeepLocPro's 6 classes), plus a PCA scatter of each protein's probability vector (see caveat below), a browsable gene table filterable by clicking a bar, and a localization × COG cross-tab. |
| **Genome-Scale Metabolic Model** | gapseq-reconstructed metabolic network — **KB18 only**. Subsystem/pathway overview bar chart (click a bar to filter), a reactions table with MetaCyc links, and a force-directed metabolite↔reaction network graph for the selected pathway. |
| **Defense Systems** | Same DefenseFinder data as the community overview, filtered to this strain. |
| **CRISPR & Cas Systems** | Same CRISPRCasFinder data, filtered to this strain. |
| **Biosynthetic Gene Clusters** | Same antiSMASH region-viewer, filtered to this strain. |
| **AMR & Specialty Genes** | Same curated + keyword AMR data, filtered to this strain. |

### A note on the Localization PCA plot
The "Protein Projection by Predicted Localization" scatter is a PCA of each protein's 6-class DeepLocPro probability vector, colored by the predicted (argmax) class. It's built from the model's own output probabilities, not independent sequence embeddings — so clusters mostly restate the predicted class rather than validating it. It's useful for spotting **ambiguous** calls (points sitting between two color clusters have split probability mass) and for seeing the overall class mix at a glance, but it isn't independent evidence for the localization calls themselves.

---

## 6. Worked example: *Acutalibacter muris* KB18

KB18 is the OMM12 reference strain and has the most complete dataset, so it's the best strain to explore every feature on:

- **Genome:** 3,802,913 bp, 1 contig, 54.6% GC, NCBI assembly `GCF_016697365.1` (GenBank CP065321).
- **Genes (Prokka 1.13):** 3,818 CDS, 54 tRNA, 6 rRNA, 1 tmRNA, 3,879 total genes.
- **Defense Systems:** 14 systems across 23 genes, including anti-CRISPR system `acriia21`.
- **CRISPR:** 1 array (76 bp repeat, 1 spacer); 0 Cas gene clusters detected.
- **Biosynthetic Gene Clusters:** 7 predicted regions, 151 genes total.
- **AMR:** 1 curated hit (`vanR`, vancomycin-resistance-associated regulator, 73.19% identity to reference).
- **Localization:** 2,754 proteins classified by DeepLocPro.
- **MobileOG:** 230 mobile genetic element hits.
- **eggNOG / GEM / TYGS tree:** the only strain with all three — this is why KB18 is the one strain where the eggNOG annotation box, the Genome-Scale Metabolic Model sub-section, and a real whole-genome Phylogenetic Tree all actually render content instead of a "not run yet" message.

To see this yourself: open **Bacteria Details → Acutalibacter muris KB18**, then step through the content index on the left — every sub-section will have real data, which makes it the fastest way to see what a "fully populated" strain page looks like before checking a partially-populated one.

---

## 7. Data pipeline (how raw tool output becomes what the app shows)

The app never runs any bioinformatics tool itself — it only reads pre-computed flat files from `data/`. Several tools' native output isn't directly usable (nested JSON, zipped result bundles, inconsistent per-run filenames), so small Python extraction scripts flatten them into the TSV/CSV format the app's R loaders expect:

- `build_crispr_bgc_pipeline/extract_crisprcasfinder.py` — unzips CRISPRCasFinder output, flattens `Crisprs_REPORT` + `Cas_REPORT` into `{bacterium}_crispr_arrays.tsv` / `{bacterium}_cas_systems.tsv`.
- `build_crispr_bgc_pipeline/extract_antismash.py` — unzips antiSMASH output, extracts region-level (`{bacterium}_antismash_regions.tsv`) and per-gene (`{bacterium}_antismash_genes.tsv`) data from the results JSON.
- `build_crispr_bgc_pipeline/extract_deeploc.py` — flattens raw DeepLoc output (`{bacterium}_deeploc/results_*.csv`) into `{bacterium}_localization.csv`.
- `build_crispr_bgc_pipeline/eggnog_to_app.py` — turns eggNOG-mapper `.emapper.annotations` (in `data/{bacterium}_eggnog/`) into `{bacterium}_eggnog.gff` and a KEGG-module table `{bacterium}_kegg.tsv` (module names/categories from KEGG br:ko00002, cached in `data/kegg_modules.tsv`).
- `preparse_models.R` — converts every `data/{bacterium}_model.xml` (gapseq) into the `.rds` the app reads.
- `to_run/HOW_TO_RUN.md` — step-by-step instructions and ready-made input files for the analyses still missing.
- `build_omm12_16s_tree/` — the full from-scratch 16S phylogenetics pipeline (extraction → alignment → distance → NJ → bootstrap) behind the OMM12 Relatedness Tree's real phylogeny.

All extraction scripts are idempotent and case-insensitive about strain-name matching in uploaded filenames (uploaded zips have arrived with inconsistent casing and even misspelled suffixes, e.g. `_crsipercasfinder.zip`) — safe to re-run any time new raw output is dropped into `data/`, they only touch the strains whose raw file is present.

---

## 8. File naming convention

Every per-strain data file follows `data/{bacterium_filename}_{suffix}`, where `bacterium_filename` is the lowercase-with-strain-code identifier used throughout the app (e.g. `enterocloster_clostridioformis_YL32`), matching the `www/images/{bacterium_filename}.png` picture and the `bacteria_filenames` list at the top of `app.R`. Community-wide aggregation functions in `app.R` (`build_all_*_data()`) scan for these files by that convention and automatically pick up new strains' data on next app restart — no code changes needed when new tool output is dropped in, as long as the filename suffix matches what that section's loader expects (see the Tools Guide in-app, or the loader functions near the top of `app.R`, for the exact expected suffix per dataset).
