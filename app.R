# app.R - Complete OMM12 Resource with Enhanced Phylogenetic Tree

#search button for a specific gene

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
library(ggplot2)
library(plotly)
library(dplyr)
library(readr)
library(ape)
library(ggtree)
library(DT)
library(stringr)
library(tidyr)
library(wordcloud2)
library(ggiraph)
library(ggtreeExtra)
library(shinySearchbar)
library(shinyWidgets)
library(networkD3)
library(xml2)
# top of app.R
library(shinycssloaders)
# Gene Neighborhood / Synteny view: proper gene-arrow diagrams
# (install.packages("gggenes") if not already present)
library(gggenes)

# CGView.js runs client-side in the browser and fetches the .gbk file itself
# via a plain HTTP request, so it needs a URL it can actually reach. Shiny
# only serves the www/ folder by default -- data/ (where the real _prokka.gbk
# files live for all bacteria) is server-side only. This exposes it at the
# "data/" URL path without duplicating every .gbk file into www/.
addResourcePath("data", "data")

# Small null-coalescing helper (input$x is NULL before the user touches a
# control, or briefly during re-render) -- used by the Primer Design inputs.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Case-insensitive literal text match. User search terms and gene product names
# contain brackets/plus signs (e.g. "Na(+)/H(+) antiporter"), which break
# grepl() when treated as regular expressions.
contains_ci <- function(x, pattern) grepl(tolower(pattern), tolower(x), fixed = TRUE)

# Undo GFF3 percent-encoding in attribute values (%2C = comma, %3B = ;, %3D = =, %26 = &, %25 = %)
gff_unescape <- function(x) {
  x <- gsub("%2C", ",", x, fixed = TRUE); x <- gsub("%3B", ";", x, fixed = TRUE)
  x <- gsub("%3D", "=", x, fixed = TRUE); x <- gsub("%26", "&", x, fixed = TRUE)
  gsub("%25", "%", x, fixed = TRUE)
}

# --- Define the names of your 12 bacteria here ---
bacteria_names <- c(
  "Acutalibacter muris KB18",
  "Akkermansia muciniphila YL44",
  "Bacteroides caecimuris I48",
  "Bifidobacterium animalis YL2",
  "Blautia pseudococcoides YL58",   # type strain of B. pseudococcoides (LPSN); called B. coccoides in Brugiroux et al. 2016
  "[Clostridium] innocuum I46",   # NCBI name: brackets = not a true Clostridium, not yet reclassified
  "Enterocloster clostridioformis YL32",
  "Enterococcus faecalis KB1",
  "Flavonifractor plautii YL31",
  "Limosilactobacillus reuteri I49",
  "Muribaculum intestinale YL27",
  "Turicimonas muris YL45"
)


# Filenames on disk keep the strain code in uppercase (e.g. "..._YL44", "..._I48"),
# matching how Prokka/the data files are actually named. A plain tolower() of
# bacteria_names doesn't match that, so it's spelled out explicitly here.
bacteria_filenames <- c(
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
  "turicimonas_muris_YL45"
)

# --- DEFINE YOUR ICONS AND DATA TYPES ---
icon_paths <- list(
  "genome" = "icons/genome_icon1.png",
  "phylogeny" = "icons/phylogeny_icon.png",
  "cog" = "icons/cog_icon.png",
  "metabolic" = "icons/metabolic_icon.png",
  "annotation" = "icons/annotation_icon.png",
  "mobileog" = "icons/mobileog_icon.png",
  "localization" = "icons/localization_icon.png"
)

data_types <- c("genome", "phylogeny", "cog", "metabolic", "annotation","mobileog","localization") 

# --- DATA LOADING FUNCTIONS ---

# Function to load Prokka GFF file and extract annotation information
load_prokka_gff <- function(bacterium_filename) {
  # Prefer the COG-enriched GFF; fall back to the plain Prokka GFF if the
  # COG-merge step hasn't been done yet for this bacterium (COG_category will
  # just come back NA/"Unknown" for every gene in that case, handled downstream).
  gff_file <- file.path("data", paste0(bacterium_filename, "_prokka_1.gff"))
  plain_file <- file.path("data", paste0(bacterium_filename, "_prokka.gff"))
  used_enriched <- file.exists(gff_file)
  if (!used_enriched) {
    if (file.exists(plain_file)) {
      message("COG-enriched GFF not found, using plain Prokka GFF instead: ", plain_file)
      gff_file <- plain_file
    } else {
      message("GFF file not found: ", gff_file)
      return(NULL)
    }
  }
  gff_data <- tryCatch({
    read_tsv(gff_file, comment = "#", col_names = FALSE, show_col_types = FALSE) %>%
      select(seqid = 1, source = 2, type = 3, start = 4, end = 5,
             score = 6, strand = 7, phase = 8, attributes = 9) %>%
      mutate(
        gene = str_extract(attributes, "gene=([^;]+)") %>% str_remove("gene="),
        locus_tag = str_extract(attributes, "locus_tag=([^;]+)") %>% str_remove("locus_tag="),
        # GFF3 percent-encodes commas etc. in attribute values (e.g. "delta(4%2C6)")
        product = str_extract(attributes, "product=([^;]+)") %>% str_remove("product=") %>% gff_unescape(),
        gene = ifelse(is.na(gene), locus_tag, gene),
        COG_category = str_extract(attributes, "COG_category=([^;]+)") %>% str_remove("COG_category=")
      )
  }, error = function(e) { message("Error reading GFF file: ", e$message); return(NULL) })

  # Guard: on at least one bacterium (acutalibacter_muris_kb18) the COG-merge
  # step that produces the "_prokka_1.gff" clobbers column 1 (seqid, meant to
  # be the genome's contig/chromosome accession) with each feature's own
  # locus_tag instead -- every gene ends up looking like it's on its own
  # private "contig", which breaks anything keyed on Contig (nucleotide
  # extraction for Primer Design, Gene Neighborhood/synteny grouping). Detect
  # that pattern (seqid == locus_tag for most rows) and repair seqid using the
  # plain Prokka GFF's (correct) seqid column, joined on locus_tag.
  if (used_enriched && !is.null(gff_data) && nrow(gff_data) > 0 && file.exists(plain_file)) {
    corrupted_frac <- mean(gff_data$seqid == gff_data$locus_tag, na.rm = TRUE)
    if (!is.na(corrupted_frac) && corrupted_frac > 0.5) {
      message("Corrupted seqid column detected in ", gff_file, " -- repairing from ", plain_file)
      plain_seqid <- tryCatch({
        read_tsv(plain_file, comment = "#", col_names = FALSE, show_col_types = FALSE) %>%
          select(seqid = 1, type = 3, attributes = 9) %>%
          mutate(locus_tag = str_extract(attributes, "locus_tag=([^;]+)") %>% str_remove("locus_tag=")) %>%
          filter(!is.na(locus_tag), seqid != locus_tag) %>%
          distinct(locus_tag, .keep_all = TRUE) %>%
          select(locus_tag, real_seqid = seqid)
      }, error = function(e) NULL)
      if (!is.null(plain_seqid) && nrow(plain_seqid) > 0) {
        gff_data <- gff_data %>%
          left_join(plain_seqid, by = "locus_tag") %>%
          mutate(seqid = ifelse(!is.na(real_seqid), real_seqid, seqid)) %>%
          select(-real_seqid)
      }
    }
  }
  return(gff_data)
}

# Function to load Subcellular Localization CSV file
load_localization_csv <- function(bacterium_filename) {
  csv_file <- file.path("data", paste0(bacterium_filename, "_localization.csv"))
  if (!file.exists(csv_file)) { return(NULL) }
  localization_data <- tryCatch({
    read_csv(csv_file, show_col_types = FALSE)
  }, error = function(e) { message("Error reading Localization CSV file: ", e$message); return(NULL) })
  if (is.null(localization_data) || !("Localization" %in% colnames(localization_data))) { return(NULL) }
  return(localization_data)
}

# Build a schematic bacterial-cell cross-section (SVG) with DeepLocPro
# protein counts placed at each compartment. Covers all 6 DeepLocPro classes;
# outer membrane / periplasm only apply to Gram-negative envelopes and cell
# wall & surface mainly to Gram-positive ones, so both are shown with a caption
# rather than guessing the Gram status per bacterium (not reliably in our data).
build_localization_cell_diagram <- function(counts_named) {
  get_n <- function(key) {
    v <- counts_named[[key]]
    if (is.null(v) || is.na(v)) 0 else v
  }
  cyto      <- get_n("Cytoplasmic")
  cyto_mem  <- get_n("Cytoplasmic Membrane")
  periplasm <- get_n("Periplasmic")
  cell_wall <- get_n("Cell wall & surface")
  outer_mem <- get_n("Outer Membrane")
  extra     <- get_n("Extracellular")

  cx <- 190; cy <- 190
  legend_row <- function(y, color, label, n) {
    sprintf(
      '<rect x="380" y="%d" width="16" height="16" rx="3" fill="%s" stroke="#0d1b2a" stroke-width="0.5"/>
       <text x="404" y="%d" font-size="13" fill="#1b263b" font-family="sans-serif">%s: <tspan font-weight="700">%d</tspan></text>',
      y, color, y + 13, label, n
    )
  }

  svg <- sprintf('
  <svg viewBox="0 0 620 380" xmlns="http://www.w3.org/2000/svg" style="width:100%%; max-width:620px; height:auto;">
    <rect x="0" y="0" width="620" height="380" fill="#f4f6f9"/>
    <text x="%d" y="20" text-anchor="middle" font-size="13" fill="#666" font-family="sans-serif">Extracellular space</text>

    <circle cx="%d" cy="%d" r="165" fill="none" stroke="#c3ccd9" stroke-width="1" stroke-dasharray="4,3"/>
    <circle cx="%d" cy="%d" r="150" fill="none" stroke="#2c3e57" stroke-width="10" stroke-dasharray="6,4" opacity="0.85"/>
    <circle cx="%d" cy="%d" r="128" fill="none" stroke="#e3c766" stroke-width="14" opacity="0.9"/>
    <circle cx="%d" cy="%d" r="102" fill="none" stroke="#e8edf3" stroke-width="16"/>
    <circle cx="%d" cy="%d" r="78"  fill="none" stroke="#c9a227" stroke-width="14"/>
    <circle cx="%d" cy="%d" r="55"  fill="#0d1b2a"/>
    <text x="%d" y="%d" text-anchor="middle" font-size="12" fill="#ffffff" font-family="sans-serif" font-weight="600">Cytoplasm</text>
    <text x="%d" y="%d" text-anchor="middle" font-size="18" fill="#ffffff" font-family="sans-serif" font-weight="700">%d</text>

    %s
  </svg>',
    cx, cx, cy, cx, cy, cx, cy, cx, cy, cx, cy, cx, cy,
    cx, cy - 4, cx, cy + 16, cyto,
    paste(
      legend_row(30,  "#0d1b2a", "Cytoplasm", cyto),
      legend_row(56,  "#c9a227", "Cytoplasmic membrane", cyto_mem),
      legend_row(82,  "#e8edf3", "Periplasmic", periplasm),
      legend_row(108, "#e3c766", "Cell wall &amp; surface", cell_wall),
      legend_row(134, "#2c3e57", "Outer membrane", outer_mem),
      legend_row(160, "#f4f6f9", "Extracellular", extra),
      sep = "\n"
    )
  )
  svg
}

# Function to load Prokka TSV file (summary of annotations)
load_prokka_tsv <- function(bacterium_filename) {
  tsv_file <- file.path("data", paste0(bacterium_filename, "_prokka.tsv"))
  if (!file.exists(tsv_file)) { return(NULL) }
  tsv_data <- tryCatch({
    read_tsv(tsv_file, show_col_types = FALSE)
  }, error = function(e) { message("Error reading TSV file: ", e$message); return(NULL) })
  # Most strains' Prokka .tsv files contain an extra "gene" row (no product)
  # next to every CDS/tRNA/rRNA row; drop it so each gene is listed once.
  if (!is.null(tsv_data) && "ftype" %in% colnames(tsv_data)) {
    tsv_data <- tsv_data[is.na(tsv_data$ftype) | tsv_data$ftype != "gene", ]
  }
  return(tsv_data)
}

# Function to load COG CSV
load_cog_tsv <- function(bacterium_filename) {
  tsv_file <- file.path("data", paste0(bacterium_filename, "_cog.tsv"))
  if (!file.exists(tsv_file)) { 
    message("COG TSV file not found: ", tsv_file)
    if (file.exists(file.path("data", basename(tsv_file)))) {
      tsv_file <- file.path("data", basename(tsv_file))
    } else {
      return(NULL)
    }
  }
  
  cog_data <- tryCatch({
    read_tsv(tsv_file, show_col_types = FALSE)
  }, error = function(e) {
    message("Error reading COG TSV file: ", e$message)
    return(NULL)
  })
  
  if (!("COG_category" %in% colnames(cog_data))) {
    message("COG_category column missing in ", tsv_file)
    return(NULL)
  }
  
  return(cog_data)
}

load_feature_table <- function(bacterium_filename) {
  f <- file.path("data", paste0(bacterium_filename, "_feature_table.txt"))
  if (!file.exists(f)) { message("Feature table not found: ", f); return(NULL) }
  
  d <- tryCatch(
    read_tsv(f, show_col_types = FALSE, na = c("", "NA")),
    error = function(e) { message("Error reading feature table: ", e$message); NULL }
  )
  if (is.null(d)) return(NULL)
  
  # NCBI prefixes the first header with "# " -> column is named "# feature".
  names(d) <- trimws(sub("^#\\s*", "", names(d)))   # "# feature" -> "feature"
  
  if (!all(c("feature", "product_accession") %in% colnames(d))) {
    message("feature table missing expected columns"); return(NULL)
  }
  
  d %>%
    filter(feature == "CDS", !is.na(product_accession)) %>%
    transmute(
      product_accession = trimws(product_accession),
      ft_locus_tag = locus_tag,
      ft_gene      = symbol,
      ft_start     = start,
      ft_gene_name = Gene_name,
      ft_end       = end,
      ft_strand    = strand
    ) %>%
    distinct(product_accession, .keep_all = TRUE)
}

# Function to calculate annotation statistics
get_annotation_stats <- function(gff_data) {
  if (is.null(gff_data) || nrow(gff_data) == 0) {
    return(data.frame(type = "No data", Count = 0, stringsAsFactors = FALSE))
  }
  if (nrow(gff_data) == 0) {
    return(data.frame(type = "No data", Count = 0, stringsAsFactors = FALSE))
  }
  # Check if 'type' column exists
  if (!"type" %in% colnames(gff_data)) {
    return(data.frame(type = "No data", Count = 0, stringsAsFactors = FALSE))
  }
  stats <- gff_data %>%
    filter(!is.na(type), type != "", type != "gene") %>%
    group_by(type) %>%
    summarise(Count = n(), .groups = "drop") %>%
    arrange(desc(Count))
  if (nrow(stats) == 0) {
    stats <- data.frame(type = "No data", Count = 0, stringsAsFactors = FALSE)
  }
  return(stats)
}

# Function to load MobileOG CSV file
load_mobileog_csv <- function(bacterium_filename) {
  csv_file <- file.path("data", paste0(bacterium_filename, "_mobileog.csv"))
  if (!file.exists(csv_file)) { return(NULL) }
  mobileog_data <- tryCatch({
    read_csv(csv_file, show_col_types = FALSE)
  }, error = function(e) { message("Error reading MobileOG CSV file: ", e$message); return(NULL) })
  if (is.null(mobileog_data) || !("Major mobileOG Category" %in% colnames(mobileog_data))) { return(NULL) }
  return(mobileog_data)
}

# Function to load MobileOG CSV file
load_kegg_tsv <- function(bacterium_filename) {
  # Build file path
  tsv_file <- file.path("data", paste0(bacterium_filename, "_kegg.tsv"))
  
  # Return NULL if file does not exist
  if (!file.exists(tsv_file)) return(NULL)
  
  # Read TSV with data.table::fread and convert to UTF-8
  kegg_data <- tryCatch({
    df <- data.table::fread(tsv_file, sep = "\t", header = TRUE, encoding = "Latin-1")
    df[] <- lapply(df, function(x) iconv(x, from = "Latin1", to = "UTF-8", sub = ""))
    as.data.frame(df)          # <-- add this
  }, error = function(e) {
    message("Error reading KEGG TSV file: ", e$message)
    return(NULL)
  })
  
  # Check required columns
  required_cols <- c("Main_Pathway", "Sub_pathway", "Pathway")
  if (!all(required_cols %in% colnames(kegg_data))) {
    message("KEGG TSV file is missing required columns: ", paste(required_cols, collapse = ", "))
    return(NULL)
  }
  
  return(kegg_data)
}

# Defense systems (DefenseFinder) — systems summary
# Accepts both "{bacterium}_defense_systems.tsv" (this app's expected name) and
# "{bacterium}_defense_finder_systems.tsv" (DefenseFinder's own default output
# name from `defense-finder run`, e.g. for enterocloster_clostridioformis_YL32).
load_defense_systems <- function(bacterium_filename) {
  candidates <- file.path("data", paste0(bacterium_filename,
                           c("_defense_systems.tsv", "_defense_finder_systems.tsv")))
  f <- candidates[file.exists(candidates)][1]
  if (is.na(f) || is.null(f)) return(NULL)
  d <- tryCatch(read_tsv(f, show_col_types = FALSE),
                error = function(e) { message("Error reading defense systems: ", e$message); NULL })
  if (is.null(d) || !("type" %in% colnames(d))) return(NULL)
  d
}

# Defense genes (DefenseFinder) — gene-level detail
# Same dual-naming fallback as load_defense_systems() above.
load_defense_genes <- function(bacterium_filename) {
  candidates <- file.path("data", paste0(bacterium_filename,
                           c("_defense_genes.tsv", "_defense_finder_genes.tsv")))
  f <- candidates[file.exists(candidates)][1]
  if (is.na(f) || is.null(f)) return(NULL)
  d <- tryCatch(read_tsv(f, show_col_types = FALSE),
                error = function(e) { message("Error reading defense genes: ", e$message); NULL })
  d
}

# CRISPR arrays (CRISPRCasFinder) — flattened from result.json by
# build_crispr_bgc_pipeline/extract_crisprcasfinder.py (raw CRISPRCasFinder
# output is a zip bundle, not something this app parses live). Evidence_Level
# is CRISPRCasFinder's own 1-4 confidence scale (1 = weak/possible, 4 =
# strong/confirmed) -- shown as-is rather than filtered, so low-confidence
# calls stay visible but clearly labeled.
load_crispr_arrays <- function(bacterium_filename) {
  f <- file.path("data", paste0(bacterium_filename, "_crispr_arrays.tsv"))
  if (!file.exists(f)) return(NULL)
  d <- tryCatch(read_tsv(f, show_col_types = FALSE),
                error = function(e) { message("Error reading CRISPR arrays: ", e$message); NULL })
  d
}

# Cas gene clusters (CRISPRCasFinder), same source pipeline as above. A
# bacterium can have CRISPR arrays with no adjacent Cas cluster (orphan
# arrays -- the interference machinery is elsewhere or absent) or vice versa,
# so these are tracked and displayed independently, not assumed paired.
load_cas_systems <- function(bacterium_filename) {
  f <- file.path("data", paste0(bacterium_filename, "_cas_systems.tsv"))
  if (!file.exists(f)) return(NULL)
  d <- tryCatch(read_tsv(f, show_col_types = FALSE),
                error = function(e) { message("Error reading Cas systems: ", e$message); NULL })
  d
}

crispr_file_exists <- function(bacterium_filename) {
  file.exists(file.path("data", paste0(bacterium_filename, "_crispr_arrays.tsv")))
}

# eggNOG-mapper functional annotation -- richer/complementary to the plain
# Prokka+COG annotation used elsewhere (KEGG-style orthologous groups, Pfam
# domains, a preferred gene name), currently only run for the reference
# strain (acutalibacter_muris_kb18). Parsed straight from the GFF attribute
# column, same str_extract approach as load_prokka_gff().
load_eggnog_annotation <- function(bacterium_filename) {
  f <- file.path("data", paste0(bacterium_filename, "_eggnog.gff"))
  if (!file.exists(f)) return(NULL)
  tryCatch({
    read_tsv(f, comment = "#", col_names = FALSE, show_col_types = FALSE) %>%
      select(seqid = 1, attributes = 9) %>%
      transmute(
        Locus_Tag = str_extract(attributes, "ID=([^;]+)") %>% str_remove("ID="),
        COG_category = str_extract(attributes, "em_COG_cat=([^;]*)") %>% str_remove("em_COG_cat="),
        Description = str_extract(attributes, "em_desc=([^;]*)") %>% str_remove("em_desc="),
        OGs = str_extract(attributes, "em_OGs=([^;]*)") %>% str_remove("em_OGs="),
        PFAMs = str_extract(attributes, "em_PFAMs=([^;]*)") %>% str_remove("em_PFAMs="),
        Preferred_Name = str_extract(attributes, "em_Preferred_name=([^;]*)") %>% str_remove("em_Preferred_name="),
        Score = as.numeric(str_extract(attributes, "em_score=([^;]*)") %>% str_remove("em_score=")),
        Evalue = str_extract(attributes, "em_evalue=([^;]*)") %>% str_remove("em_evalue=")
      ) %>%
      filter(!is.na(Locus_Tag), Locus_Tag != "")
  }, error = function(e) { message("Error reading eggNOG annotation: ", e$message); NULL })
}

# NEW: Function to load similarity data
load_similarity_data <- function(bacterium_filename) {
  similarity_file <- file.path("data", paste0(bacterium_filename, "_similarity.csv"))
  if (!file.exists(similarity_file)) { 
    message("Similarity file not found: ", similarity_file)
    return(NULL) 
  }
  similarity_data <- tryCatch({
    read_csv(similarity_file, show_col_types = FALSE)
  }, error = function(e) { 
    message("Error reading similarity CSV file: ", e$message)
    return(NULL) 
  })
  return(similarity_data)
}

load_gem_model <- function(bacterium_filename) {
  rds <- file.path("data", paste0(bacterium_filename, "_model.rds"))
  if (file.exists(rds)) return(readRDS(rds))                    # fast path
  xml <- file.path("data", paste0(bacterium_filename, "_model.xml"))
  if (file.exists(xml)) {
    # parse_gem_model() lives in preparse_models.R, a separate offline script
    # -- it is NOT sourced into the running app (same "run offline, drop in
    # the output" pattern used for AMRFinderPlus/DefenseFinder/etc.), so we
    # can't parse live here. Tell the maintainer to run it once instead of
    # crashing on an undefined function.
    message("GEM model .xml found for ", bacterium_filename,
            " but no .rds cache -- run preparse_models.R once to generate it.")
    return(NULL)
  }
  message("GEM model not found: ", bacterium_filename); NULL
}

# Which tool/version produced each analysis, for citation/reproducibility.
# Recorded manually in data/provenance.tsv (tool names are known from the
# pipelines used; leave version columns blank until confirmed).
load_provenance <- function(bacterium_filename) {
  f <- file.path("data", "provenance.tsv")
  if (!file.exists(f)) return(NULL)
  df <- tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df) || !("bacterium_name_clean" %in% colnames(df))) return(NULL)
  row <- df %>% filter(bacterium_name_clean == !!bacterium_filename)
  if (nrow(row) == 0) return(NULL)
  row[1, , drop = FALSE]
}

# --- HTML CONTENT GENERATION LOGIC ---

get_general_info_html <- function(bacterium_name_clean) {
  data_file_path <- "data/genome_info.tsv"
  if (!file.exists(data_file_path)) { return("No genome data file found.") }
  
  genome_info_df <- tryCatch({
    read_tsv(data_file_path, show_col_types = FALSE)
  }, error = function(e) {
    message("Error reading genome_info.tsv: ", e$message)
    return(NULL)
  })
  if (is.null(genome_info_df)) return("No genome data file found.")
  
  if (!("bacterium_name_clean" %in% colnames(genome_info_df))) {
    possible <- grep("name|filename|id", colnames(genome_info_df), value = TRUE, ignore.case = TRUE)
    if (length(possible) == 0) {
      return("Genome info file present but expected columns not found.")
    } else {
      genome_info_df <- genome_info_df %>% rename(bacterium_name_clean = !!sym(possible[1]))
    }
  }
  
  entry <- genome_info_df %>% filter(.data$bacterium_name_clean == bacterium_name_clean)
  if (nrow(entry) == 0) {
    return(paste("No general information available for:", bacterium_name_clean))
  }
  row <- entry[1, , drop = FALSE]
  
  safe_col <- function(df, colname) {
    if (colname %in% colnames(df)) return(as.character(df[[colname]][1]))
    return(NA_character_)
  }
  
  ncbi_taxonomy_link <- safe_col(row, "ncbi_taxonomy_link")
  ncbi_taxonomy_id   <- safe_col(row, "ncbi_taxonomy_id")
  size_bp            <- safe_col(row, "size_bp")
  gc_content         <- safe_col(row, "gc_content")
  num_genes          <- safe_col(row, "num_genes")
  ncbi_link          <- safe_col(row, "ncbi_link")
  accession          <- safe_col(row, "accession")
  assembly_level     <- safe_col(row, "assembly_level")
  ncbi_refseq        <- safe_col(row, "ncbi_refseq")
  genome_link        <- safe_col(row, "genome_link")
  rRNA_link          <- safe_col(row, "rRNA_link")
  rRNA_sequences     <- safe_col(row, "rRNA_sequences")
  lifestyle          <- safe_col(row, "lifestyle")
  shape              <- safe_col(row, "shape")
  synonyms           <- safe_col(row, "synonyms")
  type_field         <- safe_col(row, "type")
  taxonomy           <- safe_col(row, "taxonomic_classification")
  
  tryCatch({
    html_content <- tags$table(
      class = "table table-striped table-condensed",
      tags$tbody(
        tags$tr(tags$td(tags$strong("Taxonomy:")), tags$td(if (!is.na(taxonomy)) taxonomy else "Not recorded")),
        tags$tr(tags$td(tags$strong("NCBI Taxonomy ID:")), tags$td(if (!is.na(ncbi_taxonomy_link) && nzchar(ncbi_taxonomy_link)) tags$a(href = ncbi_taxonomy_link, target = "_blank", ncbi_taxonomy_id) else ncbi_taxonomy_id)),
        tags$tr(tags$td(tags$strong("Genome Size:")), tags$td(if (!is.na(size_bp)) paste0(size_bp, " bp") else "Not recorded")),
        tags$tr(tags$td(tags$strong("GC Content:")), tags$td(if (!is.na(gc_content)) paste0(gc_content, " %") else "Not recorded")),
        tags$tr(tags$td(tags$strong("Number of genes (NCBI annotation):")), tags$td(if (!is.na(num_genes)) num_genes else "Not recorded")),
        tags$tr(tags$td(tags$strong("Accession number:")), tags$td(if (!is.na(ncbi_link) && nzchar(ncbi_link)) tags$a(href = ncbi_link, target = "_blank", accession) else accession)),
        tags$tr(tags$td(tags$strong("Assembly Level:")), tags$td(if (!is.na(assembly_level)) assembly_level else "Not recorded")),
        tags$tr(tags$td(tags$strong("NCBI assembly accession:")), tags$td(if (!is.na(genome_link) && nzchar(genome_link)) tags$a(href = genome_link, target = "_blank", ncbi_refseq) else ncbi_refseq)),
        tags$tr(tags$td(tags$strong("16S sequence:")), tags$td(if (!is.na(rRNA_link) && nzchar(rRNA_link)) tags$a(href = rRNA_link, target = "_blank", rRNA_sequences) else rRNA_sequences)),
        tags$tr(tags$td(tags$strong("Temperature class:")), tags$td(if (!is.na(lifestyle)) lifestyle else "Not recorded")),
        tags$tr(tags$td(tags$strong("Shape:")), tags$td(if (!is.na(shape)) shape else "Not recorded")),
        tags$tr(tags$td(tags$strong("Synonyms:")), tags$td(if (!is.na(synonyms)) synonyms else "Not recorded")),
        tags$tr(tags$td(tags$strong("Oxygen requirement:")), tags$td(if (!is.na(type_field)) type_field else "Not recorded")),
        tags$br(),
        downloadButton("download_genome_data", "Download Genome", class = "btn-sm")),
      
    )
    return(html_content)
  }, error = function(e) {
    return(tags$p(style="color: red;", paste("Error accessing genome data columns (check genome_info.tsv):", e$message)))
  })
}

# Builds the "Data Provenance" table shown under General Information: which
# tool produced each analysis, plus the actual file-modification date of that
# analysis output (objective, computed directly from the file -- not guessed).
# Source line under the taxonomy diagram, linking to the strain's NCBI Taxonomy page
taxonomy_source_note <- function(bacterium_filename) {
  info <- tryCatch(read_tsv("data/genome_info.tsv", show_col_types = FALSE,
                            col_types = cols(.default = col_character())), error = function(e) NULL)
  row <- if (!is.null(info)) info[info$bacterium_name_clean == bacterium_filename, ] else NULL
  link <- if (!is.null(row) && nrow(row) == 1 && !is.na(row$ncbi_taxonomy_link)) row$ncbi_taxonomy_link else
    "https://www.ncbi.nlm.nih.gov/datasets/taxonomy/"
  tags$p(style = "font-size: 11px; color: #777; margin-top: 6px;",
         "Lineage: ", tags$a(href = link, target = "_blank", "NCBI Taxonomy"), " (accessed 30 Sep 2026)")
}

get_provenance_html <- function(bacterium_filename) {
  prov <- load_provenance(bacterium_filename)

  file_date <- function(suffix) {
    f <- file.path("data", paste0(bacterium_filename, suffix))
    if (file.exists(f)) format(file.info(f)$mtime, "%Y-%m-%d") else NA_character_
  }
  dates <- list(
    genome   = file_date("_prokka.gff"),
    cog      = file_date("_cog.tsv"),
    loc      = file_date("_localization.csv"),
    mobileog = file_date("_mobileog.csv"),
    model    = file_date("_model.xml")
  )

  clean <- function(x) if (is.null(x) || length(x) == 0 || is.na(x) || x == "") NA_character_ else x

  row_html <- function(label, tool, version, date) {
    tool_val <- clean(tool)
    version_val <- clean(version)
    date_val <- clean(date)
    text <- if (is.na(tool_val) && is.na(date_val)) {
      "Not yet available"
    } else {
      t <- if (is.na(tool_val)) "Tool not recorded" else tool_val
      if (!is.na(version_val)) t <- paste0(t, " v", version_val)
      if (!is.na(date_val)) t <- paste0(t, " — generated ", date_val)
      t
    }
    tags$tr(tags$td(tags$strong(label)), tags$td(text))
  }

  if (is.null(prov)) {
    return(tags$p(style = "color:#888; font-size: 13px;",
                   "No provenance record found for this bacterium."))
  }

  notes_val <- clean(prov$notes)

  tags$table(
    class = "table table-striped table-condensed",
    tags$tbody(
      row_html("Genome annotation:", prov$genome_annotation_tool, prov$genome_annotation_version, dates$genome),
      row_html("COG functional annotation:", prov$cog_annotation_tool, prov$cog_annotation_version, dates$cog),
      row_html("Subcellular localization:", prov$localization_tool, prov$localization_version, dates$loc),
      row_html("Mobile genetic elements:", prov$mobileog_tool, prov$mobileog_db_version, dates$mobileog),
      row_html("Metabolic model:", prov$metabolic_model_tool, prov$metabolic_model_version, dates$model),
      if (!is.na(notes_val)) tags$tr(tags$td(tags$strong("Notes:")), tags$td(notes_val)) else NULL
    )
  )
}

# Which tree files exist for a bacterium: whole-genome GBDP tree (always
# "_tree.phy" if present) and, for some bacteria, a separate 16S rRNA gene
# tree ("_16_tree.phy"). Used to decide what toggle(s) to offer in the UI.
phylogeny_available_trees <- function(bacterium_filename) {
  kinds <- character(0)
  if (file.exists(file.path("data", paste0(bacterium_filename, "_tree.phy")))) {
    kinds <- c(kinds, "genome")
  }
  if (file.exists(file.path("data", paste0(bacterium_filename, "_16_tree.phy")))) {
    kinds <- c(kinds, "16s")
  }
  kinds
}

# Enhanced interactive phylogeny plot function with zoom and tooltips.
# tree_kind: "genome" (whole-genome GBDP tree, "_tree.phy") or "16s" (16S
#   rRNA gene tree, "_16_tree.phy").
# branch_mode: "phylogram" uses the tree's actual branch lengths (genetic
#   distance, as computed by FastME/GBDP — the scientifically meaningful
#   view); "cladogram" spaces tips evenly regardless of distance, which is
#   easier to read when a tree has many closely-spaced tips.
get_phylogeny_plot_enhanced <- function(bacterium_filename, similarity_data = NULL,
                                         tree_kind = "genome", branch_mode = "phylogram") {
  suffix <- if (identical(tree_kind, "16s")) "_16_tree.phy" else "_tree.phy"
  tree_path <- file.path("data", paste0(bacterium_filename, suffix))
  if (!file.exists(tree_path)) return(NULL)
  
  tree <- tryCatch({
    read.tree(tree_path)
  }, error = function(e) {
    message("Error reading tree file: ", e$message)
    return(NULL)
  })
  
  if (is.null(tree)) return(NULL)

  # Find the bacterium name (moved up so we can use it to spot "this" tip)
  bact_index <- which(bacteria_filenames == bacterium_filename)
  bact_name <- if (length(bact_index) > 0) bacteria_names[bact_index[1]] else bacterium_filename

  # Identify which tip (if any) is the actual OMM12 strain rather than a
  # reference/comparison taxon pulled in for placement context. Try matching
  # the strain code (e.g. "KB18") as a whole word in the tip label first; if
  # that's not found anywhere (query genomes are sometimes labeled by their
  # own bare assembly accession instead of strain name), fall back to a
  # unique GCA_/GCF_ accession-style tip. If neither is unambiguous, leave it
  # unhighlighted rather than guessing.
  tip_labels <- tree$tip.label
  strain_code <- tail(strsplit(trimws(bact_name), "\\s+")[[1]], 1)
  strain_match <- grepl(paste0("(^|[^A-Za-z0-9])", strain_code, "([^A-Za-z0-9]|$)"),
                         tip_labels, ignore.case = TRUE)
  accession_match <- grepl("^GC[AF]_[0-9]", tip_labels)
  current_tip <- NA_character_
  if (sum(strain_match) == 1) {
    current_tip <- tip_labels[strain_match]
  } else if (sum(accession_match) == 1) {
    current_tip <- tip_labels[accession_match]
  }

  # Create metadata dataframe for tips
  metadata <- data.frame(
    label = tree$tip.label,
    is_current = tree$tip.label == current_tip,
    stringsAsFactors = FALSE
  )

  # Add similarity data if available
  if (!is.null(similarity_data)) {
    # Try to match by different possible column names
    join_col <- if ("strain_name" %in% colnames(similarity_data)) {
      "strain_name"
    } else if ("label" %in% colnames(similarity_data)) {
      "label"
    } else if ("name" %in% colnames(similarity_data)) {
      "name"
    } else {
      colnames(similarity_data)[1]  # Use first column as fallback
    }
    
    join_by_cols <- setNames("label", join_col)
    metadata <- metadata %>%
      left_join(similarity_data, by = join_by_cols)
  }
  
  # Create tooltip text dynamically based on available columns
  metadata <- metadata %>%
    mutate(
      tooltip_text = {
        base_text <- ifelse(is_current,
                             paste0("<b>&#9733; This strain:</b> ", label),
                             paste0("<b>Strain:</b> ", label))

        if (!is.null(similarity_data)) {
          # Add all numeric/character columns from similarity data
          extra_cols <- setdiff(colnames(.), c("label", "tooltip_text", "is_current"))
          for (col in extra_cols) {
            if (!is.na(.[[col]][1])) {
              col_display <- gsub("_", " ", col)
              col_display <- tools::toTitleCase(col_display)
              
              if (is.numeric(.[[col]])) {
                base_text <- paste0(base_text, "<br><b>", col_display, ":</b> ", 
                                    round(.[[col]], 2))
              } else {
                base_text <- paste0(base_text, "<br><b>", col_display, ":</b> ", .[[col]])
              }
            }
          }
        }
        base_text
      }
    )
  
  # Create the tree plot with ggtree
  bl_mode <- if (identical(branch_mode, "cladogram")) "none" else "branch.length"
  tree_label <- if (identical(tree_kind, "16s")) "16S rRNA Gene Tree" else "Genome (GBDP) Tree"

  p <- ggtree(tree, layout = "rectangular", branch.length = bl_mode) %<+% metadata +
    geom_tiplab(
      aes(label = ifelse(is_current, paste0("★ ", label), label),
          tooltip = tooltip_text, data_id = label,
          fontface = ifelse(is_current, "bold", "plain")),
      size = 3.5,
      hjust = -0.1,
    ) +
    geom_tippoint(
      aes(tooltip = tooltip_text, data_id = label, color = is_current, size = is_current),
      alpha = 0.85
    ) +
    scale_color_manual(values = c(`TRUE` = "#c9a227", `FALSE` = "#2c3e57"), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 6, `FALSE` = 4), guide = "none") +
    # Bootstrap support values, parsed straight from the Newick file's internal
    # node labels (previously discarded — only a meaningless auto node ID was shown)
    geom_text2(
      aes(subset = !isTip & !is.na(label) & label != "", label = label),
      size = 3, color = "#8a8f98", hjust = -0.3, vjust = -0.6
    ) +
    geom_nodepoint(
      aes(data_id = node,
          tooltip = ifelse(!is.na(label) & label != "",
                            paste0("<b>Bootstrap support:</b> ", label, "%"),
                            "<b>Bootstrap support:</b> n/a")),
      size = 3,
      color = "#5b7c99",
      alpha = 0.6
    ) +
    theme(
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "#666666", hjust = 0.5)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(title = paste("Phylogenetic Tree of", bact_name),
         subtitle = paste0(tree_label, " — ",
                            if (identical(branch_mode, "cladogram")) "equal tip spacing (topology only)"
                            else "branch length = genetic distance"))

  return(p)
}

# --- CROSS-GENOME GENE INDEX (built once at app startup, shared by all sessions) ---
# Flattens every bacterium's Prokka GFF (gene/product/COG, same source the
# per-organism Genome Map and COG tabs use) into one table tagged with which
# organism each row came from, so a user can search across the whole OMM12
# community instead of one bacterium at a time.
build_all_genes_index <- function() {
  rows <- lapply(seq_along(bacteria_filenames), function(i) {
    fn <- bacteria_filenames[i]
    gff_data <- tryCatch(load_prokka_gff(fn), error = function(e) NULL)
    if (is.null(gff_data) || nrow(gff_data) == 0) return(NULL)

    gff_data %>%
      filter(!is.na(type), type %in% c("CDS", "tRNA", "rRNA", "tmRNA", "misc_RNA")) %>%
      filter(!is.na(start), !is.na(end)) %>%
      transmute(
        Bacterium    = bacteria_names[i],
        Filename     = fn,
        Locus_Tag    = locus_tag,
        Gene         = ifelse(is.na(gene) | gene == "", locus_tag, gene),
        Product      = ifelse(is.na(product) | product == "", type, product),
        Type         = type,
        Contig       = seqid,
        Start        = start,
        End          = end,
        Strand       = strand,
        COG_Category = ifelse(is.na(COG_category) | COG_category == "", NA_character_, COG_category)
      )
  })
  dplyr::bind_rows(rows)
}

all_genes_index <- tryCatch(build_all_genes_index(), error = function(e) {
  message("Failed to build cross-genome gene index: ", e$message)
  data.frame()
})

# --- PUBLICATION STATISTICS (aggregated, loaded once at startup) ---
# Summary tables made by build_publication_stats/make_publication_summary.py
# from the screened OMM12 publication list. The per-publication list itself
# (Dimensions.ai free-version export + Scopus citations) stays local because
# those terms don't allow redistributing it; only these totals are shipped.
PUB_GROUP_COLORS <- c("Uses OMM12" = "#4B1F7A", "Other publications" = "#F4A460")

read_pubstats <- function(name, col_types = NULL) {
  f <- file.path("data", paste0("omm12_pubstats_", name, ".tsv"))
  if (!file.exists(f)) return(NULL)
  d <- tryCatch(read_tsv(f, show_col_types = FALSE, col_types = col_types),
                error = function(e) { message("Error reading ", f, ": ", e$message); NULL })
  if (!is.null(d) && "group" %in% names(d))
    d$group <- factor(d$group, levels = names(PUB_GROUP_COLORS))
  d
}

pubstats_summary <- local({
  d <- read_pubstats("summary", cols(.default = col_character()))
  if (is.null(d)) list() else as.list(setNames(d$value, d$key))
})
pubstats_years     <- read_pubstats("years")
pubstats_countries <- read_pubstats("countries")
pubstats_journals  <- read_pubstats("journals")
pubstats_available <- length(pubstats_summary) > 0

format_pub_date <- function(x) {
  d <- suppressWarnings(as.Date(x))
  if (is.null(x) || is.na(d)) return(ifelse(is.null(x), "", x))
  format(d, "%d %b %Y")
}

# The "Fields of Research (ANZSRC 2020)" column packs multiple semicolon-
# separated codes per row, mixing broad 2-digit parent categories (e.g. "31
# Biological Sciences") with specific 4-digit subfields (e.g. "3107
# Microbiology"). Keep only the specific ones for a chart that isn't just
# dominated by the two or three parent buckets every row inherits.
publication_top_fields <- function(df, n = 10) {
  empty <- data.frame(Field = character(0), Count = integer(0), stringsAsFactors = FALSE)
  if (is.null(df) || !("Fields of Research (ANZSRC 2020)" %in% colnames(df))) return(empty)
  vals <- df[["Fields of Research (ANZSRC 2020)"]]
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(empty)
  parts <- trimws(unlist(strsplit(vals, ";\\s*")))
  parts <- parts[grepl("^\\d{4}\\s", parts)]
  if (length(parts) == 0) return(empty)
  tab <- sort(table(parts), decreasing = TRUE)
  tab <- head(tab, n)
  data.frame(Field = names(tab), Count = as.integer(tab), stringsAsFactors = FALSE)
}

# --- PAIRWISE BACTERIA COMPARISON ---
# genome_info_data is loaded once here (rather than re-read on every call the
# way get_general_info_html() does) so the Compare page can look up either
# side's row cheaply.
genome_info_data <- tryCatch({
  f <- file.path("data", "genome_info.tsv")
  if (file.exists(f)) read_tsv(f, show_col_types = FALSE) else NULL
}, error = function(e) {
  message("Failed to load genome_info.tsv: ", e$message)
  NULL
})

compute_cog_counts <- function(bacterium_filename) {
  data <- load_cog_tsv(bacterium_filename)
  if (is.null(data) || nrow(data) == 0 || !("COG_category" %in% colnames(data))) {
    return(data.frame(COG_category = character(0), Count = integer(0), stringsAsFactors = FALSE))
  }
  data %>%
    mutate(COG_category = toupper(COG_category)) %>%
    mutate(COG_category = ifelse(is.na(COG_category) | COG_category %in% c("", "S", "UNKNOWN"), "S", COG_category)) %>%
    group_by(COG_category) %>%
    summarise(Count = n(), .groups = "drop")
}

compute_localization_counts <- function(bacterium_filename) {
  data <- load_localization_csv(bacterium_filename)
  if (is.null(data) || nrow(data) == 0 || !("Localization" %in% colnames(data))) {
    return(data.frame(Localization = character(0), Count = integer(0), stringsAsFactors = FALSE))
  }
  data %>%
    filter(!is.na(Localization), Localization != "") %>%
    count(Localization, name = "Count")
}

compute_mobileog_counts <- function(bacterium_filename) {
  data <- load_mobileog_csv(bacterium_filename)
  if (is.null(data) || nrow(data) == 0 || !("Major mobileOG Category" %in% colnames(data))) {
    return(data.frame(Category = character(0), Count = integer(0), stringsAsFactors = FALSE))
  }
  data %>%
    filter(!is.na(`Major mobileOG Category`), `Major mobileOG Category` != "") %>%
    count(`Major mobileOG Category`, name = "Count") %>%
    rename(Category = `Major mobileOG Category`)
}

compute_defense_counts <- function(bacterium_filename) {
  data <- load_defense_systems(bacterium_filename)
  if (is.null(data) || nrow(data) == 0 || !("type" %in% colnames(data))) {
    return(data.frame(Type = character(0), Count = integer(0), stringsAsFactors = FALSE))
  }
  data %>%
    count(type, name = "Count") %>%
    rename(Type = type)
}

# --- DEFENSE SYSTEMS OVERVIEW (community-wide) ---
# DefenseFinder output (data/{bacterium}_defense_systems.tsv, or DefenseFinder's own default
# "{bacterium}_defense_finder_systems.tsv" naming -- see load_defense_systems() above) aggregated
# across all 12 OMM12 members. "Defense" systems (RM, CBASS, Wadjet, CRISPR-Cas, etc.) protect the
# host against phages/MGEs; "Antidefense" systems (Anti-CRISPR, Anti-RM, Anti-RecBCD, ...)
# counteract host defenses and are themselves often carried on prophages or other mobile elements,
# so their presence is itself a signal of past/ongoing phage interaction. Any OMM12 member with no
# DefenseFinder run on file yet is shown as "not yet run" (defense_bacteria_missing below) rather
# than silently counted as zero.
build_all_defense_data <- function() {
  rows <- lapply(bacteria_filenames, function(bf) {
    df <- load_defense_systems(bf)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    idx <- match(bf, bacteria_filenames)
    df$Bacterium <- bacteria_names[idx]
    df$Bacterium_filename <- bf
    df
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) return(out)
  out %>%
    mutate(type = trimws(type), subtype = trimws(subtype), activity = trimws(activity))
}
all_defense_systems_data <- tryCatch(build_all_defense_data(), error = function(e) data.frame())

defense_bacteria_with_data <- if (nrow(all_defense_systems_data) > 0) {
  unique(all_defense_systems_data$Bacterium_filename)
} else character(0)
defense_bacteria_missing <- setNames(bacteria_names, bacteria_filenames)[
  setdiff(bacteria_filenames, defense_bacteria_with_data)
]

# --- CRISPR ARRAYS & CAS SYSTEMS OVERVIEW (community-wide) ---
# CRISPRCasFinder output (data/{bacterium}_crispr_arrays.tsv and
# _cas_systems.tsv, flattened offline from the tool's zip bundle by
# build_crispr_bgc_pipeline/extract_crisprcasfinder.py -- see that folder for
# the extraction script). CRISPR arrays and their associated Cas gene
# clusters are tracked as separate three-state pairs (ran-with-data /
# ran-clean-zero-hits / not-run-yet), same pattern as AMRFinderPlus below,
# since a bacterium can genuinely have arrays with no adjacent Cas cluster
# (orphan array) or a real, clean zero result for either.
build_all_crispr_arrays_data <- function() {
  rows <- lapply(bacteria_filenames, function(bf) {
    d <- load_crispr_arrays(bf)
    if (is.null(d) || nrow(d) == 0) return(NULL)
    idx <- match(bf, bacteria_filenames)
    d$Bacterium <- bacteria_names[idx]
    d$Bacterium_filename <- bf
    d
  })
  dplyr::bind_rows(rows)
}
all_crispr_arrays_data <- tryCatch(build_all_crispr_arrays_data(), error = function(e) data.frame())

build_all_cas_systems_data <- function() {
  rows <- lapply(bacteria_filenames, function(bf) {
    d <- load_cas_systems(bf)
    if (is.null(d) || nrow(d) == 0) return(NULL)
    idx <- match(bf, bacteria_filenames)
    d$Bacterium <- bacteria_names[idx]
    d$Bacterium_filename <- bf
    d
  })
  dplyr::bind_rows(rows)
}
all_cas_systems_data <- tryCatch(build_all_cas_systems_data(), error = function(e) data.frame())

crispr_bacteria_ran <- bacteria_filenames[vapply(bacteria_filenames, crispr_file_exists, logical(1))]
crispr_bacteria_with_arrays <- if (nrow(all_crispr_arrays_data) > 0) {
  unique(all_crispr_arrays_data$Bacterium_filename)
} else character(0)
crispr_bacteria_zero_arrays <- setNames(bacteria_names, bacteria_filenames)[
  setdiff(crispr_bacteria_ran, crispr_bacteria_with_arrays)
]
crispr_bacteria_missing <- setNames(bacteria_names, bacteria_filenames)[
  setdiff(bacteria_filenames, crispr_bacteria_ran)
]

# --- BIOSYNTHETIC GENE CLUSTERS / BGC (community-wide) ---
# antiSMASH output (data/{bacterium}_antismash_regions.tsv, flattened offline
# from the tool's zip bundle by build_crispr_bgc_pipeline/extract_antismash.py).
# One row per predicted BGC region: antiSMASH's rule-based product/category
# call, plus (when it found one) the closest known cluster from the MIBiG
# database via knownclusterblast with its percent similarity. Low similarity
# (well under ~50%) usually means "shares a few genes by chance" rather than
# a real match to that specific known product -- shown as-is, not filtered,
# so treat a low-percent hit as "unclassified/novel", not a confirmed call.
build_all_antismash_data <- function() {
  rows <- lapply(bacteria_filenames, function(bf) {
    f <- file.path("data", paste0(bf, "_antismash_regions.tsv"))
    if (!file.exists(f)) return(NULL)
    d <- tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) {
      message("Error reading antiSMASH regions: ", e$message); NULL
    })
    if (is.null(d) || nrow(d) == 0) return(NULL)
    idx <- match(bf, bacteria_filenames)
    d$Bacterium <- bacteria_names[idx]
    d$Bacterium_filename <- bf
    d
  })
  dplyr::bind_rows(rows)
}
all_antismash_data <- tryCatch(build_all_antismash_data(), error = function(e) data.frame())

antismash_file_exists <- function(bacterium_filename) {
  file.exists(file.path("data", paste0(bacterium_filename, "_antismash_regions.tsv")))
}
antismash_bacteria_ran <- bacteria_filenames[vapply(bacteria_filenames, antismash_file_exists, logical(1))]
antismash_bacteria_with_data <- if (nrow(all_antismash_data) > 0) {
  unique(all_antismash_data$Bacterium_filename)
} else character(0)
antismash_bacteria_zero_regions <- setNames(bacteria_names, bacteria_filenames)[
  setdiff(antismash_bacteria_ran, antismash_bacteria_with_data)
]
antismash_bacteria_missing <- setNames(bacteria_names, bacteria_filenames)[
  setdiff(bacteria_filenames, antismash_bacteria_ran)
]
antismash_categories <- if (nrow(all_antismash_data) > 0) {
  sort(unique(unlist(strsplit(all_antismash_data$Category[!is.na(all_antismash_data$Category)], ";"))))
} else character(0)

# --- AMR & SPECIALTY GENES (text-mining screen against existing annotations) ---
# There's no AMRFinderPlus/CARD-RGI/ResFinder/VFDB run anywhere in data/ for any of the
# 12 bacteria (checked: only a 0-byte leftover AMRFinderPlus stub file exists). Running one
# of those tools, or downloading and safely transcribing a curated reference database (CARD,
# VFDB) into this sandbox, isn't something that can be done reliably here. Instead, this is a
# transparent keyword screen against the gene symbols and product descriptions Prokka already
# assigned (all_genes_index, ultimately sourced from UniProt/Swiss-Prot via Prokka's own
# reference databases) -- flagging genes whose name or product text matches well-known
# antimicrobial-resistance or virulence-factor gene family patterns.
# This is NOT equivalent to a curated CARD/ResFinder/VFDB alignment call: no percent identity,
# no coverage, no E-value, and no result at all for any real resistance gene Prokka happened to
# annotate as "hypothetical protein" (a likely source of false negatives). It's the same class
# of method BacDive itself flags with a "text-mining / automatically generated" icon on its own
# strain pages -- a useful lead, not a validated determination.
specialty_gene_rules <- data.frame(
  category    = c("AMR","AMR","AMR","AMR","AMR","AMR","AMR","AMR","AMR",
                   "Virulence","Virulence","Virulence","Virulence"),
  subcategory = c("Aminoglycoside resistance", "Beta-lactam resistance", "Tetracycline resistance",
                   "Macrolide/lincosamide/streptogramin resistance", "Glycopeptide (vancomycin) resistance",
                   "Chloramphenicol resistance", "Sulfonamide/trimethoprim resistance",
                   "Multidrug efflux", "Colistin/polymyxin resistance",
                   "Toxin", "Adhesion/invasion factor", "Secretion system", "Iron acquisition/siderophore"),
  gene_regex = c("^(aac|aph|ant|aad)\\(?", "^bla", "^tet[a-zA-Z0-9]",
                 "^(erm|mef|lnu|lsa|vga|msr)", "^van[a-zA-Z]",
                 "^(cat|cml)", "^(sul|dfr)",
                 "^(mdt|emr|acr|mex)", "^mcr",
                 NA, NA, NA, NA),
  product_regex = c("aminoglycoside", "beta-?lactamase", "tetracycline (resistance|efflux)",
                     "(macrolide|lincosamide|streptogramin) resistance", "vancomycin resistance",
                     "chloramphenicol (acetyltransferase|resistance)", "(sulfonamide resistance|dihydrofolate reductase)",
                     "multidrug (resistance|efflux|transporter)", "(colistin|polymyxin) resistance",
                     "(hemolysin|haemolysin|enterotoxin|cytotoxin)", "(adhesin|fimbrial|pilus assembly)",
                     "type (iii|iv|vi) secretion", "(siderophore|enterobactin)"),
  stringsAsFactors = FALSE
)

classify_specialty_genes <- function() {
  if (is.null(all_genes_index) || nrow(all_genes_index) == 0) return(data.frame())
  genes <- all_genes_index %>% filter(Type == "CDS")
  results <- list()
  for (i in seq_len(nrow(specialty_gene_rules))) {
    rule <- specialty_gene_rules[i, ]
    gene_match <- if (!is.na(rule$gene_regex)) {
      grepl(rule$gene_regex, genes$Gene, ignore.case = TRUE, perl = TRUE)
    } else rep(FALSE, nrow(genes))
    product_match <- if (!is.na(rule$product_regex)) {
      grepl(rule$product_regex, genes$Product, ignore.case = TRUE, perl = TRUE)
    } else rep(FALSE, nrow(genes))
    hit <- gene_match | product_match
    if (any(hit)) {
      hits <- genes[hit, ]
      hits$Category <- rule$category
      hits$Subcategory <- rule$subcategory
      results[[length(results) + 1]] <- hits
    }
  }
  if (length(results) == 0) return(data.frame())
  dplyr::bind_rows(results) %>%
    distinct(Bacterium, Filename, Locus_Tag, Category, Subcategory, .keep_all = TRUE)
}

specialty_genes_data <- tryCatch(classify_specialty_genes(), error = function(e) {
  message("Failed to classify specialty genes: ", e$message); data.frame()
})

# AMR phenotype = real lab-tested susceptibility (MIC / disk-diffusion), which can't be
# predicted from genome sequence alone -- this only ever shows real data, never a guess.
# Expected file: data/{bacterium}_amr_phenotype.tsv with columns like antibiotic, method,
# result, mic_value. Same empty-state pattern as the Similarity Data / Provenance sections.
load_amr_phenotype <- function(bacterium_filename) {
  f <- file.path("data", paste0(bacterium_filename, "_amr_phenotype.tsv"))
  if (!file.exists(f)) return(NULL)
  tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
}

# Curated AMR/virulence/stress gene calls from a real AMRFinderPlus run (run offline by the
# user -- see data/{bacterium}_amrfinder.tsv, standard AMRFinderPlus tabular output). This is
# the real thing (percent identity, percent coverage, matched reference accession) and takes
# priority over the keyword screen above wherever it's available; the keyword screen stays in
# place as a lower-confidence fallback for bacteria that don't have an AMRFinderPlus run yet.
load_amrfinder <- function(bacterium_filename) {
  # AMRFinderPlus's -o flag writes exactly the filename you give it, no extension
  # added -- accept both "{bacterium}_amrfinder.tsv" and the extension-less
  # "{bacterium}_amrfinder" AMRFinderPlus produces by default (`amrfinder -o {bacterium}_amrfinder`).
  candidates <- file.path("data", paste0(bacterium_filename, c("_amrfinder.tsv", "_amrfinder")))
  f <- candidates[file.exists(candidates)][1]
  if (is.na(f) || is.null(f)) return(NULL)
  d <- tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) {
    message("Error reading AMRFinderPlus output: ", e$message); NULL
  })
  if (is.null(d) || !("Element symbol" %in% colnames(d))) return(NULL)
  d
}

amrfinder_file_exists <- function(bacterium_filename) {
  candidates <- file.path("data", paste0(bacterium_filename, c("_amrfinder.tsv", "_amrfinder")))
  any(file.exists(candidates))
}

build_all_amrfinder_data <- function() {
  rows <- lapply(bacteria_filenames, function(bf) {
    df <- load_amrfinder(bf)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    idx <- match(bf, bacteria_filenames)
    df$Bacterium <- bacteria_names[idx]
    df$Bacterium_filename <- bf
    df
  })
  dplyr::bind_rows(rows)
}
all_amrfinder_data <- tryCatch(build_all_amrfinder_data(), error = function(e) {
  message("Failed to build AMRFinderPlus dataset: ", e$message); data.frame()
})

# A bacterium can be in one of three states: no AMRFinderPlus run on file at all
# (amrfinder_bacteria_missing -- falls back to the keyword screen), a run that
# completed and found zero AMR elements (amrfinder_bacteria_zero_hits -- a real,
# clean result, not "unknown"), or a run with one or more hits (feeds into
# all_amrfinder_data above). Collapsing "not run" and "run, found nothing" into
# one bucket would misrepresent a genuine negative result as missing data.
amrfinder_bacteria_ran <- bacteria_filenames[vapply(bacteria_filenames, amrfinder_file_exists, logical(1))]
amrfinder_bacteria_with_data <- if (nrow(all_amrfinder_data) > 0) {
  unique(all_amrfinder_data$Bacterium_filename)
} else character(0)
amrfinder_bacteria_zero_hits <- setNames(bacteria_names, bacteria_filenames)[
  setdiff(amrfinder_bacteria_ran, amrfinder_bacteria_with_data)
]
amrfinder_bacteria_missing <- setNames(bacteria_names, bacteria_filenames)[
  setdiff(bacteria_filenames, amrfinder_bacteria_ran)
]

# --- SEQUENCE SEARCH (MicrobesOnline-style "paste a sequence, search the proteomes") ---
# Same lightweight method validated for data/ortholog_best_hits.tsv: k=4 amino-acid
# k-mer seed index + containment score (|shared kmers| / min(|query kmers|, |target kmers|)).
# Unlike the offline ortholog table (precomputed for the 12x11 fixed protein set), this has
# to work for ANY pasted query, so it's computed live in R. Per-genome k-mer indices are
# built lazily (only the first time a given bacterium is searched) and cached in-memory for
# the life of the R process, so repeat searches and "all 12 bacteria" searches after the
# first stay fast.

# Simple FASTA parser (id -> sequence), vectorized over lines rather than looping char-by-char.
parse_faa_r <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  header_idx <- grep("^>", lines)
  if (length(header_idx) == 0) return(NULL)
  ids <- sub("^>(\\S+).*$", "\\1", lines[header_idx])
  starts <- header_idx + 1L
  ends <- c(header_idx[-1] - 1L, length(lines))
  seqs <- mapply(function(s, e) {
    if (s > e) return("")
    paste(lines[s:e], collapse = "")
  }, starts, ends, SIMPLIFY = TRUE)
  seqs <- toupper(unname(seqs))
  names(seqs) <- ids
  seqs
}

# Vectorized k-mer extraction for one sequence.
get_kmers_r <- function(seq, k = 4) {
  n <- nchar(seq)
  if (is.na(n) || n < k) return(character(0))
  starts <- seq_len(n - k + 1L)
  unique(substring(seq, starts, starts + k - 1L))
}

# Lazy, memoized per-genome k-mer index cache (environment used as a mutable global cache).
.seq_search_cache <- new.env(parent = emptyenv())

get_genome_kmer_index <- function(bacterium_filename, k = 4) {
  cache_key <- paste0(bacterium_filename, "_k", k)
  cached <- mget(cache_key, envir = .seq_search_cache, ifnotfound = list(NULL))[[1]]
  if (!is.null(cached)) return(cached)

  faa_path <- file.path("data", paste0(bacterium_filename, "_protein.faa"))
  seqs <- parse_faa_r(faa_path)
  if (is.null(seqs) || length(seqs) == 0) {
    index <- list(ids = character(0), seqs = character(0), kmers = list())
    assign(cache_key, index, envir = .seq_search_cache)
    return(index)
  }
  kmer_list <- lapply(seqs, get_kmers_r, k = k)
  index <- list(ids = names(seqs), seqs = unname(seqs), kmers = kmer_list)
  assign(cache_key, index, envir = .seq_search_cache)
  index
}

# Search a pasted query sequence against one or more of the 12 proteomes.
# Returns a ranked data.frame of hits (empty data.frame if none pass thresholds).
search_sequence_r <- function(query_seq, target_filenames, k = 4, min_shared = 3,
                               min_score = 0.10, top_n = 25) {
  query_seq <- toupper(gsub("[^A-Za-z]", "", query_seq))
  query_kmers <- get_kmers_r(query_seq, k = k)
  if (length(query_kmers) == 0) return(data.frame())

  bact_lookup <- setNames(bacteria_names, bacteria_filenames)
  all_hits <- list()

  for (bf in target_filenames) {
    idx <- get_genome_kmer_index(bf, k = k)
    if (length(idx$kmers) == 0) next

    shared_counts <- vapply(idx$kmers, function(tk) {
      if (length(tk) == 0) return(0L)
      length(intersect(query_kmers, tk))
    }, integer(1))
    target_lengths <- vapply(idx$kmers, length, integer(1))
    min_len <- pmin(length(query_kmers), target_lengths)
    scores <- ifelse(min_len > 0, shared_counts / min_len, 0)

    keep <- which(shared_counts >= min_shared & scores >= min_score)
    if (length(keep) == 0) next

    all_hits[[bf]] <- data.frame(
      Bacterium_filename = bf,
      Locus_Tag = idx$ids[keep],
      Score = round(scores[keep], 4),
      SharedKmers = shared_counts[keep],
      TargetLength = nchar(idx$seqs[keep]),
      stringsAsFactors = FALSE
    )
  }

  if (length(all_hits) == 0) return(data.frame())

  result <- dplyr::bind_rows(all_hits) %>%
    arrange(desc(Score), desc(SharedKmers))
  if (nrow(result) > top_n) result <- result[seq_len(top_n), ]

  result$Bacterium <- unname(bact_lookup[result$Bacterium_filename])
  result <- result %>%
    left_join(
      all_genes_index %>% select(Filename, Locus_Tag, Gene, Product),
      by = c("Bacterium_filename" = "Filename", "Locus_Tag" = "Locus_Tag")
    ) %>%
    select(Bacterium, Gene, Product, `Locus Tag` = Locus_Tag, `Hit length (aa)` = TargetLength,
           Score, `Shared k-mers` = SharedKmers)

  result
}

# --- PRIMER DESIGN (Primer3-lite, native R) ---
# Candidate PCR primers scored against approximate rules of thumb, not full
# nearest-neighbor thermodynamics: Tm via the common empirical formula
# Tm = 64.9 + 41*(GC_count - 16.4)/length (reasonable for ~18-25 nt oligos,
# but not a substitute for real thermodynamic Tm prediction), GC% in a target
# window, no long homopolymer runs, and a naive 3'-end self-complementarity
# check as a crude hairpin-risk flag. There's no true secondary-structure
# prediction and no forward/reverse primer-dimer cross-check. Treat results
# as a starting shortlist to verify with a dedicated tool (Primer3, IDT
# OligoAnalyzer) before ordering, exactly like the caveat on Sequence Search.

revcomp_dna <- function(seq) {
  paste(rev(chartr("ACGTacgt", "TGCAtgca", strsplit(seq, "")[[1]])), collapse = "")
}

# Lazily-cached whole-genome FASTA per bacterium (reuses the same cache
# environment as the protein k-mer index, just a different key namespace),
# used to slice out a gene's nucleotide sequence via its GFF coordinates.
get_genome_contig_seqs <- function(bacterium_filename) {
  cache_key <- paste0(bacterium_filename, "_genomefasta")
  cached <- mget(cache_key, envir = .seq_search_cache, ifnotfound = list(NULL))[[1]]
  if (!is.null(cached)) return(cached)
  fasta_path <- file.path("data", paste0(bacterium_filename, "_genome.fasta"))
  seqs <- parse_faa_r(fasta_path)  # generic FASTA parser, works fine for genomic contigs too
  if (is.null(seqs)) seqs <- character(0)
  assign(cache_key, seqs, envir = .seq_search_cache)
  seqs
}

# Nucleotide sequence of one gene (by bacterium + locus tag), pulled from the
# whole-genome FASTA using the same Start/End/Strand/Contig already parsed
# into all_genes_index from the Prokka GFF (1-based, inclusive coordinates).
# Reverse-complemented automatically for "-" strand genes.
get_gene_nucleotide_seq <- function(bacterium_filename, locus_tag) {
  if (is.null(all_genes_index) || nrow(all_genes_index) == 0) return(NULL)
  row <- all_genes_index %>% filter(Filename == bacterium_filename, Locus_Tag == locus_tag)
  if (nrow(row) == 0) return(NULL)
  row <- row[1, ]
  contigs <- get_genome_contig_seqs(bacterium_filename)
  if (length(contigs) == 0 || !(row$Contig %in% names(contigs))) return(NULL)
  contig_seq <- contigs[[row$Contig]]
  start <- as.integer(row$Start); end <- as.integer(row$End)
  if (is.na(start) || is.na(end) || start < 1 || end > nchar(contig_seq) || start > end) return(NULL)
  seq <- toupper(substring(contig_seq, start, end))
  if (identical(row$Strand, "-")) seq <- revcomp_dna(seq)
  seq
}

primer_tm <- function(primer) {
  n <- nchar(primer)
  gc <- lengths(regmatches(primer, gregexpr("[GC]", primer)))
  64.9 + 41 * (gc - 16.4) / n
}

primer_gc_percent <- function(primer) {
  n <- nchar(primer)
  gc <- lengths(regmatches(primer, gregexpr("[GC]", primer)))
  100 * gc / n
}

primer_has_poly_run <- function(primer, run_len = 4) {
  grepl(paste0("([ACGT])\\1{", run_len - 1, ",}"), primer)
}

# Crude hairpin-risk proxy: does the primer's own 3' tail have a
# reverse-complement match elsewhere in the primer? Not a real dG/structure
# prediction, just a cheap red flag.
primer_has_3prime_self_complement <- function(primer, k = 5) {
  n <- nchar(primer)
  if (n < k) return(FALSE)
  tail_seq <- substring(primer, n - k + 1, n)
  tail_rc <- revcomp_dna(tail_seq)
  grepl(tail_rc, substring(primer, 1, n - k), fixed = TRUE)
}

find_candidate_primers <- function(seq, region = c("start", "end"),
                                    len_range = 18:25, gc_range = c(40, 60),
                                    tm_target = 60, scan_window = 150) {
  region <- match.arg(region)
  n <- nchar(seq)
  out <- list()
  for (len in len_range) {
    if (n < len) next
    if (region == "start") {
      starts <- seq_len(min(n - len + 1, scan_window))
    } else {
      last_start <- n - len + 1
      starts <- seq(max(1, last_start - scan_window + 1), last_start)
    }
    for (s in starts) {
      primer <- substring(seq, s, s + len - 1)
      if (grepl("[^ACGT]", primer)) next
      gc <- primer_gc_percent(primer)
      if (gc < gc_range[1] || gc > gc_range[2]) next
      if (primer_has_poly_run(primer)) next
      tm <- primer_tm(primer)
      out[[length(out) + 1]] <- data.frame(
        start = s, end = s + len - 1, length = len, sequence = primer,
        gc = round(gc, 1), tm = round(tm, 1),
        hairpin_risk = primer_has_3prime_self_complement(primer),
        score = abs(tm - tm_target),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(out) == 0) return(data.frame())
  dplyr::bind_rows(out)
}

design_primers <- function(seq, product_size_range = c(100, 1000), n_pairs = 5,
                            len_range = 18:25, gc_range = c(40, 60), tm_target = 60) {
  seq <- toupper(gsub("[^ACGTacgt]", "", seq))
  n <- nchar(seq)
  if (n < 80) {
    return(list(error = "Sequence is too short for primer design (need at least ~80 bp).",
                pairs = data.frame(), seq_length = n))
  }

  fwd <- find_candidate_primers(seq, "start", len_range, gc_range, tm_target)
  # Reverse primers: look near the position that gives a product within the
  # requested size range (not at the very end of the gene -- for genes longer
  # than the max product size that could never give a valid pair). Truncating
  # from the 3' end keeps coordinates unchanged.
  rev_seq <- substring(seq, 1, min(n, product_size_range[2]))
  rev_window <- max(150, min(400, product_size_range[2] - product_size_range[1] + 1))
  rev <- find_candidate_primers(rev_seq, "end", len_range, gc_range, tm_target,
                                scan_window = rev_window)
  if (nrow(fwd) == 0 || nrow(rev) == 0) {
    return(list(error = paste("No primer candidates passed the length/GC% filters.",
                               "Try widening the GC% range or primer length range."),
                pairs = data.frame(), seq_length = n))
  }
  fwd <- fwd %>% arrange(score) %>% head(30)
  rev <- rev %>% arrange(score) %>% head(30)

  pairs <- list()
  for (i in seq_len(nrow(fwd))) {
    f <- fwd[i, ]
    for (j in seq_len(nrow(rev))) {
      r <- rev[j, ]
      if (r$start <= f$end) next
      product_size <- r$end - f$start + 1
      if (product_size < product_size_range[1] || product_size > product_size_range[2]) next
      pair_score <- f$score + r$score + abs(f$tm - r$tm) * 1.5 +
        (if (isTRUE(f$hairpin_risk)) 3 else 0) + (if (isTRUE(r$hairpin_risk)) 3 else 0)
      pairs[[length(pairs) + 1]] <- data.frame(
        `Forward primer (5'->3')` = f$sequence, `Fwd length` = f$length,
        `Fwd Tm (C)` = f$tm, `Fwd GC%` = f$gc,
        `Reverse primer (5'->3')` = revcomp_dna(r$sequence), `Rev length` = r$length,
        `Rev Tm (C)` = r$tm, `Rev GC%` = r$gc,
        `Product size (bp)` = product_size, `Tm diff` = round(abs(f$tm - r$tm), 1),
        `Hairpin risk` = (isTRUE(f$hairpin_risk) || isTRUE(r$hairpin_risk)),
        pair_score = round(pair_score, 2),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    }
  }
  if (length(pairs) == 0) {
    return(list(error = paste("No valid primer pairs found within that product-size range.",
                               "Try widening the product size range."),
                pairs = data.frame(), seq_length = n))
  }
  result <- dplyr::bind_rows(pairs) %>% arrange(pair_score) %>% head(n_pairs) %>%
    select(-pair_score)
  list(error = NULL, pairs = result, seq_length = n)
}

# --- GENE-NEIGHBORHOOD / SYNTENY VIEW ---
# Orthologs are resolved via data/ortholog_best_hits.tsv: real protein-sequence
# comparison (k-mer seed index + containment score, k=4 amino acid k-mers)
# computed once offline across all 12 proteomes, with reciprocal-best-hit
# flagged. This replaces an earlier version that aligned bacteria by exact
# gene-symbol text match, which is wrong -- many true orthologs have different
# annotated symbols (or none at all, "hypothetical protein"), and conversely a
# shared symbol doesn't guarantee orthology (paralogs, Prokka's automatic
# "_1"/"_2" suffixing). Validated against 25,666 reciprocal hits where both
# sides happened to also have an annotated gene symbol: 61% matched exactly,
# and most of the "mismatches" on inspection were just naming variants
# (gdh/gdhA, mreB/mreB_1) rather than genuine errors.
# Neighbors are still aligned by genomic RANK position around the anchor, not
# by bp distance, since intergenic spacing differs wildly between species.

ortholog_hits <- tryCatch({
  f <- file.path("data", "ortholog_best_hits.tsv")
  if (file.exists(f)) read_tsv(f, show_col_types = FALSE) else NULL
}, error = function(e) {
  message("Failed to load ortholog_best_hits.tsv: ", e$message)
  NULL
})

# --- PAN-GENOME / ORTHOLOG MATRIX ---
# Orthogroups = connected components of the reciprocal-best-hit graph built
# from ortholog_hits (Union-Find over all 41,553 genes across the 12
# proteomes, computed offline in Python -- same underlying method as the Gene
# Neighborhood page). This is a lightweight stand-in for a proper OrthoFinder
# run: it will undercount true core genes (a genuinely conserved gene can
# still fail to clear the reciprocal-best-hit bar between two very divergent
# members) and can't distinguish orthologs from very close paralogs as
# cleanly as a phylogeny-aware tool would.
orthogroup_summary_data <- tryCatch({
  f <- file.path("data", "orthogroup_summary.tsv")
  if (file.exists(f)) read_tsv(f, show_col_types = FALSE) else NULL
}, error = function(e) { message("Failed to load orthogroup_summary.tsv: ", e$message); NULL })

orthogroup_genes_data <- tryCatch({
  f <- file.path("data", "orthogroups.tsv")
  if (file.exists(f)) read_tsv(f, show_col_types = FALSE) else NULL
}, error = function(e) { message("Failed to load orthogroups.tsv: ", e$message); NULL })

# Enrich the summary with a representative gene symbol + sample product per
# orthogroup, so the browser table isn't just a wall of "OG_00001" IDs.
orthogroup_summary_full <- if (!is.null(orthogroup_summary_data) && !is.null(orthogroup_genes_data)) {
  repr <- orthogroup_genes_data %>%
    group_by(orthogroup_id) %>%
    summarise(
      representative_gene = {
        g <- gene[!is.na(gene) & gene != ""]
        if (length(g) == 0) NA_character_ else names(sort(table(g), decreasing = TRUE))[1]
      },
      sample_product = {
        p <- product[!is.na(product) & product != ""]
        if (length(p) == 0) NA_character_ else p[1]
      },
      .groups = "drop"
    )
  orthogroup_summary_data %>% left_join(repr, by = "orthogroup_id")
} else NULL

# Per-bacterium genome composition by orthogroup classification (core / soft-
# core / shell / unique), for the stacked composition chart.
pangenome_composition <- if (!is.null(orthogroup_genes_data) && !is.null(orthogroup_summary_data)) {
  orthogroup_genes_data %>%
    left_join(orthogroup_summary_data %>% select(orthogroup_id, classification), by = "orthogroup_id") %>%
    count(bacterium, classification, name = "n")
} else NULL

# --- OMM12 RELATEDNESS TREE (all 12 members, from k-mer ortholog similarity) ---
# A stand-in for a real phylogeny, covering all 12 OMM12 members (the actual
# GBDP/16S tree only exists for the reference strain, KB18). For each ordered
# pair of bacteria, we take the mean k-mer containment score across their
# reciprocal-best-hit orthologs (data/ortholog_best_hits.tsv) as a proteome-wide
# similarity — conceptually close to an AAI (average amino-acid identity)
# estimate, though computed with a cheaper k-mer method rather than real
# alignment. Every one of the 132 ordered pairs has at least 72 reciprocal
# orthologs behind it, so no pair needs a fallback. Symmetrized similarity is
# converted to a distance (1 - similarity) and clustered with UPGMA
# (hclust, method = "average") to produce a dendrogram. This reflects overall
# proteome similarity, not a validated evolutionary model — it will not
# exactly match a real 16S/GBDP tree, and is presented with that caveat.
build_omm12_relatedness_matrix <- function() {
  if (is.null(ortholog_hits) || nrow(ortholog_hits) == 0) return(NULL)
  recip <- ortholog_hits %>% filter(reciprocal == TRUE)
  if (nrow(recip) == 0) return(NULL)

  pair_scores <- recip %>%
    group_by(source_bacterium, target_bacterium) %>%
    summarise(mean_score = mean(score), n_orthologs = dplyr::n(), .groups = "drop")

  n <- length(bacteria_filenames)
  sim <- matrix(NA_real_, n, n, dimnames = list(bacteria_filenames, bacteria_filenames))
  for (i in seq_len(nrow(pair_scores))) {
    a <- pair_scores$source_bacterium[i]; b <- pair_scores$target_bacterium[i]
    if (a %in% bacteria_filenames && b %in% bacteria_filenames) sim[a, b] <- pair_scores$mean_score[i]
  }

  sim_sym <- matrix(1, n, n, dimnames = dimnames(sim))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      vals <- c(sim[i, j], sim[j, i])
      vals <- vals[!is.na(vals)]
      sim_sym[i, j] <- if (length(vals) > 0) mean(vals) else 0
    }
  }
  list(similarity = sim_sym, pair_scores = pair_scores)
}
omm12_relatedness <- tryCatch(build_omm12_relatedness_matrix(), error = function(e) {
  message("Failed to build OMM12 relatedness matrix: ", e$message); NULL
})

get_omm12_relatedness_plot <- function() {
  if (is.null(omm12_relatedness)) return(NULL)
  sim <- omm12_relatedness$similarity
  dist_mat <- as.dist(1 - sim)
  hc <- hclust(dist_mat, method = "average")
  tree <- ape::as.phylo(hc)
  tree$tip.label <- bacteria_names[match(tree$tip.label, bacteria_filenames)]

  metadata <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE) %>%
    mutate(tooltip_text = paste0("<b>", label, "</b>"))

  ggtree(tree, layout = "rectangular") %<+% metadata +
    geom_tiplab(aes(label = label, tooltip = tooltip_text, data_id = label),
                size = 3.6, hjust = -0.05) +
    geom_tippoint(aes(tooltip = tooltip_text, data_id = label),
                  color = "#c9a227", size = 4, alpha = 0.9) +
    theme(
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "#666666", hjust = 0.5)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(title = "OMM12 Relatedness (all 12 members)",
         subtitle = "UPGMA clustering of k-mer ortholog similarity — not a validated phylogeny")
}

# --- OMM12 16S rRNA GENE PHYLOGENY (all 12 members, real tree) ---
# A genuine single-locus molecular phylogeny covering all 12 OMM12 members,
# built entirely offline (the build environment has no internet access to
# install MAFFT/IQ-TREE/RAxML/etc., so none of the usual tools were
# available -- see build_omm12_16s_tree/ at the repo root for the full
# from-scratch pipeline and its unit tests). Method: one representative 16S
# rRNA gene copy per genome (the longest annotated copy, from each
# bacterium's Prokka GFF), aligned with a custom Gotoh affine-gap-penalty
# multiple sequence alignment (center-star method, medoid-anchored),
# pairwise distances corrected with the Jukes-Cantor (1969) substitution
# model, a tree built by Neighbor-Joining (Saitou & Nei 1987), and support
# from 200 bootstrap replicates (column resampling, fixed RNG seed for
# reproducibility). Every component (pairwise aligner, JC correction,
# Neighbor-Joining) was validated against known-correct reference cases
# before being run on real data. This has real alignment + a substitution
# model + bootstrap support -- unlike the k-mer/UPGMA proteome-similarity
# heatmap alongside it, which is a different (whole-proteome AAI-like)
# signal, not a phylogeny. As with any single-gene tree, it has limited
# power to resolve the deepest, most ancient splits between different phyla
# -- nodes with low bootstrap support (well under 70%) should be read with
# that caution in mind. Raw sequences, the alignment, and a full methods
# report are in data/omm12_16s_{raw,msa,tree,pipeline_report}.* for anyone
# who wants to check the work.
load_omm12_16s_tree <- function() {
  f <- file.path("data", "omm12_16s_tree.nwk")
  if (!file.exists(f)) return(NULL)
  tryCatch(ape::read.tree(f), error = function(e) {
    message("Error reading OMM12 16S tree: ", e$message); NULL
  })
}
omm12_16s_tree <- tryCatch(load_omm12_16s_tree(), error = function(e) {
  message("Failed to load OMM12 16S tree: ", e$message); NULL
})

get_omm12_16s_tree_plot <- function() {
  tree <- omm12_16s_tree
  if (is.null(tree)) return(NULL)
  tree$tip.label <- bacteria_names[match(tree$tip.label, bacteria_filenames)]

  metadata <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE) %>%
    mutate(tooltip_text = paste0("<b>", label, "</b>"))

  ggtree(tree, layout = "rectangular") %<+% metadata +
    geom_tiplab(aes(label = label, tooltip = tooltip_text, data_id = label),
                size = 3.6, hjust = -0.05) +
    geom_tippoint(aes(tooltip = tooltip_text, data_id = label),
                  color = "#c9a227", size = 4, alpha = 0.9) +
    # Bootstrap support values, parsed straight from the Newick file's
    # internal node labels (200 replicates -- see the pipeline report).
    geom_text2(
      aes(subset = !isTip & !is.na(label) & label != "", label = label),
      size = 3, color = "#8a8f98", hjust = -0.3, vjust = -0.6
    ) +
    geom_nodepoint(
      aes(data_id = node,
          tooltip = ifelse(!is.na(label) & label != "",
                            paste0("<b>Bootstrap support:</b> ", label, "%"),
                            "<b>Bootstrap support:</b> n/a")),
      size = 3,
      color = "#5b7c99",
      alpha = 0.6
    ) +
    theme(
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "#666666", hjust = 0.5)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(title = "OMM12 16S rRNA Gene Phylogeny (all 12 members)",
         subtitle = "Custom alignment · Jukes-Cantor distance · Neighbor-Joining · 200 bootstrap replicates")
}

# --- GENE/PROTEIN TREE (per orthogroup) ---
# BV-BRC-style "build a tree from a specific gene family" view, scoped to one
# orthogroup at a time (picked from the Pan-genome browser) rather than a
# whole-genome tree. Pairwise similarity between two members of the orthogroup
# is computed live from the same per-genome k-mer indices built for Sequence
# Search (get_genome_kmer_index(), defined below) -- this guarantees a
# complete pairwise matrix even for member pairs that never showed up as each
# other's reciprocal best hit in the offline ortholog table (orthogroups are
# connected components of the RBH graph, so two members can be linked only
# transitively through a third bacterium).
compute_pairwise_kmer_score <- function(bf1, lt1, bf2, lt2, k = 4) {
  idx1 <- get_genome_kmer_index(bf1, k = k)
  idx2 <- get_genome_kmer_index(bf2, k = k)
  k1 <- idx1$kmers[[lt1]]
  k2 <- idx2$kmers[[lt2]]
  if (is.null(k1) || is.null(k2) || length(k1) == 0 || length(k2) == 0) return(NA_real_)
  length(intersect(k1, k2)) / min(length(k1), length(k2))
}

build_orthogroup_tree <- function(orthogroup_id_val, max_members = 14) {
  if (is.null(orthogroup_genes_data)) {
    return(list(error = "Orthogroup data not available.", plot = NULL, n_members = 0))
  }
  members <- orthogroup_genes_data %>% filter(orthogroup_id == orthogroup_id_val)
  if (nrow(members) < 3) {
    return(list(error = paste0("This orthogroup has ", nrow(members),
                                " member(s) — need at least 3 to draw a tree."),
                plot = NULL, n_members = nrow(members)))
  }
  truncated <- nrow(members) > max_members
  if (truncated) members <- members[seq_len(max_members), ]

  n <- nrow(members)
  bact_label <- bacteria_names[match(members$bacterium, bacteria_filenames)]
  gene_label <- ifelse(is.na(members$gene) | members$gene == "", members$locus_tag, members$gene)
  labels <- make.unique(paste0(bact_label, " (", gene_label, ")"))

  sim <- matrix(1, n, n, dimnames = list(labels, labels))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i >= j) next
      s <- compute_pairwise_kmer_score(members$bacterium[i], members$locus_tag[i],
                                        members$bacterium[j], members$locus_tag[j])
      if (is.na(s)) s <- 0.05
      sim[i, j] <- s; sim[j, i] <- s
    }
  }

  hc <- hclust(as.dist(1 - sim), method = "average")
  tree <- ape::as.phylo(hc)

  metadata <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE) %>%
    mutate(tooltip_text = paste0("<b>", label, "</b>"))

  p <- ggtree(tree, layout = "rectangular") %<+% metadata +
    geom_tiplab(aes(label = label, tooltip = tooltip_text, data_id = label),
                size = 3.3, hjust = -0.05) +
    geom_tippoint(aes(tooltip = tooltip_text, data_id = label),
                  color = "#c9a227", size = 3.5, alpha = 0.9) +
    theme(
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "none",
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.5))) +
    labs(title = paste("Gene tree:", orthogroup_id_val))

  list(error = if (truncated) paste0("Showing the first ", max_members, " of ", nrow(orthogroup_genes_data %>% filter(orthogroup_id == orthogroup_id_val)), " members.") else NULL,
       plot = p, n_members = n)
}

# Real gene symbols per bacterium (excludes cases where Gene just fell back to
# the Prokka locus tag), used to populate the per-bacterium anchor gene picker.
get_gene_choices_for <- function(bacterium_filename) {
  if (is.null(all_genes_index) || nrow(all_genes_index) == 0) return(character(0))
  genes <- all_genes_index %>% filter(Filename == bacterium_filename, Gene != Locus_Tag)
  sort(unique(genes$Gene))
}

# Look up the best ortholog (highest score, preferring a reciprocal best hit)
# of one locus tag in a target bacterium.
lookup_ortholog <- function(source_filename, source_locus_tag, target_filename) {
  if (is.null(ortholog_hits) || nrow(ortholog_hits) == 0) return(NULL)
  hits <- ortholog_hits %>%
    filter(source_bacterium == source_filename, source_locus_tag == !!source_locus_tag,
           target_bacterium == target_filename) %>%
    arrange(desc(reciprocal), desc(score))
  if (nrow(hits) == 0) return(NULL)
  hits[1, ]
}

build_neighborhood_data <- function(anchor_filename, anchor_locus_tag, window = 3) {
  empty <- data.frame()
  if (is.null(anchor_filename) || is.null(anchor_locus_tag) || !nzchar(anchor_locus_tag) ||
      is.null(all_genes_index) || nrow(all_genes_index) == 0) {
    return(empty)
  }
  window <- max(1, min(8, as.integer(window)))

  neighborhood_around <- function(fn, locus_tag) {
    genes <- all_genes_index %>% filter(Filename == fn) %>% arrange(Contig, Start)
    if (nrow(genes) == 0) return(NULL)
    anchor_pos_global <- which(genes$Locus_Tag == locus_tag)
    if (length(anchor_pos_global) == 0) return(NULL)
    anchor_contig <- genes$Contig[anchor_pos_global[1]]
    contig_genes <- genes %>% filter(Contig == anchor_contig)
    anchor_pos <- which(contig_genes$Locus_Tag == locus_tag)[1]
    lo <- max(1, anchor_pos - window)
    hi <- min(nrow(contig_genes), anchor_pos + window)
    neigh <- contig_genes[lo:hi, ]
    neigh$RelPos <- (lo:hi) - anchor_pos
    neigh$IsAnchor <- neigh$RelPos == 0
    neigh
  }

  pieces <- lapply(seq_along(bacteria_filenames), function(i) {
    fn <- bacteria_filenames[i]
    if (fn == anchor_filename) {
      locus_tag <- anchor_locus_tag
      is_ortholog <- TRUE
      recip <- TRUE
      oscore <- NA_real_
    } else {
      hit <- lookup_ortholog(anchor_filename, anchor_locus_tag, fn)
      if (is.null(hit)) return(NULL)
      locus_tag <- hit$target_locus_tag[1]
      recip <- isTRUE(hit$reciprocal[1])
      oscore <- hit$score[1]
      is_ortholog <- TRUE
    }
    neigh <- neighborhood_around(fn, locus_tag)
    if (is.null(neigh)) return(NULL)
    neigh$BactIndex <- i
    neigh$OrthologReciprocal <- recip
    neigh$OrthologScore <- oscore
    neigh
  })
  dplyr::bind_rows(pieces)
}

# Shared prep for both the gggenes plot and its companion detail table:
# resolves the neighborhood, tags each gene's role (anchor / reciprocal
# ortholog / one-directional ortholog / plain neighbor), and orders the
# bacterium facets with the anchor's own genome first.
build_neighborhood_display_data <- function(anchor_filename, anchor_locus_tag, window = 3) {
  neigh <- build_neighborhood_data(anchor_filename, anchor_locus_tag, window)
  if (nrow(neigh) == 0) return(NULL)

  anchor_bact_idx <- match(anchor_filename, bacteria_filenames)

  neigh <- neigh %>%
    mutate(
      xmin = RelPos - 0.45,
      xmax = RelPos + 0.45,
      forward = Strand == "+",
      is_anchor_bact = BactIndex == anchor_bact_idx,
      role = case_when(
        RelPos == 0 & is_anchor_bact  ~ "Anchor gene",
        RelPos == 0 & OrthologReciprocal ~ "Ortholog (reciprocal best hit)",
        RelPos == 0 & !OrthologReciprocal ~ "Ortholog (one-directional hit)",
        TRUE ~ "Neighboring gene"
      ),
      role = factor(role, levels = c("Anchor gene", "Ortholog (reciprocal best hit)",
                                      "Ortholog (one-directional hit)", "Neighboring gene")),
      label = ifelse(RelPos == 0, paste0("* ", Gene), Gene)
    )

  # Facet/row order: anchor bacterium first, then the rest in their usual order.
  bact_levels <- c(
    neigh$Bacterium[neigh$is_anchor_bact][1],
    setdiff(bacteria_names[sort(unique(neigh$BactIndex))], neigh$Bacterium[neigh$is_anchor_bact][1])
  )
  neigh$Bacterium <- factor(neigh$Bacterium, levels = bact_levels)
  neigh
}

# Publication-style gene-arrow diagram via gggenes, one row (facet) per
# bacterium. gggenes' geom_gene_arrow()/geom_gene_label() draw real pointed
# gene glyphs with in-arrow labels instead of the hand-rolled polygons this
# used before -- the standard R package for exactly this kind of figure.
# Note: gggenes' arrow geom is a custom grob that plotly::ggplotly() can't
# reliably convert, so this is rendered as a static ggplot (renderPlot), with
# a companion table below carrying the per-gene detail that used to live in
# hover tooltips.
render_gene_neighborhood_plot <- function(anchor_filename, anchor_locus_tag, window = 3) {
  neigh <- build_neighborhood_display_data(anchor_filename, anchor_locus_tag, window)
  if (is.null(neigh)) return(NULL)

  anchor_gene_label <- neigh$Gene[neigh$RelPos == 0 & neigh$is_anchor_bact][1]
  n_bact <- length(unique(neigh$BactIndex))

  # gggenes' API is being used without an R session available to verify it
  # against the installed package version, and gggenes may not even be
  # installed -- fall back to a plain-text error plot rather than crashing
  # this output if geom_gene_arrow()/geom_gene_label() don't behave as
  # expected here.
  p <- tryCatch({
    ggplot(neigh, aes(xmin = xmin, xmax = xmax, y = 1, fill = role, forward = forward)) +
      geom_gene_arrow(arrowhead_height = unit(4, "mm"), arrowhead_width = unit(2.5, "mm"),
                       arrow_body_height = unit(3, "mm"), color = "white", linewidth = 0.25) +
      geom_gene_label(aes(label = label), align = "centre", height = unit(3, "mm"),
                       padding.x = unit(1, "mm"), grow = FALSE, fontface = "plain") +
      facet_wrap(~ Bacterium, ncol = 1, scales = "fixed", strip.position = "left") +
      scale_fill_manual(
        values = c("Anchor gene" = "#c9a227", "Ortholog (reciprocal best hit)" = "#1b263b",
                   "Ortholog (one-directional hit)" = "#5b7c99", "Neighboring gene" = "#c3ccd9"),
        name = NULL
      ) +
      scale_x_continuous(breaks = seq(-window, window, by = 1),
                          labels = seq(-window, window, by = 1)) +
      labs(x = "Genes from anchor (rank order, not to scale)", y = NULL,
           title = paste("Gene neighborhood:", anchor_gene_label),
           caption = "* marks the gene aligned to the anchor in each row (the anchor itself, or its ortholog)") +
      theme_genes() +
      theme(
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 10),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing.y = unit(6, "pt"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(color = "#666666", size = 9, hjust = 0)
      )
  }, error = function(e) {
    message("gggenes rendering failed: ", e$message)
    ggplot() +
      annotate("text", x = 0, y = 0, size = 4, color = "#a33", hjust = 0.5,
               label = paste0("Could not render with gggenes (", e$message, ").\n",
                               "Check that the gggenes package is installed: install.packages(\"gggenes\").")) +
      theme_void()
  })

  list(plot = p, n_bacteria = n_bact)
}

# Gene-arrow diagram for one antiSMASH-predicted BGC region, mirroring
# antiSMASH's own region viewer: real genomic coordinates on the x-axis (not
# rank order, unlike the Gene Neighborhood plot -- everything here is
# already on one contiguous stretch of one bacterium's genome), arrows
# colored by antiSMASH's own gene_kind classification (core biosynthetic /
# additional biosynthetic / transport / regulatory / other), from
# data/{bacterium}_antismash_genes.tsv (build_crispr_bgc_pipeline/extract_antismash.py).
render_bgc_region_plot <- function(bacterium_filename, region_id) {
  f <- file.path("data", paste0(bacterium_filename, "_antismash_genes.tsv"))
  if (!file.exists(f)) return(NULL)
  genes <- tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(genes) || nrow(genes) == 0) return(NULL)
  genes <- genes %>% filter(Region == region_id)
  if (nrow(genes) == 0) return(NULL)

  genes <- genes %>%
    mutate(
      forward = Strand == "+",
      role = ifelse(is.na(Gene_Kind) | Gene_Kind == "", "other", Gene_Kind),
      label = ifelse(!is.na(Gene_Name) & Gene_Name != "" & Gene_Name != paste0(Locus_Tag, "_gene"),
                      Gene_Name, Locus_Tag)
    )

  role_levels <- c("biosynthetic", "biosynthetic-additional", "transport", "regulatory", "resistance", "other")
  role_labels <- c(biosynthetic = "Core biosynthetic", `biosynthetic-additional` = "Additional biosynthetic",
                    transport = "Transport-related", regulatory = "Regulatory",
                    resistance = "Resistance", other = "Other")
  role_colors <- c(biosynthetic = "#8b1a1a", `biosynthetic-additional` = "#e8817e",
                    transport = "#5b9bd5", regulatory = "#2e8b57",
                    resistance = "#bdbdbd", other = "#999999")
  genes$role <- factor(genes$role, levels = role_levels, labels = role_labels)
  names(role_colors) <- role_labels[names(role_colors)]

  p <- tryCatch({
    ggplot(genes, aes(xmin = Start, xmax = End, y = 1, fill = role, forward = forward)) +
      geom_gene_arrow(arrowhead_height = unit(4, "mm"), arrowhead_width = unit(2.5, "mm"),
                       arrow_body_height = unit(3, "mm"), color = "white", linewidth = 0.25) +
      geom_gene_label(aes(label = label), align = "centre", height = unit(3, "mm"),
                       padding.x = unit(1, "mm"), grow = FALSE, fontface = "plain") +
      scale_fill_manual(values = role_colors, name = NULL, drop = TRUE) +
      labs(x = "Position (bp)", y = NULL, title = paste("BGC", region_id)) +
      theme_genes() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
      )
  }, error = function(e) {
    message("BGC gggenes rendering failed: ", e$message)
    NULL
  })
  p
}

# --- CUSTOM JS FOR SCROLLING ---
js_scroll_code <- tags$script(HTML('
  // Custom message handler for smooth scrolling
  Shiny.addCustomMessageHandler("scrollTo", function(message) {
    var targetId = message.id;
    var targetOffset = $("#" + targetId).offset() ? $("#" + targetId).offset().top : 0;
    var adjustedOffset = targetOffset - 50;
    $("html, body").animate({ scrollTop: adjustedOffset }, 500);
  });

  // Interactive circular genome viewer (CGView.js), replacing the static PNG.
  // Two gotchas this works around:
  //  - The "Circular" pane is a hidden (display:none) Bootstrap tab until
  //    clicked, and CGView/D3 measure zero width/height and draw nothing if
  //    initialized while hidden -- even after the tab is later shown. So we
  //    defer drawing until the tab is actually visible.
  //  - #cgview_container is inside a Shiny uiOutput, so it may not exist in
  //    the DOM yet the instant the message arrives -- so we retry briefly.
  var cgvViewer = null;
  var cgvDrawnUrl = null;
  var cgvPendingUrl = null;

  function cgvDraw(url) {
    var container = document.getElementById("cgview_container");
    if (!container || container.offsetParent === null) return; // not in DOM, or tab hidden
    if (cgvDrawnUrl === url) return; // already showing this one

    container.innerHTML = "<p style=\'color:#666; padding: 20px;\'>Loading genome map…</p>";
    fetch(url)
      .then(function(response) {
        if (!response.ok) { throw new Error("File not found: " + url); }
        return response.text();
      })
      .then(function(gbText) {
        container.innerHTML = "";
        var seqFile = new CGParse.SequenceFile(gbText);
        var cgvJSON = seqFile.toCGViewJSON();
        cgvViewer = new CGView.Viewer("#cgview_container", { width: container.clientWidth || 800, height: 600 });
        cgvViewer.io.loadJSON(cgvJSON);
        cgvViewer.draw();
        cgvDrawnUrl = url;
      })
      .catch(function(err) {
        console.error("CGView load error:", err);
        container.innerHTML = "<p style=\'color:#a33;\'>Could not load the circular genome map for this bacterium (GenBank file not available yet).</p>";
        cgvDrawnUrl = null;
      });
  }

  function cgvWaitForContainer(attemptsLeft) {
    var container = document.getElementById("cgview_container");
    if (container) { cgvDraw(cgvPendingUrl); return; }
    if (attemptsLeft > 0) { setTimeout(function() { cgvWaitForContainer(attemptsLeft - 1); }, 250); }
  }

  Shiny.addCustomMessageHandler("loadCGView", function(message) {
    cgvPendingUrl = message.gbkUrl;
    cgvDrawnUrl = null; // force a redraw even if the tab is already visible
    cgvWaitForContainer(20); // retry for up to ~5s while the UI renders
  });

  // Draw (or catch up) whenever the Circular tab itself is opened
  $(document).on("shown.bs.tab", \'a[data-toggle="tab"]\', function(e) {
    if ($(e.target).text().trim() === "Circular" && cgvPendingUrl) {
      cgvDraw(cgvPendingUrl);
    }
  });
'))

# --- Tools Guide content (static reference data for the "Tools Guide" page
# and the short section blurbs that link back to it) ---
tools_guide_categories <- list(
  list(title = "Genome & Functional Annotation", status = "primary", color = "#1b263b", entries = list(
    list(icon = "map-location-dot", name = "Genome Map (Prokka)",
         identifies = "Every gene on the genome (CDS, tRNA, rRNA) with its coordinates -- shown to you as the Genome Viewer (linear + circular maps) and the Annotation feature table.",
         use = "The base gene catalog everything else in this resource is built on.",
         version = "Prokka 1.13", link = "https://doi.org/10.1093/bioinformatics/btu153",
         link_label = "Seemann 2014, Bioinformatics"),
    list(icon = "layer-group", name = "COG Categories",
         identifies = "Which broad functional category (metabolism, cell processes, etc.) each gene belongs to -- shown as a bar chart of gene counts per category.",
         use = "A quick read on what a genome is functionally built to do.",
         version = "COGclassifier (version not recorded in output)", link = "https://github.com/moshi4/COGclassifier",
         link_label = "COGclassifier on GitHub"),
    list(icon = "circle-nodes", name = "eggNOG Annotation",
         identifies = "Orthologous groups, gene names, and Pfam domains via deeper orthology-based annotation -- shown as a supplemental table alongside the COG chart.",
         use = "Adds detail Prokka's basic gene calls don't capture.",
         version = "eggNOG-mapper 2.1.12", link = "https://doi.org/10.1093/molbev/msab293",
         link_label = "Cantalapiedra et al. 2021, Mol Biol Evol"),
    list(icon = "diagram-project", name = "KEGG Pathways",
         identifies = "Which metabolic/signaling pathways (glycolysis, amino acid synthesis, etc.) a gene belongs to -- shown as a collapsible Main pathway -> Sub-pathway -> gene table accordion.",
         use = "Shows which biochemical pathways a strain can actually run.",
         version = "KEGG Orthology (KO) reference database", link = "https://www.genome.jp/kegg/",
         link_label = "Kanehisa et al. 2023, Nucleic Acids Res"),
    list(icon = "share-nodes", name = "MobileOG (Mobile Genetic Elements)",
         identifies = "Transposons, integrons, and other mobile genetic elements -- shown as a bar chart of hits per category.",
         use = "Flags genes that can jump between genomes or be horizontally transferred.",
         version = "mobileOG-db, DIAMOND homology + keyword search (DB version not recorded)",
         link = "https://doi.org/10.1128/aem.00991-22", link_label = "Brown et al. 2022, Appl Environ Microbiol"),
    list(icon = "location-dot", name = "Subcellular Localization",
         identifies = "Predicted location of each protein (cytoplasm, membrane, secreted) -- shown as a bar chart, a PCA scatter of each protein's probability vector, a filterable gene table, and a localization x COG cross-tab.",
         use = "Helps spot surface-exposed or secreted proteins of interest.",
         version = "DeepLocPro (version not recorded in output)", link = "https://doi.org/10.1093/bioinformatics/btae677",
         link_label = "Moreno et al. 2024, Bioinformatics"),
    list(icon = "atom", name = "Genome-Scale Metabolic Model",
         identifies = "The full set of metabolic reactions a strain's genome can carry out -- shown as a pathway/subsystem bar chart, a reactions table with MetaCyc links, and a force-directed metabolite<->reaction network. KB18 only so far.",
         use = "Predicts growth/nutrient requirements and potential metabolic cross-feeding with other members.",
         version = "gapseq (version not recorded in output)", link = "https://doi.org/10.1186/s13059-021-02295-1",
         link_label = "Zimmermann et al. 2021, Genome Biology")
  )),
  list(title = "Defense & Resistance", status = "danger", color = "#7a2222", entries = list(
    list(icon = "shield-halved", name = "Defense Systems (DefenseFinder)",
         identifies = "Anti-phage immune systems (restriction-modification, abortive infection, etc.) -- shown as bar charts (types per bacterium, defense vs. antidefense) plus a browsable table.",
         use = "How well a strain resists phage infection and foreign DNA invasion.",
         version = "DefenseFinder (version not recorded in output)", link = "https://doi.org/10.1038/s41467-022-30269-9",
         link_label = "Tesson et al. 2022, Nat Commun"),
    list(icon = "scissors", name = "CRISPR & Cas Systems (CRISPRCasFinder)",
         identifies = "CRISPR arrays (spacer \"memory\" of past phage encounters) and their Cas machinery -- shown as summary stats plus two browsable tables (arrays, Cas gene clusters).",
         use = "A record of phage exposure history plus the strain's adaptive defense system.",
         version = "CRISPRCasFinder 4.2.30", link = "https://doi.org/10.1093/nar/gky425",
         link_label = "Couvin et al. 2018, Nucleic Acids Res"),
    list(icon = "pills", name = "AMR & Specialty Genes (AMRFinderPlus)",
         identifies = "Antimicrobial resistance and virulence genes -- shown as a curated hits table (%identity, %coverage, reference accession) plus a keyword-screen fallback table for cross-checking.",
         use = "Flags antibiotic-resistance risk for a gut commensal strain.",
         version = "AMRFinderPlus (version not recorded in output)", link = "https://doi.org/10.1038/s41598-021-91456-0",
         link_label = "Feldgarden et al. 2021, Sci Rep")
  )),
  list(title = "Biosynthesis", status = "success", color = "#1b7a3d", entries = list(
    list(icon = "flask", name = "Biosynthetic Gene Clusters (antiSMASH)",
         identifies = "Clusters of genes that jointly produce a secondary metabolite (bacteriocins, siderophores, polyketides, etc.) -- shown as antiSMASH's own region-viewer style: a gene-arrow diagram per region, colored by gene role, with MIBiG known-cluster matches.",
         use = "Shows what specialized compounds a strain can make -- e.g. antimicrobials or signaling molecules that shape the gut community.",
         version = "antiSMASH 8.0.4", link = "https://doi.org/10.1093/nar/gkaf334",
         link_label = "Blin et al. 2025, Nucleic Acids Res")
  )),
  list(title = "Comparative & Evolutionary", status = "info", color = "#2c5f7c", entries = list(
    list(icon = "layer-group", name = "Pan-genome / Orthogroups",
         identifies = "Orthologous gene families shared across the 12 members -- shown as a stacked bar chart (core/soft-core/shell/unique) and a browsable orthogroup table with a gene tree per group.",
         use = "Shows what's core (shared by everyone) vs accessory (strain-specific).",
         version = "This app's own implementation (reciprocal-best-hit sequence search)", link = NULL),
    list(icon = "sitemap", name = "OMM12 Relatedness Tree (16S phylogeny)",
         identifies = "Evolutionary relationships among all 12 members -- shown as a real aligned 16S tree with bootstrap support, plus a separate pairwise proteome-similarity heatmap for comparison.",
         use = "How closely strains are related, which helps interpret shared functions.",
         version = "This app's own from-scratch pipeline: Gotoh affine-gap alignment, Jukes-Cantor (1969) distance, Neighbor-Joining, 200x bootstrap",
         link = "https://doi.org/10.1093/oxfordjournals.molbev.a040454", link_label = "Saitou & Nei 1987 (NJ method), Mol Biol Evol"),
    list(icon = "code-branch", name = "Phylogenetic Tree (per strain, KB18 only)",
         identifies = "A whole-genome + 16S rRNA GBDP phylogeny with branch support, for a single strain -- shown as an interactive tree on that strain's Bacteria Details page.",
         use = "The most rigorous single-strain tree in this app. Other strains show a note that none has been built yet, rather than a blank chart.",
         version = "TYGS (Type (Strain) Genome Server), tree search via FastME 2.1.4",
         link = "https://doi.org/10.1038/s41467-019-10210-3", link_label = "Meier-Kolthoff & Goker 2019, Nat Commun"),
    list(icon = "code-compare", name = "Compare Bacteria",
         identifies = "Side-by-side gene/trait differences between two chosen members -- shown as paired COG, localization, and mobile-element charts plus a shared-orthologs panel.",
         use = "Quickly spot functional differences between two strains.",
         version = "This app's own implementation", link = NULL),
    list(icon = "diagram-project", name = "Gene Neighborhood",
         identifies = "The genes flanking a gene of interest, across strains -- shown as stacked gggenes arrow diagrams, one row per strain.",
         use = "Checks whether local gene order (synteny/operon structure) is conserved.",
         version = "This app's own implementation (gggenes-based visualization)", link = NULL)
  )),
  list(title = "Practical Tools", status = "warning", color = "#8a5a12", entries = list(
    list(icon = "magnifying-glass", name = "Cross-Genome Search",
         identifies = "Which community members carry a gene, by symbol or product keyword (e.g. \"choloylglycine hydrolase\", the bile salt hydrolase) -- shown as a results table across all 12 genomes.",
         use = "Find a gene of interest community-wide in one keyword search.",
         version = "This app's own implementation", link = NULL),
    list(icon = "dna", name = "Sequence Search",
         identifies = "Matches for a pasted protein sequence against the OMM12 proteomes -- shown as a ranked hits table with alignment stats.",
         use = "BLAST-style lookup when you have an actual sequence, not just a gene name.",
         version = "This app's own implementation", link = NULL),
    list(icon = "ruler-horizontal", name = "Primer Design",
         identifies = "Candidate PCR primer pairs for a chosen gene -- shown as a table of forward/reverse pairs with product size, melting temp, and GC content.",
         use = "A practical wet-lab tool for strain-specific detection or qPCR.",
         version = "Native R reimplementation inspired by Primer3's approach (not the original Primer3 binary)",
         link = "https://doi.org/10.1093/nar/gks596", link_label = "Untergasser et al. 2012, Nucleic Acids Res")
  ))
)

tool_entry_ui <- function(e) {
  tags$div(
    style = "margin-bottom: 16px; padding-left: 4px; border-left: 3px solid #e3e7ec;",
    tags$div(style = "font-weight: 600; color: #1b263b; font-size: 14.5px; padding-left: 10px;",
             icon(e$icon), " ", e$name),
    tags$div(style = "color: #444; font-size: 13px; padding-left: 10px; margin-top: 3px;",
             tags$b("Identifies: "), e$identifies),
    tags$div(style = "color: #444; font-size: 13px; padding-left: 10px;",
             tags$b("Use: "), e$use),
    tags$div(style = "color: #777; font-size: 12px; padding-left: 10px; margin-top: 3px;",
             tags$b("Version: "), e$version %||% "unknown",
             if (!is.null(e$link)) {
               tagList(" -- ", tags$a(href = e$link, target = "_blank", rel = "noopener noreferrer",
                                       e$link_label %||% "reference"))
             } else NULL)
  )
}

render_tools_guide_page <- function() {
  tagList(lapply(tools_guide_categories, function(cat) {
    fluidRow(
      box(title = cat$title, width = 12, solidHeader = TRUE, status = cat$status,
          lapply(cat$entries, tool_entry_ui))
    )
  }))
}

# One-line blurb (identify + use) shown at the top of a section, linking back
# to the full Tools Guide page. `tool_name` must match an entry's `name`
# above (across any category) or this throws early, at app start, rather
# than silently showing nothing.
section_blurb <- function(tool_name) {
  entry <- NULL
  for (cat in tools_guide_categories) {
    for (e in cat$entries) {
      if (e$name == tool_name) { entry <- e; break }
    }
    if (!is.null(entry)) break
  }
  if (is.null(entry)) stop("section_blurb: no Tools Guide entry named '", tool_name, "'")
  # Plain anchor + inline onclick (not a Shiny actionLink) so this can be
  # reused on many pages at once without colliding Shiny input IDs -- it just
  # clicks the "Tools Guide" sidebar entry client-side, no server round-trip.
  p(style = "color:#555; font-size: 13px; margin-bottom: 14px;",
    icon(entry$icon), " ", tags$b("Identifies: "), entry$identifies, " ",
    tags$b("Use: "), entry$use, " ",
    tags$a(href = "#", "More in Tools Guide →",
           style = "font-size: 12.5px; margin-left: 4px;",
           onclick = "$('a[data-value=\"tools_guide\"]').click(); return false;"))
}

# --- UI ---
# --- DATA STATUS TABLE (About page) ---
# Built at app startup from which files exist in data/, so it updates itself
# whenever new tool output is added and the app/image is rebuilt.
data_status_ui <- function() {
  checks <- list(
    "Annotation" = "_prokka.gff", "COG" = "_cog.tsv", "KEGG" = "_kegg.tsv",
    "eggNOG" = "_eggnog.gff", "Localization" = "_localization.csv",
    "MobileOG" = "_mobileog.csv",
    "Defense" = c("_defense_systems.tsv", "_defense_finder_systems.tsv"),
    "CRISPR" = "_crispr_arrays.tsv", "BGC" = "_antismash_regions.tsv",
    "AMR" = "_amrfinder.tsv", "Metabolic model" = "_model.rds",
    "Genome tree" = "_tree.phy"
  )
  has <- function(bf, sfx) any(file.exists(file.path("data", paste0(bf, sfx))))
  header <- tags$tr(
    tags$th("Strain"),
    unname(lapply(names(checks), function(n) tags$th(n, style = "text-align:center;")))
  )
  rows <- lapply(seq_along(bacteria_filenames), function(i) {
    bf <- bacteria_filenames[i]
    tags$tr(
      tags$td(tags$em(bacteria_names[i])),
      unname(lapply(checks, function(sfx) {
        ok <- has(bf, sfx)
        tags$td(if (ok) "\u2713" else "\u2013",
                style = paste0("text-align:center;color:", if (ok) "#2e7d32" else "#bbbbbb", ";"))
      }))
    )
  })
  div(style = "overflow-x:auto;",
      tags$table(class = "table table-condensed table-bordered", style = "font-size:12px;",
                 tags$thead(header), tags$tbody(rows)))
}

ui <- dashboardPage(
  skin = "black-light",
  
  dashboardHeader(
    title = tags$a(
      href = "#", style = "display:flex; align-items:center; gap:8px;",
      tags$span("OMM12 Resource")
    ),
    titleWidth = 250
  ),
  
  dashboardSidebar(
    width = 250,
    sidebarMenu(
      id = "sidebar_menu",
      menuItem("About", tabName = "about", icon = icon("info-circle")),
      menuItem("Tools Guide", tabName = "tools_guide", icon = icon("circle-question")),
      menuItem("OMM12 Resource", tabName = "omm12", icon = icon("dna")),
      menuItem("Cross-Genome Search", tabName = "cross_search", icon = icon("magnifying-glass")),
      menuItem("Compare Bacteria", tabName = "compare_bacteria", icon = icon("code-compare")),
      menuItem("Gene Neighborhood", tabName = "gene_neighborhood", icon = icon("diagram-project")),
      menuItem("Pan-genome", tabName = "pan_genome", icon = icon("layer-group")),
      menuItem("OMM12 Relatedness Tree", tabName = "omm12_relatedness", icon = icon("sitemap")),
      menuItem("Sequence Search", tabName = "sequence_search", icon = icon("dna")),
      menuItem("Primer Design", tabName = "primer_design", icon = icon("ruler-horizontal")),
      menuItem("Defense Systems", tabName = "defense_overview", icon = icon("shield-halved")),
      menuItem("CRISPR & Cas Systems", tabName = "crispr_systems", icon = icon("scissors")),
      menuItem("Biosynthetic Gene Clusters", tabName = "bgc_clusters", icon = icon("flask")),
      menuItem("AMR & Specialty Genes", tabName = "specialty_genes", icon = icon("pills")),
      shinyjs::hidden(
        menuItem("Bacteria Details", tabName = "details_page", icon = icon("microscope"))
      )
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$title("OMM12 Resource — Genomic & Functional Browser"),
      tags$meta(name = "description", content = paste(
        "Genomic, functional, and phylogenetic browser for the 12 bacterial strains of the",
        "Oligo-Mouse-Microbiota-12 (OMM12) synthetic gut community.")),
      tags$link(rel = "icon", type = "image/x-icon", href = "favicon.ico"),
      tags$link(rel = "icon", type = "image/png", sizes = "32x32", href = "favicon.png"),
      tags$link(rel = "stylesheet", type = "text/css",href = "custom.css?v=5"),
      tags$link(rel = "stylesheet", href = "https://cdn.jsdelivr.net/npm/cgview/dist/cgview.css"),
      tags$script(src = "https://cdn.jsdelivr.net/npm/d3@7"),
      tags$script(src = "https://cdn.jsdelivr.net/npm/cgview/dist/cgview.min.js"),
      tags$script(src = "https://cdn.jsdelivr.net/npm/cgparse/dist/cgparse.min.js"),
      js_scroll_code
    ),
    useShinyjs(),
    
    tabItems(
      tabItem(
        tabName = "about",
        h2("About This Resource"),
        fluidRow(
          box(
            title = "About This Resource", width = 12, solidHeader = TRUE, status = "primary",
            p("This resource provides detailed genomic, functional, and phylogenetic information for the 12 ",
              "bacterial strains of the Oligo-Mouse-Microbiota-12 (OMM12) synthetic gut community: annotated ",
              "genomes, COG and KEGG functional profiles, predicted subcellular localization, mobile genetic ",
              "elements, phylogenetic relationships, and — where available — genome-scale metabolic models.")
          )
        ),
        fluidRow(
          box(
            title = "What is OMM12?", width = 12, solidHeader = TRUE, status = "info",
            p("OMM12 (Oligo-Mouse-Microbiota-12) is a synthetic, defined bacterial community of 12 strains ",
              "representing the major bacterial phyla of the mouse gut microbiota. It was designed by ",
              tags$a(href = "https://www.nature.com/articles/nmicrobiol2016215", target = "_blank",
                     "Brugiroux et al."),
              " and first described in \"Genome-guided design of a defined mouse microbiota that confers ",
              "colonization resistance against Salmonella enterica serovar Typhimurium\" (Nature Microbiology, ",
              "2016; DOI: ",
              tags$a(href = "https://doi.org/10.1038/nmicrobiol.2016.215", target = "_blank",
                     "10.1038/nmicrobiol.2016.215"), ")."),
            p("Once established in germ-free mice, OMM12 is stable across consecutive mouse generations and ",
              "reproducible across different gnotobiotic facilities, and it confers colonization resistance ",
              "against enteric pathogens such as ", tags$em("Salmonella"), " Typhimurium. That combination of ",
              "phylogenetic breadth, reproducibility, and public strain availability (all 12 strains are ",
              "deposited at the ",
              tags$a(href = "https://www.dsmz.de/collection/catalogue/microorganisms/special-groups-of-organisms/dzif-sammlung/maus-mikrobiomliste/oligo-mouse-microbiota",
                     target = "_blank", "DSMZ culture collection"),
              ") has made it a widely used, tractable model for studying host-microbiome interactions under ",
              "fully controlled conditions."),
            tags$h4("The 12 community members", style = "margin-top: 20px;"),
            tags$ol(lapply(bacteria_names, function(n) tags$li(n)))
          )
        ),
        fluidRow(
          box(
            title = "Publication Statistics", width = 12, solidHeader = TRUE, status = "warning",
            uiOutput("pub_caption"),
            uiOutput("pub_stats_summary"),
            fluidRow(
              column(width = 6, plotlyOutput("pub_year_chart", height = "360px")),
              column(width = 6, plotlyOutput("pub_citations_chart", height = "360px"))
            ),
            tags$hr(style = "border-top: 1px solid #eee;"),
            radioButtons("pub_map_group", NULL, inline = TRUE,
                         choices = c("All publications" = "all",
                                     "Only publications that use OMM12" = "used")),
            plotlyOutput("pub_world_map", height = "460px"),
            tags$p(style = "color:#888; font-size: 12px;",
                   "Countries from the authors' affiliations; a publication with authors from several ",
                   "countries counts once for each country."),
            tags$hr(style = "border-top: 1px solid #eee;"),
            plotlyOutput("pub_journals_chart", height = "420px")
          )
        )
      ),

      tabItem(
        tabName = "tools_guide",
        h2("Tools Guide"),
        fluidRow(
          box(width = 12, solidHeader = TRUE, status = "primary",
              p(style = "color:#555;",
                "What each section of this resource identifies, and what it's useful for -- grouped by category. ",
                "Every section also has a short version of this at the top linking back here."))
        ),
        render_tools_guide_page()
      ),

      tabItem(
        tabName = "omm12",
        h2("OMM12 Bacterial Collection"),
        fluidRow(
          lapply(1:4, function(i) {
            column(
              width = 3, class = "bacteria-item",
              div(class = "bacteria-item-content",
                  tags$img(src = paste0(bacteria_filenames[i], ".png"),
                           alt = bacteria_names[i],
                           id = paste0("bacteria", i, "-img"),
                           class = "bacteria-image"),
                  tags$p(bacteria_names[i], class = "bacteria-name")
              )
            )
          })
        ),
        fluidRow(
          lapply(5:8, function(i) {
            column(
              width = 3, class = "bacteria-item",
              div(class = "bacteria-item-content",
                  tags$img(src = paste0(bacteria_filenames[i], ".png"),
                           alt = bacteria_names[i],
                           id = paste0("bacteria", i, "-img"),
                           class = "bacteria-image"),
                  tags$p(bacteria_names[i], class = "bacteria-name")
              )
            )
          })
        ),
        fluidRow(
          lapply(9:12, function(i) {
            column(
              width = 3, class = "bacteria-item",
              div(class = "bacteria-item-content",
                  tags$img(src = paste0(bacteria_filenames[i], ".png"),
                           alt = bacteria_names[i],
                           id = paste0("bacteria", i, "-img"),
                           class = "bacteria-image"),
                  tags$p(bacteria_names[i], class = "bacteria-name")
              )
            )
          })
        )
      ),
      
      tabItem(
        tabName = "cross_search",
        h2("Cross-Genome Gene Search"),
        section_blurb("Cross-Genome Search"),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "info",
            title = "Search across all 12 OMM12 community members",
            p("Search by gene symbol or a keyword in the product description (e.g. \"choloylglycine hydrolase\" (bile salt hydrolase), ",
              "\"flagellin\", \"ABC transporter\") to see which community members carry it, and where."),
            fluidRow(
              column(width = 6,
                     textInput("cross_search_query", NULL,
                               placeholder = "e.g. recA, cbh, choloylglycine hydrolase, ABC transporter...",
                               width = "100%")
              ),
              column(width = 4,
                     selectInput("cross_search_type", NULL,
                                 choices = c("Gene symbol or product" = "any",
                                             "Gene symbol only" = "gene",
                                             "Product/description only" = "product"),
                                 selected = "any", width = "100%")
              )
            ),
            uiOutput("cross_search_summary")
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "primary", title = "Results",
            p(style = "color:#666; font-size: 12px;",
              "Click a row to jump to that bacterium's Genome Map."),
            DTOutput("cross_search_results")
          )
        )
      ),

      tabItem(
        tabName = "compare_bacteria",
        h2("Compare Two Bacteria"),
        section_blurb("Compare Bacteria"),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "info",
            title = "Select two community members",
            fluidRow(
              column(width = 5,
                     selectInput("compare_bact_a", "Bacterium A",
                                 choices = setNames(bacteria_filenames, bacteria_names),
                                 selected = bacteria_filenames[1], width = "100%")
              ),
              column(width = 2,
                     div(style = "text-align:center; padding-top: 32px; font-weight:600; color:#888;", "vs")
              ),
              column(width = 5,
                     selectInput("compare_bact_b", "Bacterium B",
                                 choices = setNames(bacteria_filenames, bacteria_names),
                                 selected = bacteria_filenames[2], width = "100%")
              )
            )
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "primary", title = "Genome Overview",
            DTOutput("compare_genome_stats_table")
          )
        ),
        fluidRow(
          box(width = 6, solidHeader = TRUE, status = "info", title = "COG Functional Profile",
              plotlyOutput("compare_cog_chart", height = "380px")),
          box(width = 6, solidHeader = TRUE, status = "info", title = "Subcellular Localization Profile",
              plotlyOutput("compare_localization_chart", height = "380px"))
        ),
        fluidRow(
          box(width = 6, solidHeader = TRUE, status = "warning", title = "Mobile Genetic Element Categories",
              plotlyOutput("compare_mobileog_chart", height = "380px")),
          box(width = 6, solidHeader = TRUE, status = "warning", title = "Shared Genes (orthologs)",
              p(style = "color:#666; font-size: 12px;",
                "Reciprocal-best-hit orthologs found by comparing actual protein sequences (k-mer similarity) ",
                "between the two genomes, not by matching gene-symbol text. A lightweight approximation of true ",
                "orthology, so it will miss some genuine but highly diverged orthologs."),
              uiOutput("compare_shared_summary"),
              DTOutput("compare_shared_genes_table")
          )
        )
      ),

      tabItem(
        tabName = "gene_neighborhood",
        h2("Gene Neighborhood / Synteny View"),
        section_blurb("Gene Neighborhood"),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "info",
            title = "Compare a gene's genomic neighborhood across OMM12 members",
            p(style = "color:#666; font-size: 13px;",
              "Pick a gene in one bacterium (the anchor), and see the genes immediately upstream and downstream ",
              "of it — and of its ortholog — in every other OMM12 member. Orthologs are found by comparing actual ",
              "protein sequences (k-mer-based similarity with reciprocal-best-hit checking across all 12 ",
              "proteomes), not by matching gene-symbol text. Neighbors are aligned by rank position around the ",
              "anchor (not physical distance), similar to MicrobesOnline's neighborhood view. This is a ",
              "lightweight approximation of true orthology — good enough to be directionally useful, but not a ",
              "substitute for a proper OrthoFinder/DIAMOND-based analysis."),
            fluidRow(
              column(width = 4,
                     selectInput("neighborhood_bact", "Anchor bacterium",
                                 choices = setNames(bacteria_filenames, bacteria_names),
                                 selected = bacteria_filenames[1], width = "100%")
              ),
              column(width = 4,
                     selectizeInput("neighborhood_gene", "Anchor gene", choices = NULL,
                                    options = list(placeholder = "Type a gene symbol…", maxOptions = 50),
                                    width = "100%")
              ),
              column(width = 3,
                     numericInput("neighborhood_window", "Genes each side", value = 3, min = 1, max = 8, step = 1)
              )
            )
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "primary", title = "Neighborhood",
            uiOutput("neighborhood_summary"),
            shinycssloaders::withSpinner(uiOutput("neighborhood_plot_ui"), type = 6, color = "#1b263b")
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "warning", title = "Gene details",
            p(style = "color:#666; font-size: 12px;",
              "Every gene shown in the diagram above, with its product and ortholog status."),
            DTOutput("neighborhood_genes_table")
          )
        )
      ),

      tabItem(
        tabName = "pan_genome",
        h2("Pan-genome / Ortholog Matrix"),
        section_blurb("Pan-genome / Orthogroups"),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "info", title = "Core, soft-core, shell & unique genes",
            p(style = "color:#666; font-size: 13px;",
              "Orthogroups are connected components of the reciprocal-best-hit graph built from real ",
              "protein-sequence comparison across all 12 proteomes (same method as the Gene Neighborhood and ",
              "Compare Bacteria pages) — a lightweight stand-in for a proper OrthoFinder run, not a replacement ",
              "for one. It will undercount true core genes between highly divergent members, so treat these ",
              "numbers as directional rather than definitive. ",
              tags$b("Core"), " = present in all 12 members, ", tags$b("soft-core"), " = 9–11, ",
              tags$b("shell"), " = 2–8, ", tags$b("unique"), " = found in only one member (by this method)."),
            uiOutput("pangenome_summary_stats")
          )
        ),
        fluidRow(
          box(width = 12, solidHeader = TRUE, status = "primary",
              title = "Genome composition by conservation level",
              plotlyOutput("pangenome_composition_chart", height = "420px"))
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "warning", title = "Browse orthogroups",
            fluidRow(
              column(width = 4,
                     selectInput("pangenome_class_filter", "Classification",
                                 choices = c("All" = "all", "Core" = "core", "Soft-core" = "soft-core",
                                             "Shell" = "shell", "Unique" = "unique"),
                                 selected = "core", width = "100%")
              )
            ),
            p(style = "color:#666; font-size: 12px;", "Click a row to see its full gene membership below."),
            DTOutput("pangenome_orthogroup_table")
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "warning",
            title = uiOutput("pangenome_selected_og_title"),
            DTOutput("pangenome_orthogroup_members_table")
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "warning",
            title = "Gene tree for this orthogroup",
            p(style = "color:#666; font-size: 12px;",
              "A quick UPGMA dendrogram of just this gene family's members, from live pairwise k-mer ",
              "similarity (same method as Sequence Search) — a BV-BRC-style \"Gene Tree\" scoped to one ",
              "orthogroup, not a validated phylogeny. Needs an orthogroup with 3+ members selected above."),
            shinycssloaders::withSpinner(uiOutput("pangenome_orthogroup_tree_ui"), type = 6, color = "#1b263b")
          )
        )
      ),

      tabItem(
        tabName = "omm12_relatedness",
        h2("OMM12 Relatedness Tree"),
        section_blurb("OMM12 Relatedness Tree (16S phylogeny)"),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "info",
            title = "How this differs from the Phylogeny tab",
            p(style = "color:#666; font-size: 13px;",
              "The genome-wide GBDP/16S phylogenetic tree (Phylogeny tab, under a bacterium's ",
              "details page) currently only exists for the reference strain, Acutalibacter muris ",
              "KB18. This page covers all 12 members instead, with a real single-locus (16S rRNA ",
              "gene) phylogeny: one representative 16S copy per genome, aligned with a custom ",
              "affine-gap-penalty multiple alignment, Jukes-Cantor-corrected pairwise distances, a ",
              "Neighbor-Joining tree, and support from 200 bootstrap replicates. ",
              tags$b("This is a real alignment + substitution model + bootstrap tree"),
              ", not a rough sketch — but it's still a single gene, so it has limited power to ",
              "resolve the deepest splits between very different phyla; treat nodes with bootstrap ",
              "support well under 70% with appropriate caution. Downloads: ",
              tags$a(href = "data/omm12_16s_raw.fasta", target = "_blank", "16S sequences (FASTA)"), ", ",
              tags$a(href = "data/omm12_16s_msa.fasta", target = "_blank", "alignment (FASTA)"), ", ",
              tags$a(href = "data/omm12_16s_tree.nwk", target = "_blank", "tree (Newick)"), ", ",
              tags$a(href = "data/omm12_16s_pipeline_report.txt", target = "_blank", "methods report"),
              ". The heatmap alongside the tree is a ", tags$em("different"),
              " signal — whole-proteome k-mer similarity (AAI-like), not a phylogeny — kept as a ",
              "complementary view.")
          )
        ),
        fluidRow(
          box(width = 6, solidHeader = TRUE, status = "primary", title = "16S rRNA gene phylogeny",
              shinycssloaders::withSpinner(girafeOutput("omm12_relatedness_tree", height = "500px"),
                                            type = 6, color = "#1b263b")),
          box(width = 6, solidHeader = TRUE, status = "primary", title = "Pairwise proteome similarity (k-mer, approximate)",
              p(style = "color:#666; font-size: 12px;",
                "Mean k-mer containment score across reciprocal-best-hit orthologs between each pair — ",
                "a different, whole-proteome signal, not a phylogeny."),
              shinycssloaders::withSpinner(plotlyOutput("omm12_relatedness_heatmap", height = "460px"),
                                            type = 6, color = "#1b263b"))
        )
      ),

      tabItem(
        tabName = "sequence_search",
        h2("Sequence Search"),
        section_blurb("Sequence Search"),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "info",
            title = "Search a protein sequence against the OMM12 proteomes",
            p(style = "color:#666; font-size: 13px;",
              "Paste a protein sequence (raw residues, or FASTA with a ", tags$code(">header"),
              " line — the header is stripped automatically) and search it against one or all 12 ",
              "OMM12 proteomes, similar to MicrobesOnline's Sequence Search. Matches use the same ",
              "k-mer containment method validated for the ortholog table elsewhere in this app (k=4 ",
              "amino-acid k-mers, score = shared k-mers / smaller sequence's k-mer count). ",
              tags$b("This is a lightweight approximation, not a real BLAST/DIAMOND alignment"),
              " — there's no gapped alignment, E-value, percent identity, or coverage statistic, ",
              "so treat results as a starting point for finding likely homologs, not a definitive call."),
            fluidRow(
              column(width = 12,
                     textAreaInput("seq_search_query", "Query sequence (protein, one-letter amino acid code)",
                                   value = "", rows = 6, resize = "vertical", width = "100%",
                                   placeholder = "Paste a protein sequence here, e.g.\n>my_query\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV...")
              )
            ),
            fluidRow(
              column(width = 4,
                     selectInput("seq_search_target", "Search against",
                                 choices = c("All 12 bacteria" = "all",
                                             setNames(bacteria_filenames, bacteria_names)),
                                 selected = "all", width = "100%")
              ),
              column(width = 3,
                     numericInput("seq_search_top_n", "Max hits to show", value = 25, min = 1, max = 200, step = 1)
              ),
              column(width = 3, style = "padding-top: 25px;",
                     actionButton("seq_search_run", "Search", icon = icon("magnifying-glass"),
                                  class = "btn-primary", width = "100%")
              )
            )
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "primary", title = "Results",
            shinycssloaders::withSpinner(uiOutput("seq_search_summary"), type = 6, color = "#1b263b"),
            DTOutput("seq_search_results")
          )
        )
      ),

      tabItem(
        tabName = "primer_design",
        h2("Primer Design"),
        section_blurb("Primer Design"),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "info",
            title = "Design candidate PCR primers for a gene",
            p(style = "color:#666; font-size: 13px;",
              "Pick a gene from any OMM12 member (nucleotide sequence pulled from its genome using the ",
              "Prokka-annotated coordinates), or paste your own sequence, and get a ranked shortlist of ",
              "candidate forward/reverse primer pairs — similar in spirit to BV-BRC's Primer Design ",
              "Service. ", tags$b("This uses simplified rules of thumb, not full Primer3-style ",
              "thermodynamics:"), " Tm is estimated from a standard empirical GC-content formula (not ",
              "nearest-neighbor calculation), and hairpin/self-complementarity checks are naive pattern ",
              "matches, not real secondary-structure prediction. There's also no forward/reverse primer-",
              "dimer cross-check. Verify any pair with a dedicated tool (Primer3, IDT OligoAnalyzer) ",
              "before ordering.")
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "primary", title = "Target sequence",
            radioButtons("primer_source", NULL,
                         choices = c("Pick a gene from an OMM12 member" = "gene",
                                     "Paste my own sequence" = "custom"),
                         selected = "gene", inline = TRUE),
            conditionalPanel(
              condition = "input.primer_source == 'gene'",
              fluidRow(
                column(width = 4,
                       selectInput("primer_bact", "Bacterium",
                                   choices = setNames(bacteria_filenames, bacteria_names),
                                   selected = bacteria_filenames[1], width = "100%")),
                column(width = 4,
                       selectizeInput("primer_gene", "Gene", choices = NULL,
                                      options = list(placeholder = "Type a gene symbol…", maxOptions = 50),
                                      width = "100%"))
              ),
              uiOutput("primer_gene_length_note")
            ),
            conditionalPanel(
              condition = "input.primer_source == 'custom'",
              textAreaInput("primer_custom_seq", "Nucleotide sequence (DNA, one-letter code)",
                            value = "", rows = 5, resize = "vertical", width = "100%",
                            placeholder = "Paste a DNA sequence here, raw or FASTA with a >header line…")
            )
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "warning", title = "Design parameters",
            fluidRow(
              column(width = 2, numericInput("primer_len_min", "Min length", value = 18, min = 12, max = 30)),
              column(width = 2, numericInput("primer_len_max", "Max length", value = 25, min = 12, max = 35)),
              column(width = 2, numericInput("primer_tm_target", "Target Tm (C)", value = 60, min = 45, max = 75)),
              column(width = 2, numericInput("primer_gc_min", "Min GC%", value = 40, min = 0, max = 100)),
              column(width = 2, numericInput("primer_gc_max", "Max GC%", value = 60, min = 0, max = 100))
            ),
            fluidRow(
              column(width = 2, numericInput("primer_product_min", "Min product (bp)", value = 100, min = 40)),
              column(width = 2, numericInput("primer_product_max", "Max product (bp)", value = 800, min = 60)),
              column(width = 2, numericInput("primer_n_pairs", "Pairs to show", value = 5, min = 1, max = 20)),
              column(width = 3, style = "padding-top: 25px;",
                     actionButton("primer_design_run", "Design Primers", icon = icon("dna"),
                                  class = "btn-primary", width = "100%"))
            )
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "primary", title = "Candidate primer pairs",
            shinycssloaders::withSpinner(uiOutput("primer_design_summary"), type = 6, color = "#1b263b"),
            DTOutput("primer_design_results")
          )
        )
      ),

      tabItem(
        tabName = "defense_overview",
        h2("Defense Systems Overview"),
        section_blurb("Defense Systems (DefenseFinder)"),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "info",
            title = "Anti-phage / anti-MGE defense systems across OMM12",
            p(style = "color:#666; font-size: 13px;",
              "Systems predicted with ", tags$a(href = "https://github.com/mdmparis/defense-finder",
                                                 target = "_blank", "DefenseFinder"),
              " for each OMM12 member. ", tags$b("Defense"), " systems (restriction-modification, ",
              "CBASS, Wadjet, CRISPR-Cas, BREX, and others) protect the host against phages and ",
              "other mobile genetic elements. ", tags$b("Antidefense"), " systems (Anti-CRISPR, ",
              "Anti-RM, Anti-RecBCD) counteract host defenses and are themselves often carried on ",
              "prophages, so their presence is itself a signal of past or ongoing phage interaction."),
            uiOutput("defense_overview_missing_note"),
            uiOutput("defense_overview_summary_stats")
          )
        ),
        fluidRow(
          box(width = 7, solidHeader = TRUE, status = "primary",
              title = "Defense system types by bacterium",
              shinycssloaders::withSpinner(plotlyOutput("defense_overview_heatmap", height = "440px"),
                                            type = 6, color = "#1b263b")),
          box(width = 5, solidHeader = TRUE, status = "primary",
              title = "Defense vs. antidefense systems per bacterium",
              shinycssloaders::withSpinner(plotlyOutput("defense_overview_bar", height = "440px"),
                                            type = 6, color = "#1b263b"))
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "warning", title = "Browse all defense systems",
            fluidRow(
              column(width = 4,
                     selectInput("defense_overview_activity_filter", "Activity",
                                 choices = c("All" = "all", "Defense" = "Defense", "Antidefense" = "Antidefense"),
                                 selected = "all", width = "100%")
              ),
              column(width = 4,
                     selectInput("defense_overview_bact_filter", "Bacterium",
                                 choices = c("All 12 bacteria" = "all", setNames(bacteria_filenames, bacteria_names)),
                                 selected = "all", width = "100%")
              )
            ),
            DTOutput("defense_overview_table")
          )
        )
      ),

      tabItem(
        tabName = "crispr_systems",
        h2("CRISPR Arrays & Cas Systems"),
        section_blurb("CRISPR & Cas Systems (CRISPRCasFinder)"),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "info",
            title = "CRISPR immunity across OMM12",
            p(style = "color:#666; font-size: 13px;",
              "Arrays and Cas gene clusters predicted with ",
              tags$a(href = "https://crisprcas.i2bc.paris-saclay.fr/", target = "_blank", "CRISPRCasFinder"),
              " (version 4.2.30) for each OMM12 member. Each array gets CRISPRCasFinder's own ",
              tags$b("evidence level"), " from 1 (weak/possible — often a single short, low-conservation ",
              "repeat) to 4 (strong — multiple well-conserved repeats), shown as-is rather than filtered, ",
              "so treat level-1 calls with more caution than level-4 ones. A CRISPR array and a nearby Cas ",
              "gene cluster aren't guaranteed to co-occur — some bacteria carry ", tags$b("orphan arrays"),
              " with no adjacent Cas machinery (the interference genes are elsewhere, degraded, or absent), ",
              "and CRISPRCasFinder tracks the two independently, so this page does too."),
            uiOutput("crispr_missing_note"),
            uiOutput("crispr_summary_stats")
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "warning", title = "Browse CRISPR arrays",
            fluidRow(
              column(width = 4,
                     selectInput("crispr_bact_filter", "Bacterium",
                                 choices = c("All bacteria with data" = "all", setNames(bacteria_filenames, bacteria_names)),
                                 selected = "all", width = "100%"))
            ),
            DTOutput("crispr_arrays_table")
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "warning", title = "Browse Cas gene clusters",
            fluidRow(
              column(width = 4,
                     selectInput("cas_bact_filter", "Bacterium",
                                 choices = c("All bacteria with data" = "all", setNames(bacteria_filenames, bacteria_names)),
                                 selected = "all", width = "100%"))
            ),
            DTOutput("cas_systems_table")
          )
        )
      ),

      tabItem(
        tabName = "bgc_clusters",
        h2("Biosynthetic Gene Clusters (BGC)"),
        section_blurb("Biosynthetic Gene Clusters (antiSMASH)"),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "info",
            title = "Secondary-metabolite potential across OMM12",
            p(style = "color:#666; font-size: 13px;",
              "Regions predicted with ", tags$a(href = "https://antismash.secondarymetabolites.org/",
                                                  target = "_blank", "antiSMASH"),
              " (version 8.0.4) for each OMM12 member. ", tags$b("Predicted product"), " is antiSMASH's ",
              "rule-based call for what kind of gene cluster this is (RiPP, terpene, PKS, etc.) — a ",
              "structural prediction, not confirmation the compound is actually made or biologically active. ",
              "When antiSMASH's knownclusterblast module finds a resembling entry in the ", tags$b("MIBiG"),
              " database of experimentally characterized clusters, its accession, description, and percent ",
              "similarity are shown too — but a low percentage (well under ~50%) usually just means a few ",
              "genes overlap by chance (e.g. shared housekeeping genes), not a real match to that specific ",
              "known product, so read those as ", tags$em("unclassified/novel"), " rather than identified."),
            uiOutput("antismash_missing_note"),
            uiOutput("antismash_summary_stats")
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "warning", title = "Region viewer",
            p(style = "color:#666; font-size: 12px;",
              "Same style as antiSMASH's own region view: gene arrows colored by function, drawn to scale ",
              "along the region. Pick a bacterium and one of its predicted regions."),
            fluidRow(
              column(width = 4,
                     selectInput("bgc_viewer_bact", "Bacterium",
                                 choices = setNames(antismash_bacteria_with_data,
                                                     bacteria_names[match(antismash_bacteria_with_data, bacteria_filenames)]),
                                 width = "100%")),
              column(width = 5,
                     selectInput("bgc_viewer_region", "Region", choices = NULL, width = "100%"))
            ),
            uiOutput("bgc_viewer_header"),
            shinycssloaders::withSpinner(plotOutput("bgc_viewer_plot", height = "260px"),
                                          type = 6, color = "#1b263b"),
            uiOutput("bgc_viewer_mibig")
          )
        )
      ),

      tabItem(
        tabName = "specialty_genes",
        h2("AMR & Specialty Genes"),
        section_blurb("AMR & Specialty Genes (AMRFinderPlus)"),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "info",
            title = "How this screen works",
            p(style = "color:#666; font-size: 13px;",
              "Real AMRFinderPlus results are being added strain by strain (see the curated section ",
              "below); CARD-RGI, ResFinder, and VFDB runs aren't available yet. For any bacterium ",
              "without a curated run, this is a ", tags$b("transparent keyword screen"),
              " against the gene symbols and product descriptions Prokka already assigned (from its ",
              "own UniProt/Swiss-Prot reference databases) — flagging genes whose name or product text ",
              "matches well-known antimicrobial-resistance gene families (", tags$em("bla, tet, erm, van,",
              " aac/aph/ant, cat, sul/dfr, mcr"), ", multidrug efflux) or virulence-associated terms ",
              "(toxins, adhesins, secretion systems, siderophores). ",
              tags$b("This is not equivalent to a curated CARD/ResFinder/VFDB alignment call"),
              " — no percent identity, coverage, or E-value, and no result at all for a real resistance ",
              "gene Prokka happened to annotate as \"hypothetical protein\". Treat hits as leads worth ",
              "confirming with a dedicated tool, not a validated resistance determination — the same ",
              "spirit as BacDive's own \"text-mining / automatically generated\" data flag.")
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "success",
            title = "Curated AMR calls (AMRFinderPlus)",
            uiOutput("amrfinder_missing_note"),
            uiOutput("amrfinder_summary_stats")
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "success", title = "Browse curated AMR calls",
            fluidRow(
              column(width = 4,
                     selectInput("amrfinder_bact_filter", "Bacterium",
                                 choices = c("All bacteria with data" = "all", setNames(bacteria_filenames, bacteria_names)),
                                 selected = "all", width = "100%"))
            ),
            DTOutput("amrfinder_table")
          )
        ),
        fluidRow(
          box(width = 12, solidHeader = TRUE, status = "primary",
              title = "Keyword screen summary (all 12 bacteria, approximate)",
              uiOutput("specialty_genes_summary_stats"))
        ),
        fluidRow(
          box(width = 12, solidHeader = TRUE, status = "primary",
              title = "Keyword screen: hits by bacterium and subcategory",
              shinycssloaders::withSpinner(plotlyOutput("specialty_genes_heatmap", height = "440px"),
                                            type = 6, color = "#1b263b"))
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "warning", title = "Browse keyword-screen hits",
            fluidRow(
              column(width = 4,
                     selectInput("specialty_genes_category_filter", "Category",
                                 choices = c("All" = "all", "AMR" = "AMR", "Virulence" = "Virulence"),
                                 selected = "all", width = "100%")),
              column(width = 4,
                     selectInput("specialty_genes_bact_filter", "Bacterium",
                                 choices = c("All 12 bacteria" = "all", setNames(bacteria_filenames, bacteria_names)),
                                 selected = "all", width = "100%"))
            ),
            DTOutput("specialty_genes_table")
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = TRUE, status = "warning", title = "AMR phenotype (lab-tested)",
            p(style = "color:#666; font-size: 13px;",
              "Genome-based screens above can only ever suggest candidate resistance genes — actual ",
              "susceptibility (MIC or disk-diffusion results) has to come from real testing, not ",
              "prediction, so this section only ever shows real data, never a guess."),
            fluidRow(
              column(width = 4,
                     selectInput("amr_phenotype_bact", "Bacterium",
                                 choices = setNames(bacteria_filenames, bacteria_names),
                                 selected = bacteria_filenames[1], width = "100%"))
            ),
            uiOutput("amr_phenotype_content")
          )
        )
      ),

      tabItem(tabName = "details_page",
              fluidRow(
                column(2, id = "navigation_panel_column", style = "position: sticky; top: 50px; height: calc(100vh - 50px); overflow-y: auto; padding-top: 15px;",
                       
                       tags$div(
                         h4(icon("list-alt"), "Content Index", style = "border-bottom: 1px solid #eee; padding-bottom: 5px;"),
                         tags$ul(class = "list-unstyled", 
                                 tags$li(actionLink("nav_genomic_info", tags$strong("Genomic Information"), style = "font-size: 16px;")),
                                 tags$ul(class = "list-unstyled", style = "padding-left: 15px;",
                                         tags$li(actionLink("nav_phylogeny", "Phylogeny")),
                                         tags$li(actionLink("nav_genome_map", "Genome Map")),
                                         tags$li(actionLink("nav_cog", "COG")),
                                         tags$li(actionLink("nav_mobileog", "Bacterial Mobile Genetic Elements")),
                                         tags$li(actionLink("nav_kegg", "KEGG Pathways")),
                                         tags$li(actionLink("nav_localization", "Localization")),
                                         tags$li(actionLink("nav_defense", "Defense Systems")),
                                         tags$li(actionLink("nav_crispr", "CRISPR & Cas Systems")),
                                         tags$li(actionLink("nav_bgc", "Biosynthetic Gene Clusters")),
                                         tags$li(actionLink("nav_specialty_genes", "AMR & Specialty Genes"))
                                 ),
                                 tags$li(tags$br()),
                                 #tags$li(actionLink("nav_proteomic_info", tags$strong("Proteomic Information"), style = "font-size: 16px;")),
                                 tags$ul(class = "list-unstyled", style = "padding-left: 15px;"
                                 ),
                                 tags$li(tags$br()),
                                 tags$li(actionLink("nav_metabolic_info", tags$strong("Metabolic Information"), style = "font-size: 16px;")),
                                 tags$ul(class = "list-unstyled", style = "padding-left: 15px;",
                                         tags$li(actionLink("nav_metabolic_model", "Metabolic Model"))
                         )
                       )
                       )
                ),
                
                column(10, uiOutput("bacteria_details_ui"))
              ),
      )
    ),

    tags$footer(
      class = "omm12-footer",
      fluidRow(
        column(
          width = 4,
          tags$div(class = "omm12-footer-brand",
            tags$div(
              tags$strong("OMM12 Resource"),
              tags$br(),
              tags$span(class = "omm12-footer-tagline",
                        "A genomic & functional browser for the Oligo-Mouse-Microbiota-12 synthetic gut community.")
            )
          )
        ),
        column(
          width = 3,
          tags$div(class = "omm12-footer-block",
            tags$strong("Developed by"),
            tags$br(),
            tags$a(href = "https://www.mls.ls.tum.de/inmb/startseite/", target = "_blank",
                   class = "omm12-footer-lablogo", title = "Stecher Lab",
                   tags$img(src = "images/branding/stecher_lab_logo.png", alt = "Stecher Lab logo")),
            tags$a(href = "https://www.mls.ls.tum.de/inmb/startseite/", target = "_blank",
                   "Stecher Lab — Chair for Intestinal Microbiome"),
            tags$br(),
            "TUM School of Life Sciences, Technical University of Munich"
          )
        ),
        column(
          width = 5,
          tags$div(class = "omm12-footer-block",
            tags$strong("How to cite"),
            tags$br(),
            tags$span(
              "Niedermeier LS*, Burrichter AG*, Matchado MS*, Stecher B. ",
              tags$em("\"The OMM12 Model After Ten Years: What did we learn about the gut ",
                      "microbiome from a defined microbial community?\""),
              " (*equal contribution)"
            ),
            tags$br(),
            "OMM12 community design: ",
            tags$a(href = "https://doi.org/10.1038/nmicrobiol.2016.215", target = "_blank",
                   "Brugiroux et al., Nature Microbiology 2016"),
            tags$br(),
            tags$span(class = "omm12-footer-tagline",
                      paste0("Data last built: ", format(Sys.Date(), "%B %Y")))
          )
        )
      )
    ),

    actionLink("scroll_to_top", label = HTML("<i class='fas fa-arrow-up'></i> Top"), class = "scroll-to-top-btn"),
    tags$script(HTML('$(function() {
var scrollButton = $(".scroll-to-top-btn");
$(window).scroll(function() {
if ($(this).scrollTop() > 100) {
scrollButton.fadeIn();
} else {
scrollButton.fadeOut();
}
});
scrollButton.on("click", function(e) {
e.preventDefault();
$("html, body").animate({ scrollTop: 0 }, "slow");
return false;
});
});
'))
  )
)

# --- SERVER ---
server <- function(input, output, session) {

  # --- ABOUT PAGE: PUBLICATION STATISTICS ---
  # Aggregated tables (pubstats_*) are loaded once at app startup (global scope).

  output$pub_caption <- renderUI({
    if (!pubstats_available) return(NULL)
    s <- pubstats_summary
    tags$p(style = "color:#666; font-size: 13px;",
      paste0(s$publications, " research articles, reviews and book chapters (", s$year_min, "–", s$year_max,
             ") found by a Dimensions.ai full-text search for the OMM12 name variants (Oligo-MM12, ",
             "OMM12, sDMDMm2, …; exports up to ", format_pub_date(s$dimensions_export_date),
             ") and screened by hand. “Uses OMM12” means use of the community, its strains or ",
             "genomes was confirmed in the full text (", s$uses_omm12, " publications); the others mention or cite ",
             "OMM12. Off-topic hits, preprints, conference abstracts and commentaries are excluded. ",
             "Citations: Scopus, retrieved ", format_pub_date(s$citations_date), ". ",
             s$year_max, " is a partial year. Data sources: Dimensions (Digital Science) and Scopus (Elsevier)."))
  })

  output$pub_stats_summary <- renderUI({
    if (!pubstats_available) {
      return(tags$p(style = "color:#888;", "Publication data not available."))
    }
    s <- pubstats_summary
    n_all <- as.numeric(s$publications); n_used <- as.numeric(s$uses_omm12)
    stat_box <- function(label, value, sub = NULL) {
      column(width = 2,
             div(style = "text-align:center; padding: 10px;",
                 tags$div(style = "font-size: 26px; font-weight: 700; color: var(--navy-800, #1b263b);", value),
                 tags$div(style = "font-size: 11px; color:#666; text-transform: uppercase; letter-spacing: 0.5px;", label),
                 if (!is.null(sub)) tags$div(style = "font-size: 11px; color:#999;", sub)
             )
      )
    }
    fluidRow(
      stat_box("Publications", s$publications),
      stat_box("Use OMM12", s$uses_omm12, paste0(round(100 * n_used / n_all), "% of all")),
      stat_box("Countries", s$countries),
      stat_box("Citations (Scopus)", format(as.numeric(s$citations_scopus_total), big.mark = ",")),
      stat_box("h-index", s$h_index_all, paste0(s$h_index_uses_omm12, " for OMM12 users")),
      stat_box("Open access", paste0(s$open_access_pct, "%"))
    )
  })

  output$pub_year_chart <- renderPlotly({
    req(!is.null(pubstats_years))
    plot_ly(pubstats_years, x = ~year, y = ~publications, color = ~group, colors = PUB_GROUP_COLORS, type = "bar",
            hovertemplate = "%{x}: %{y} publications<extra>%{fullData.name}</extra>") %>%
      layout(barmode = "stack",
             title = list(text = "Publications per year", x = 0, xanchor = "left", font = list(size = 15)),
             xaxis = list(title = "", dtick = 1), yaxis = list(title = "Publications"),
             legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.15, traceorder = "normal"),
             margin = list(t = 50))
  })

  output$pub_citations_chart <- renderPlotly({
    req(!is.null(pubstats_years))
    plot_ly(pubstats_years, x = ~year, y = ~citations_scopus, color = ~group, colors = PUB_GROUP_COLORS, type = "bar",
            hovertemplate = "Published %{x}: %{y:,} citations<extra>%{fullData.name}</extra>") %>%
      layout(barmode = "stack",
             title = list(text = "Citations (Scopus) by publication year", x = 0, xanchor = "left",
                          font = list(size = 15)),
             xaxis = list(title = "", dtick = 1), yaxis = list(title = "Citations"),
             legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.15, traceorder = "normal"),
             margin = list(t = 50))
  })

  output$pub_world_map <- renderPlotly({
    req(!is.null(pubstats_countries))
    counts <- pubstats_countries
    counts$n <- if (identical(input$pub_map_group, "used")) counts$publications_uses_omm12 else counts$publications_all
    counts <- counts[counts$n > 0, ]
    req(nrow(counts) > 0)
    ticks <- c(1, 3, 10, 30, 100, 300)
    ticks <- ticks[ticks <= max(counts$n) * 1.5]
    plot_ly(counts, type = "choropleth", locations = ~iso3, z = ~log10(n),
            text = ~paste0(country, ": ", n, " publication", ifelse(n == 1, "", "s")),
            hoverinfo = "text", zmin = 0, zmax = log10(max(counts$n)),
            colorscale = list(c(0, "#FDDEA0"), c(0.2, "#F4A460"), c(0.4, "#E07B54"),
                              c(0.6, "#B85C8A"), c(0.8, "#7B3F9E"), c(1, "#4B1F7A")),
            marker = list(line = list(color = "#ffffff", width = 0.5)),
            colorbar = list(title = "Publications", tickvals = log10(ticks), ticktext = ticks)) %>%
      layout(title = paste0("Publications by country (", length(unique(counts$iso3)), " countries)"),
             geo = list(showframe = FALSE, showcoastlines = FALSE, showcountries = TRUE,
                        countrycolor = "#ffffff", showland = TRUE, landcolor = "#e9e9ee",
                        projection = list(type = "natural earth"), lataxis = list(range = c(-58, 85))))
  })

  output$pub_journals_chart <- renderPlotly({
    req(!is.null(pubstats_journals))
    order <- unique(pubstats_journals$journal[order(-pubstats_journals$journal_total, pubstats_journals$journal)])
    top <- pubstats_journals %>% mutate(journal = factor(journal, levels = rev(order)))
    plot_ly(top, x = ~publications, y = ~journal, color = ~group, colors = PUB_GROUP_COLORS, type = "bar",
            orientation = "h",
            text = ~paste0(journal, "<br>", publications, " (", group, ") of ", journal_total, " publications<br>",
                           format(journal_citations_scopus, big.mark = ","), " citations (Scopus, all)"),
            hoverinfo = "text", textposition = "none") %>%
      layout(barmode = "stack",
             title = list(text = "Top 10 peer-reviewed journals", x = 0, xanchor = "left", font = list(size = 15)),
             xaxis = list(title = "Publications (research articles and reviews)"),
             yaxis = list(title = ""), margin = list(l = 230, t = 50),
             legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.2, traceorder = "normal"))
  })

  current_bacterium <- reactiveVal(NULL)
  active_content_type <- reactiveVal(NULL)

  
  # --- REACTIVES: STATS/COUNTS ---
  current_bacterium_localization_data <- reactive({
    req(current_bacterium())
    load_localization_csv(current_bacterium()$filename)
  })
  
  current_prokka_gff <- reactive({
    req(current_bacterium())
    load_prokka_gff(current_bacterium()$filename)
  })
  
  current_prokka_tsv <- reactive({
    req(current_bacterium())
    load_prokka_tsv(current_bacterium()$filename)
  })
  
  current_bacterium_cog_data <- reactive({
    req(current_bacterium())
    load_cog_tsv(current_bacterium()$filename)
  })
  
  current_bacterium_mobileog_data <- reactive({
    req(current_bacterium())
    load_mobileog_csv(current_bacterium()$filename)
  })
  
  current_bacterium_kegg_data <- reactive({
    req(current_bacterium())
    load_kegg_tsv(current_bacterium()$filename)
  })
  
    output$gem_content <- renderUI({
    req(active_content_type() == "metabolic_model")
    gem <- current_gem()
    if (is.null(gem) || nrow(gem$reactions) == 0) {
      return(box(width = 12, status = "warning", solidHeader = TRUE,
                 p(icon("circle-info"), " A genome-scale metabolic model is not yet available for this strain.")))
    }
    box(
      title = paste("Metabolic Model for", current_bacterium()$name),
      width = 12, solidHeader = TRUE, status = "info",
      fluidRow(column(12, tags$h5(paste0(
        "Genome-scale model reconstructed with gapseq: ",
        nrow(gem$reactions), " reactions across ",
        length(unique(gem$subsystems$subsystem)), " pathways. ",
        "Select a pathway to see its reactions and the genes that catalyse them.")))),
      fluidRow(
        column(9, p(style = "color:#666; font-size: 12.5px;",
          "Source: Zimmermann J, Burrichter A (2025). Metabolic models of OMM12 community (gapseq v1.4). Zenodo. ",
          tags$a(href = "https://doi.org/10.5281/zenodo.17358311", target = "_blank",
                 "doi:10.5281/zenodo.17358311"))),
        column(3, style = "text-align: right;",
               downloadButton("download_gem_sbml", "Download model (SBML)", class = "btn-sm"))
      ),

      tags$h5("Pathway landscape", style = "margin-top: 20px;"),
      p(style = "color:#666; font-size: 13px;",
        "Number of reactions per pathway. Click a bar to jump to that pathway below."),
      shinycssloaders::withSpinner(plotlyOutput("gem_subsystem_overview", height = "500px"),
                                    type = 6, color = "#1b263b"),

      tags$hr(style = "border-top: 1px solid #ccc;"),

      # --- STEP 4: server-side selectize (empty here; filled by the observer below) ---
      fluidRow(
        column(6, selectizeInput(
          "gem_subsystem", "Metabolic pathway:",
          choices  = gem_subsystem_choices(),
          selected = gem_subsystem_choices()[1],
          options = list(placeholder = "Type to search pathways…", maxOptions = 50)
        )),
        column(6, div(style = "padding-top: 25px;", uiOutput("gem_metacyc_link")))
      ),

      tags$hr(style = "border-top: 1px solid #ccc;"),

      # --- STEP 3: spinner-wrapped table ---
      shinycssloaders::withSpinner(DTOutput("gem_table"), type = 6, color = "#1b263b"),

      tags$hr(style = "border-top: 1px solid #ccc;"),

      tags$h5("Pathway network", style = "margin-top: 10px;"),
      p(style = "color:#666; font-size: 13px;",
        "Metabolites (gold) and reactions (navy) for the selected pathway, ",
        "linked substrate → reaction → product. Drag nodes to untangle, scroll to zoom."),
      shinycssloaders::withSpinner(forceNetworkOutput("gem_network", height = "550px"),
                                    type = 6, color = "#1b263b")
    )
  })
  
  # gapseq keeps the real hyphenated MetaCyc pathway ID as the subsystem name
  # (e.g. "PWY-1121"), so we can link straight out to MetaCyc for the readable
  # pathway description instead of maintaining an offline ID->name table.
  output$gem_metacyc_link <- renderUI({
    req(input$gem_subsystem)
    pwy_id <- input$gem_subsystem
    url <- paste0("https://metacyc.org/pathway?orgid=META&id=", pwy_id)
    tags$a(
      href = url, target = "_blank", rel = "noopener noreferrer",
      class = "btn btn-sm", style = "background-color:#c9a227; color:#0d1b2a; font-weight:600;",
      icon("external-link-alt"),
      paste0(" View \"", pwy_id, "\" on MetaCyc")
    )
  })

  gem_subsystem_choices <- reactive({
    gem <- current_gem(); req(gem)
    counts <- gem$subsystems %>%
      group_by(subsystem) %>% summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n))
    setNames(counts$subsystem, paste0(counts$subsystem, "  (", counts$n, " rxns)"))
  })
  
  gem_filtered <- reactive({
    gem <- current_gem(); req(gem, input$gem_subsystem)
    rids <- gem$subsystems$rid[gem$subsystems$subsystem == input$gem_subsystem]
    gem$reactions %>%
      filter(rid %in% rids) %>%
      transmute(
        `Reaction ID` = sub("^R_", "", rid),
        `Name`        = name,
        `Equation`    = equation,
        `Reversible`  = ifelse(reversible, "Yes", "No"),
        `Genes`       = genes
      )
  })

  # --- Pathway landscape: reaction count per subsystem, click a bar to select it ---
  gem_subsystem_counts <- reactive({
    gem <- current_gem(); req(gem)
    gem$subsystems %>%
      group_by(subsystem) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n))
  })

  output$gem_subsystem_overview <- renderPlotly({
    counts <- gem_subsystem_counts()
    req(counts, nrow(counts) > 0)
    counts <- counts %>% slice_head(n = 25)
    counts$subsystem <- factor(counts$subsystem, levels = rev(counts$subsystem))

    plot_ly(
      data = counts, x = ~n, y = ~subsystem, type = "bar", orientation = "h",
      marker = list(color = "#2c3e57"),
      hovertemplate = "%{y}<br>%{x} reactions<extra></extra>",
      source = "gem_overview_clicks"
    ) %>%
      layout(
        xaxis = list(title = "Reactions"), yaxis = list(title = ""),
        margin = list(l = 260)
      ) %>%
      plotly::event_register("plotly_click")
  })

  observeEvent(event_data("plotly_click", source = "gem_overview_clicks"), {
    ed <- event_data("plotly_click", source = "gem_overview_clicks")
    if (!is.null(ed) && "y" %in% names(ed)) {
      updateSelectizeInput(session, "gem_subsystem", selected = as.character(ed$y))
    }
  })

  # --- Pathway network: metabolite <-> reaction bipartite graph for the selected subsystem ---
  gem_network_data <- reactive({
    gem <- current_gem(); req(gem, input$gem_subsystem)
    rids <- gem$subsystems$rid[gem$subsystems$subsystem == input$gem_subsystem]
    rxns <- gem$reactions %>% filter(rid %in% rids)
    req(nrow(rxns) > 0)

    # Cap pathway size so the graph stays readable rather than an unreadable hairball
    if (nrow(rxns) > 60) rxns <- rxns %>% slice_head(n = 60)

    parse_side <- function(s) {
      s <- trimws(s)
      if (!nzchar(s)) return(character(0))
      terms <- strsplit(s, " + ", fixed = TRUE)[[1]]
      trimws(sub("^[0-9]+(\\.[0-9]+)?\\s+", "", terms))
    }

    edge_list <- lapply(seq_len(nrow(rxns)), function(i) {
      r <- rxns[i, ]
      eq <- r$equation
      if (is.na(eq) || !nzchar(eq)) return(NULL)
      arrow <- if (grepl(" <=> ", eq, fixed = TRUE)) " <=> " else " => "
      sides <- strsplit(eq, arrow, fixed = TRUE)[[1]]
      if (length(sides) != 2) return(NULL)
      reactants <- parse_side(sides[1])
      products  <- parse_side(sides[2])
      out <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
      if (length(reactants) > 0)
        out <- rbind(out, data.frame(from = reactants, to = r$rid, stringsAsFactors = FALSE))
      if (length(products) > 0)
        out <- rbind(out, data.frame(from = r$rid, to = products, stringsAsFactors = FALSE))
      out
    })
    edges <- bind_rows(edge_list)
    req(nrow(edges) > 0)

    rxn_labels <- setNames(
      ifelse(!is.na(rxns$name) & nzchar(rxns$name), rxns$name, rxns$rid),
      rxns$rid
    )
    all_ids <- unique(c(edges$from, edges$to))
    is_rxn <- all_ids %in% rxns$rid

    nodes <- data.frame(
      id    = all_ids,
      label = ifelse(is_rxn, unname(rxn_labels[all_ids]), all_ids),
      group = ifelse(is_rxn, "Reaction", "Metabolite"),
      stringsAsFactors = FALSE
    )

    links <- edges %>%
      mutate(
        source = match(from, nodes$id) - 1,
        target = match(to, nodes$id) - 1,
        value  = 1
      )

    list(nodes = nodes, links = links)
  })

  output$gem_network <- renderForceNetwork({
    gd <- gem_network_data()
    req(gd, nrow(gd$nodes) > 0, nrow(gd$links) > 0)

    forceNetwork(
      Links = gd$links, Nodes = gd$nodes,
      Source = "source", Target = "target", Value = "value",
      NodeID = "label", Group = "group",
      opacity = 0.9, zoom = TRUE, arrows = TRUE, bounded = TRUE,
      fontSize = 13, fontFamily = "sans-serif",
      linkDistance = 80, charge = -150,
      colourScale = htmlwidgets::JS(
        "d3.scaleOrdinal().domain(['Reaction','Metabolite']).range(['#0d1b2a','#c9a227'])"
      )
    )
  })
  
  kegg_filter <- reactive({
    clicked <- selected_kegg_node()
    req(clicked, nzchar(clicked))
    df <- current_bacterium_kegg_data()
    req(df)
    clicked <- trimws(clicked)
    
    keep <- (trimws(df$Main_Pathway) == clicked) |
      (trimws(df$Sub_pathway)  == clicked) |
      (trimws(df$Pathway)      == clicked)
    keep[is.na(keep)] <- FALSE
    df[keep, , drop = FALSE]
  })
  
  current_defense_systems <- reactive({
    req(current_bacterium())
    load_defense_systems(current_bacterium()$filename)
  })
  
  current_defense_genes <- reactive({
    req(current_bacterium())
    load_defense_genes(current_bacterium()$filename)
  })
  current_feature_table <- reactive({
    req(current_bacterium())
    load_feature_table(current_bacterium()$filename)
  })
  
  # counts by system type, for the donut
  # counts by system type, filtered by the chosen activity
  defense_type_counts <- reactive({
    raw <- current_defense_systems()
    req(input$defense_activity)            # the radio selection
    if (is.null(raw) || nrow(raw) == 0) {
      return(data.frame(Type = "No data", Count = 0, stringsAsFactors = FALSE))
    }
    sel_act <- input$defense_activity
    filt <- raw %>% filter(!is.na(activity), trimws(activity) == sel_act)
    
    if (nrow(filt) == 0) {
      return(data.frame(Type = "No data", Count = 0, stringsAsFactors = FALSE))
    }
    filt %>%
      filter(!is.na(type), type != "") %>%
      group_by(Type = type) %>%
      summarise(Count = n(), .groups = "drop") %>%
      arrange(desc(Count))
  })
  
  selected_defense_type <- reactiveVal(NULL)
  
  filtered_defense_data_reactive <- reactive({
    sys <- current_defense_systems()
    req(input$defense_activity)
    if (is.null(sys) || nrow(sys) == 0) {
      return(data.frame(Message = "No defense system data loaded.", stringsAsFactors = FALSE))
    }
    sel_act <- input$defense_activity
    req(selected_defense_type())
    sel_type <- selected_defense_type()
      genes <- current_defense_genes()
    ft    <- current_feature_table()
    
    if (!is.null(genes) && nrow(genes) > 0 && "sys_id" %in% colnames(genes)) {
      g <- genes %>%
        filter(type == sel_type, trimws(activity) == sel_act) %>%
        mutate(hit_id = trimws(hit_id))
      
      if (!is.null(ft)) {
        g <- g %>% left_join(ft, by = c("hit_id" = "product_accession"))
      } else {
        g <- g %>% mutate(ft_locus_tag = NA, ft_gene = NA,
                          ft_start = NA, ft_end = NA, ft_strand = NA)
      }
      
      out <- g %>%
        transmute(
          `System`     = type,
          `Subtype`    = subtype,
          `Activity`   = activity,
          `Locus Tag`  = ft_locus_tag,
          `Gene`       = ifelse(!is.na(ft_gene) & ft_gene != "", ft_gene, gene_name),
          `Protein ID` = hit_id,
          `Start`      = ft_start,
          `End`        = ft_end,
          `Strand`     = ft_strand,
          `i-Evalue`   = hit_i_eval,
          `Score`      = hit_score,
          `System ID`  = sys_id,
          `Gene Name`   = ft_gene_name
        ) %>%
        distinct()
      
      if (nrow(out) == 0) {
        return(data.frame(Message = "No genes for this selection.", stringsAsFactors = FALSE))
      }
      out
      
    } else {
      # systems-only fallback (unchanged from before)
      sys %>%
        filter(type == sel_type, trimws(activity) == sel_act) %>%
        transmute(
          `System`          = type,
          `Subtype`         = subtype,
          `Activity`        = activity,
          `Genes in system` = genes_count,
          `Profiles`        = name_of_profiles_in_sys,
          `System ID`       = sys_id
        ) %>%
        distinct()
    }
    
    if (nrow(out) == 0) {
      return(data.frame(Message = "No genes for this selection.", stringsAsFactors = FALSE))
    }
    out
  })
  
  # NEW: Add similarity data reactive
  current_bacterium_similarity_data <- reactive({
    req(current_bacterium())
    load_similarity_data(current_bacterium()$filename)
  })
  
  
  
  cog_legend <- reactive({
    legend_file <- "data/cog_legend.tsv"
    if (!file.exists(legend_file)) { return(NULL) }
    legend_data <- read_tsv(legend_file, show_col_types = FALSE)
    colnames(legend_data) <- c("Group", "Letter", "Description")
    legend_data
  })
  
  cog_colors <- reactive({
    legend <- cog_legend()
    req(legend)
    group_colors <- c(
      "INFORMATION STORAGE AND PROCESSING" = "#4C72B0",
      "CELLULAR PROCESSES AND SIGNALING" 	= "#DD8452",
      "METABOLISM" 	= "#55A868",
      "POORLY CHARACTERIZED" 	= "#8172B2",
      "UNKNOWN"="grey"
    )
    colors <- setNames(group_colors[legend$Group], legend$Letter)
    colors
  })
  
  cog_category_counts <- reactive({
    data <- current_bacterium_cog_data()
    req(data)

    expanded <- data %>%
      mutate(COG_category = toupper(COG_category)) %>%
      # Merge "S" (Function unknown) and "Unknown"/blank/NA (Not categorised)
      # into a single "S" bucket so they show as one bar in the plot.
      mutate(COG_category = ifelse(is.na(COG_category) | COG_category %in% c("", "S", "UNKNOWN"),
                                    "S", COG_category)) %>%
      rowwise()
    #mutate(COG_category_split = strsplit(COG_category, "")[[1]]) %>%
    #tidyr::unnest(COG_category_split)

    counts <- expanded %>%
      group_by(COG_category) %>%
      summarise(Count = n(), .groups = "drop") %>%
      arrange(desc(Count))

    legend <- cog_legend()
    req(legend)
    counts <- left_join(counts, legend, by = c("COG_category" = "Letter"))
    counts$Description[counts$COG_category == "S"] <- "Function unknown / not categorised"
    return(counts)
  })
  
  # 2D projection (PCA) of each protein's 6-class DeepLocPro probability vector,
  # colored by predicted localization. NOTE: this is built from the model's
  # output probabilities, not its internal sequence embeddings (unlike the
  # UMAP in the DeepLoc2 paper) -- proteins separate mostly because the
  # predicted class *is* the argmax of these same numbers, so treat this as an
  # illustrative overview rather than an independent validation of clustering.
  localization_pca_data <- reactive({
    raw_data <- current_bacterium_localization_data()
    req(raw_data, nrow(raw_data) > 1)

    prob_cols <- c("Cell_wall_surface", "Extracellular", "Cytoplasmic",
                    "Cytoplasmic_Membrane", "Outer_Membrane", "Periplasmic")
    present_cols <- intersect(prob_cols, colnames(raw_data))
    req(length(present_cols) >= 2)

    mat <- raw_data %>%
      select(all_of(present_cols)) %>%
      mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
      mutate(across(everything(), ~ ifelse(is.na(.x), 0, .x)))

    keep <- apply(mat, 2, function(col) stats::sd(col) > 0)
    mat <- mat[, keep, drop = FALSE]
    req(ncol(mat) >= 2)

    pca <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
    var_explained <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

    data.frame(
      Protein_name = raw_data$Protein_name,
      Localization = raw_data$Localization,
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      var1 = var_explained[1],
      var2 = var_explained[2],
      stringsAsFactors = FALSE
    )
  })

  output$localization_pca_plot <- renderPlotly({
    df <- localization_pca_data()
    req(df, nrow(df) > 0)

    loc_colors <- c(
      "Cytoplasmic"          = "#0d1b2a",
      "Cytoplasmic Membrane" = "#c9a227",
      "Periplasmic"          = "#2c8fa1",
      "Cell wall & surface"  = "#8b5e3c",
      "Outer Membrane"       = "#7b4b94",
      "Extracellular"        = "#c94f4f"
    )

    plot_ly(
      data = df, x = ~PC1, y = ~PC2, type = "scatter", mode = "markers",
      color = ~Localization, colors = loc_colors,
      text = ~Protein_name,
      hovertemplate = "%{text}<br>%{fullData.name}<extra></extra>",
      marker = list(size = 6, opacity = 0.75)
    ) %>%
      layout(
        xaxis = list(title = paste0("PC1 (", df$var1[1], "%)")),
        yaxis = list(title = paste0("PC2 (", df$var2[1], "%)")),
        legend = list(title = list(text = "Predicted localization"))
      )
  })

  localization_word_counts <- reactive({
    raw_data <- current_bacterium_localization_data()
    if (is.null(raw_data) || nrow(raw_data) == 0) {
      return(data.frame(word = "No Data", freq = 1, stringsAsFactors = FALSE))
    }
    counts <- raw_data %>%
      filter(!is.na(Localization), Localization != "") %>%
      group_by(word = Localization) %>%
      summarise(freq = n(), .groups = "drop") %>%
      arrange(desc(freq))
    return(counts)
  })

  output$localization_cell_diagram <- renderUI({
    counts_df <- localization_word_counts()
    req(counts_df, "word" %in% colnames(counts_df))
    counts_named <- setNames(as.list(counts_df$freq), counts_df$word)
    HTML(build_localization_cell_diagram(counts_named))
  })

  # --- Localization x COG functional group cross-tab ---
  localization_cog_crosstab <- reactive({
    loc <- current_bacterium_localization_data()
    cog <- current_bacterium_cog_data()
    legend <- cog_legend()
    req(loc, cog, legend, nrow(loc) > 0, nrow(cog) > 0)

    cog_grouped <- cog %>%
      mutate(COG_category = toupper(COG_category)) %>%
      left_join(legend, by = c("COG_category" = "Letter")) %>%
      mutate(Group = ifelse(is.na(Group) | Group == "", "UNKNOWN", Group))

    joined <- loc %>%
      filter(!is.na(Localization), Localization != "") %>%
      inner_join(cog_grouped %>% select(locus_tag, Group), by = c("Protein_name" = "locus_tag"))

    req(nrow(joined) > 0)

    joined %>%
      count(Localization, Group, name = "n") %>%
      tidyr::complete(Localization, Group, fill = list(n = 0))
  })

  output$localization_cog_heatmap <- renderPlotly({
    df <- localization_cog_crosstab()
    req(df, nrow(df) > 0)

    plot_ly(
      data = df, x = ~Group, y = ~Localization, z = ~n, type = "heatmap",
      colors = colorRamp(c("#f4f6f9", "#c9a227", "#0d1b2a")),
      hovertemplate = "%{y} × %{x}<br>%{z} proteins<extra></extra>",
      showscale = TRUE,
      source = "localization_cog_heatmap_clicks"
    ) %>%
      layout(
        xaxis = list(title = "", tickangle = -30),
        yaxis = list(title = ""),
        margin = list(b = 150, l = 160)
      ) %>%
      plotly::event_register("plotly_click")
  })

  annotation_stats <- reactive({
    gff_data <- current_prokka_gff()
    if (is.null(gff_data)) {
      return(data.frame(type = "No data", Count = 0, stringsAsFactors = FALSE))
    }
    req(gff_data)
    get_annotation_stats(gff_data)
  })
  
  mobileog_category_counts <- reactive({
    raw_data <- current_bacterium_mobileog_data()
    if (is.null(raw_data) || nrow(raw_data) == 0) {
      return(data.frame(Category = "No data", Count = 0, stringsAsFactors = FALSE))
    }
    counts <- raw_data %>%
      filter(!is.na(`Major mobileOG Category`), `Major mobileOG Category` != "") %>%
      group_by(Category = `Major mobileOG Category`) %>%
      summarise(Count = n(), .groups = "drop") %>%
      arrange(desc(Count))
    return(counts)
  })
  
  # --- REACTIVES: SELECTION & FILTERING ---
  
  selected_localization_type <- reactiveVal(NULL)
  selected_annotation_type <- reactiveVal(NULL)
  selected_mobileog_category <- reactiveVal(NULL)
  selected_cog_category <- reactiveVal(NULL)
  selected_clade <- reactiveVal(NULL)
  # Set by the Cross-Genome Search results table so the Genome Map's gene
  # search box can be pre-filled with the gene the user clicked through on.
  pending_gene_search <- reactiveVal(NULL)
  
  # COG Filtered Data
  filtered_cog_data_reactive <- reactive({
    raw_data <- current_bacterium_cog_data()
    if (is.null(raw_data) || nrow(raw_data) == 0) {
      return(data.frame(Message = "No COG data loaded.", stringsAsFactors = FALSE))
    }
    req(selected_cog_category())
    selected_letter <- toupper(selected_cog_category())

    if (selected_letter == "S") {
      # "S" bar now represents the merged Function-unknown / Not-categorised bucket
      filtered <- raw_data %>%
        filter(is.na(COG_category) | toupper(COG_category) %in% c("", "S", "UNKNOWN"))
    } else {
      filtered <- raw_data %>%
        filter(str_detect(toupper(COG_category), selected_letter))
    }

    # COGclassifier's locus_tag is the raw Prokka ID (e.g. "GFAIFJLA_00001"),
    # not useful on its own. Join against the Prokka annotation table to pull
    # in the real gene symbol and product name, same as the Localization tab.
    prokka <- current_prokka_tsv()
    if (!is.null(prokka) && nrow(prokka) > 0 && "locus_tag" %in% colnames(prokka)) {
      filtered <- filtered %>%
        left_join(prokka %>% select(locus_tag, gene, product), by = "locus_tag")
    } else {
      filtered <- filtered %>% mutate(gene = NA_character_, product = NA_character_)
    }

    filtered <- filtered %>%
      transmute(
        `Locus Tag`     = locus_tag,
        `Gene`          = ifelse(is.na(gene) | gene == "", "-", gene),
        `Protein`       = ifelse(is.na(product) | product == "", "Unknown / hypothetical protein", product),
        `COG Category`  = COG_category,
        `Description`   = Description,
        `PFAMs`         = PFAMs
      ) %>%
      distinct()

    return(filtered)
  })
  
  # Localization Filtered Data
  filtered_localization_data_reactive <- reactive({
    raw_data <- current_bacterium_localization_data()
    if (is.null(raw_data) || nrow(raw_data) == 0) { return(data.frame(Message = "No Localization data loaded.", stringsAsFactors = FALSE)) }
    req(selected_localization_type())
    selected_loc <- tolower(selected_localization_type())

    # Protein_name in the DeepLoc output is just the Prokka locus tag (e.g.
    # "GFAIFJLA_00001"), not a real protein name. Join against the Prokka
    # annotation table to pull in the actual gene symbol and product
    # description, so the table shows something meaningful instead of an ID.
    prokka <- current_prokka_tsv()
    annotated <- raw_data
    if (!is.null(prokka) && nrow(prokka) > 0 && "locus_tag" %in% colnames(prokka)) {
      annotated <- raw_data %>%
        left_join(
          prokka %>% select(locus_tag, gene, product),
          by = c("Protein_name" = "locus_tag")
        )
    } else {
      annotated <- raw_data %>% mutate(gene = NA_character_, product = NA_character_)
    }

    filtered <- annotated %>%
      filter(tolower(Localization) == selected_loc) %>%
      transmute(
        `Locus Tag`         = Protein_name,
        `Gene`              = ifelse(is.na(gene) | gene == "", "-", gene),
        `Protein`           = ifelse(is.na(product) | product == "", "Unknown / hypothetical protein", product),
        `Localization`      = Localization,
        `Cell wall/surface` = Cell_wall_surface,
        `Extracellular`     = Extracellular,
        `Cytoplasmic`       = Cytoplasmic,
        `Cytoplasmic membrane` = Cytoplasmic_Membrane
      ) %>%
      distinct()
    return(filtered)
  })
  
  # Annotation Filtered Data
  filtered_annotation_data_reactive <- reactive({
    gff_data <- current_prokka_gff()
    if (is.null(gff_data) || nrow(gff_data) == 0) {
      return(data.frame(Message = "No GFF data available for detailed filtering.", stringsAsFactors = FALSE))
    }
    req(selected_annotation_type())
    
    filtered <- gff_data %>%
      filter(type == selected_annotation_type()) %>%
      select(Contig = seqid, Feature_Type = type, Locus_Tag = locus_tag, Product = product, Start = start, End = end, Strand = strand, COG= COG_category)
    return(filtered)
  })
  
  # MobileOG Filtered Data
  filtered_mobileog_data_reactive <- reactive({
    raw_data <- current_bacterium_mobileog_data()
    if (is.null(raw_data) || nrow(raw_data) == 0) { return(data.frame(Message = "No MobileOG data loaded.", stringsAsFactors = FALSE)) }
    req(selected_mobileog_category())

    # The mobileOG hit's own "Gene Name" is the reference database protein's
    # symbol, not this bacterium's own annotation. Where the hit's query ID
    # happens to be a Prokka locus tag (true for pipelines run directly on
    # Prokka output), join back to Prokka to get this organism's real gene
    # symbol and product/description instead of an opaque ID.
    prokka <- current_prokka_tsv()
    join_col <- NA_character_
    if (!is.null(prokka) && nrow(prokka) > 0 && "locus_tag" %in% colnames(prokka)) {
      locus_tags <- prokka$locus_tag
      if ("Query Title" %in% names(raw_data) &&
          mean(raw_data$`Query Title` %in% locus_tags) > 0.5) {
        join_col <- "Query Title"
      } else if ("Contig/ORF Name" %in% names(raw_data) &&
                 mean(raw_data$`Contig/ORF Name` %in% locus_tags) > 0.5) {
        join_col <- "Contig/ORF Name"
      }
    }

    if (!is.na(join_col)) {
      by_vec <- setNames("locus_tag", join_col)
      raw_data <- raw_data %>%
        left_join(prokka %>% select(locus_tag, prokka_gene = gene, prokka_product = product),
                   by = by_vec)
    } else {
      raw_data <- raw_data %>% mutate(prokka_gene = NA_character_, prokka_product = NA_character_)
    }

    filtered <- raw_data %>%
      filter(`Major mobileOG Category` == selected_mobileog_category()) %>%
      mutate(
        `Gene`    = ifelse(!is.na(prokka_gene) & prokka_gene != "", prokka_gene,
                     ifelse(!is.na(`Gene Name`) & `Gene Name` != "" & `Gene Name` != "NA:Keyword",
                            `Gene Name`, "-")),
        `Protein` = ifelse(!is.na(prokka_product) & prokka_product != "", prokka_product,
                            "Unknown / hypothetical protein")
      )

    keep_cols <- c("Gene", "Protein", "Best Hit Accession ID", "Minor mobileOG Category",
                   "Evidence Type", "Pident", "Subject Sequence Length", "e-value")
    keep_cols <- intersect(keep_cols, names(filtered))

    filtered <- filtered %>% select(all_of(keep_cols)) %>% distinct()
    return(filtered)
  })

  # --- CROSS-GENOME SEARCH ---

  cross_search_hits <- reactive({
    req(input$cross_search_query, nzchar(trimws(input$cross_search_query)))
    df <- all_genes_index
    req(nrow(df) > 0)

    q <- trimws(input$cross_search_query)
    search_type <- input$cross_search_type
    if (is.null(search_type)) search_type <- "any"

    match_target <- switch(search_type,
      "gene"    = df$Gene,
      "product" = df$Product,
      paste(df$Gene, df$Product)
    )

    df[contains_ci(match_target, q), , drop = FALSE] %>%
      arrange(Bacterium, Gene)
  })

  output$cross_search_summary <- renderUI({
    req(input$cross_search_query, nzchar(trimws(input$cross_search_query)))
    hits <- cross_search_hits()
    n_bact <- length(unique(hits$Bacterium))
    tags$p(
      style = "color:#666; font-size: 13px; margin-top: 10px; margin-bottom: 0;",
      sprintf("Found %d matching gene(s) across %d of the 12 community members.", nrow(hits), n_bact)
    )
  })

  output$cross_search_results <- renderDT({
    if (is.null(input$cross_search_query) || !nzchar(trimws(input$cross_search_query))) {
      return(datatable(data.frame(Message = "Enter a gene symbol or keyword above to search."),
                        rownames = FALSE, options = list(dom = 't')))
    }
    hits <- cross_search_hits()
    if (nrow(hits) == 0) {
      return(datatable(data.frame(Message = "No matches found."),
                        rownames = FALSE, options = list(dom = 't')))
    }
    display <- hits %>%
      transmute(
        Bacterium, Gene, Product,
        `COG Category` = ifelse(is.na(COG_Category), "-", COG_Category),
        `Locus Tag` = Locus_Tag, Contig, Start, End, Strand
      )
    datatable(display, options = list(pageLength = 15, scrollX = TRUE),
              rownames = FALSE, selection = "single")
  })

  # Clicking a result row jumps straight to that bacterium's Genome Map with
  # the gene pre-filled in the search box (reuses the region-jump logic that
  # already exists for manual gene search there).
  observeEvent(input$cross_search_results_rows_selected, {
    sel <- input$cross_search_results_rows_selected
    req(sel)
    hits <- cross_search_hits()
    req(nrow(hits) >= sel)
    row <- hits[sel, ]
    idx <- match(row$Filename, bacteria_filenames)
    req(!is.na(idx))

    current_bacterium(list(name = bacteria_names[idx], filename = bacteria_filenames[idx], index = idx))
    selected_localization_type(NULL)
    selected_annotation_type(NULL)
    selected_mobileog_category(NULL)
    selected_cog_category(NULL)
    selected_clade(NULL)
    pending_gene_search(row$Gene)

    updateTabItems(session, "sidebar_menu", selected = "details_page")
    active_content_type("genome")
  })

  # --- Sequence Search ---
  seq_search_results_data <- reactiveVal(data.frame())
  seq_search_error <- reactiveVal(NULL)

  observeEvent(input$seq_search_run, {
    seq_search_error(NULL)
    raw <- input$seq_search_query

    if (is.null(raw) || !nzchar(trimws(raw))) {
      seq_search_error("Paste a protein sequence first.")
      seq_search_results_data(data.frame())
      return(invisible(NULL))
    }

    # Strip a FASTA header line if present, and any whitespace/newlines.
    lines <- strsplit(raw, "\n")[[1]]
    lines <- lines[!grepl("^>", lines)]
    clean_seq <- toupper(gsub("[^A-Za-z]", "", paste(lines, collapse = "")))

    if (nchar(clean_seq) < 8) {
      seq_search_error("Sequence is too short (need at least 8 residues) or contains no valid amino-acid letters.")
      seq_search_results_data(data.frame())
      return(invisible(NULL))
    }
    if (nchar(clean_seq) > 20000) {
      seq_search_error("Sequence is too long (max 20,000 residues) for this lightweight search.")
      seq_search_results_data(data.frame())
      return(invisible(NULL))
    }

    targets <- if (identical(input$seq_search_target, "all")) bacteria_filenames else input$seq_search_target
    top_n <- input$seq_search_top_n
    if (is.null(top_n) || is.na(top_n) || top_n < 1) top_n <- 25

    hits <- tryCatch(
      search_sequence_r(clean_seq, targets, top_n = top_n),
      error = function(e) {
        seq_search_error(paste("Search failed:", e$message))
        data.frame()
      }
    )

    if (nrow(hits) == 0 && is.null(seq_search_error())) {
      seq_search_error("No hits passed the similarity thresholds (score ≥ 0.10, ≥ 3 shared k-mers). Try a different sequence or search all 12 bacteria.")
    }

    seq_search_results_data(hits)
  })

  output$seq_search_summary <- renderUI({
    err <- seq_search_error()
    if (!is.null(err)) {
      return(p(style = "color:#a33; font-size: 13px;", err))
    }
    hits <- seq_search_results_data()
    if (nrow(hits) == 0) {
      return(p(style = "color:#666; font-size: 13px;",
                "Paste a sequence above and click Search."))
    }
    p(style = "color:#666; font-size: 13px;",
      paste0("Found ", nrow(hits), " hit", if (nrow(hits) != 1) "s" else "",
             ", ranked by containment score."))
  })

  output$seq_search_results <- renderDT({
    hits <- seq_search_results_data()
    if (nrow(hits) == 0) return(datatable(data.frame(), options = list(dom = "t")))
    datatable(
      hits,
      rownames = FALSE,
      selection = "none",
      options = list(pageLength = 10, scrollX = TRUE)
    ) %>%
      formatStyle("Score", background = styleColorBar(range(0, 1), "#e8f0fe"))
  })

  # --- PRIMER DESIGN ---
  observeEvent(input$primer_bact, {
    updateSelectizeInput(session, "primer_gene", choices = get_gene_choices_for(input$primer_bact),
                          server = TRUE, selected = "")
  }, ignoreInit = FALSE)

  primer_target_seq <- reactive({
    if (identical(input$primer_source, "custom")) {
      raw <- input$primer_custom_seq
      if (is.null(raw) || !nzchar(trimws(raw))) return(NULL)
      lines <- strsplit(raw, "\n")[[1]]
      lines <- lines[!grepl("^>", lines)]
      toupper(gsub("[^ACGTacgt]", "", paste(lines, collapse = "")))
    } else {
      req(input$primer_bact, input$primer_gene, nzchar(input$primer_gene))
      row <- all_genes_index %>%
        filter(Filename == input$primer_bact, Gene == input$primer_gene)
      if (nrow(row) == 0) return(NULL)
      get_gene_nucleotide_seq(input$primer_bact, row$Locus_Tag[1])
    }
  })

  output$primer_gene_length_note <- renderUI({
    seq <- tryCatch(primer_target_seq(), error = function(e) NULL)
    if (is.null(seq)) {
      return(p(style = "color:#a33; font-size: 12px;", "Could not load a nucleotide sequence for this gene."))
    }
    p(style = "color:#666; font-size: 12px;", paste0("Gene length: ", nchar(seq), " bp."))
  })

  primer_design_result <- eventReactive(input$primer_design_run, {
    seq <- tryCatch(primer_target_seq(), error = function(e) NULL)
    if (is.null(seq) || nzchar(trimws(seq)) == FALSE) {
      return(list(error = "No target sequence available — pick a gene or paste a sequence.",
                  pairs = data.frame(), seq_length = 0))
    }
    len_min <- input$primer_len_min; len_max <- input$primer_len_max
    if (is.null(len_min) || is.null(len_max) || len_min > len_max) { len_min <- 18; len_max <- 25 }
    tryCatch(
      design_primers(
        seq,
        product_size_range = c(input$primer_product_min %||% 100, input$primer_product_max %||% 800),
        n_pairs = input$primer_n_pairs %||% 5,
        len_range = len_min:len_max,
        gc_range = c(input$primer_gc_min %||% 40, input$primer_gc_max %||% 60),
        tm_target = input$primer_tm_target %||% 60
      ),
      error = function(e) list(error = paste("Primer design failed:", e$message),
                                pairs = data.frame(), seq_length = nchar(seq))
    )
  })

  output$primer_design_summary <- renderUI({
    req(primer_design_result())
    res <- primer_design_result()
    if (!is.null(res$error)) {
      return(p(style = "color:#a33; font-size: 13px;", res$error))
    }
    p(style = "color:#666; font-size: 13px;",
      paste0("Target sequence: ", res$seq_length, " bp. Showing ", nrow(res$pairs),
             " candidate primer pair", if (nrow(res$pairs) != 1) "s" else "", ", best first."))
  })

  output$primer_design_results <- renderDT({
    req(primer_design_result())
    res <- primer_design_result()
    if (is.null(res$pairs) || nrow(res$pairs) == 0) return(datatable(data.frame(), options = list(dom = "t")))
    datatable(res$pairs, rownames = FALSE, selection = "none",
              options = list(pageLength = 10, scrollX = TRUE)) %>%
      formatStyle("Hairpin risk", backgroundColor = styleEqual(c(TRUE, FALSE), c("#fbeaea", "white")))
  })

  # --- OMM12 RELATEDNESS TREE (all 12 members) ---
  output$omm12_relatedness_tree <- renderGirafe({
    # Prefer the real 16S NJ+bootstrap phylogeny; fall back to the k-mer/
    # UPGMA sketch only if the tree file is somehow missing.
    p <- get_omm12_16s_tree_plot()
    if (is.null(p)) p <- get_omm12_relatedness_plot()
    req(p)
    girafe(
      ggobj = p,
      width_svg = 8,
      height_svg = 6,
      options = list(
        opts_zoom(min = 0.5, max = 6),
        opts_hover(css = "cursor:pointer;fill:#c9a227;stroke:#c9a227;stroke-width:2px;", reactive = TRUE),
        opts_hover_inv(css = "opacity:0.35;"),
        opts_tooltip(
          css = "background-color:#1f2937;color:#f9fafb;padding:10px;border-radius:6px;font-size:13px;",
          opacity = 0.95, use_fill = FALSE, use_stroke = FALSE
        ),
        opts_toolbar(position = "topright", saveaspng = TRUE, pngname = "omm12_relatedness_tree"),
        opts_sizing(rescale = TRUE, width = 1)
      )
    )
  })

  output$omm12_relatedness_heatmap <- renderPlotly({
    req(omm12_relatedness)
    sim <- omm12_relatedness$similarity
    labels <- bacteria_names[match(rownames(sim), bacteria_filenames)]
    plot_ly(
      x = labels, y = labels, z = sim, type = "heatmap",
      colors = colorRamp(c("#f4f6f9", "#c9a227", "#0d1b2a")),
      hovertemplate = "%{y} vs %{x}<br>similarity: %{z:.3f}<extra></extra>"
    ) %>%
      layout(
        xaxis = list(title = "", tickangle = -45),
        yaxis = list(title = "", automargin = TRUE),
        margin = list(l = 10, b = 140)
      )
  })

  # --- DEFENSE SYSTEMS OVERVIEW (community-wide) ---
  output$defense_overview_missing_note <- renderUI({
    if (length(defense_bacteria_missing) == 0) return(NULL)
    p(style = "color:#a33; font-size: 12.5px; margin-top: 8px;",
      icon("triangle-exclamation"), " No DefenseFinder results yet for: ",
      paste(defense_bacteria_missing, collapse = ", "),
      ". These bacteria are excluded from the charts and counts below rather than shown as zero.")
  })

  output$defense_overview_summary_stats <- renderUI({
    d <- all_defense_systems_data
    if (nrow(d) == 0) {
      return(p(style = "color:#666;", "No defense system data available."))
    }
    n_total <- nrow(d)
    n_bact <- length(unique(d$Bacterium))
    n_types <- length(unique(d$type))
    n_antidef <- sum(d$activity == "Antidefense", na.rm = TRUE)
    fluidRow(
      column(width = 3, tags$div(class = "small-box-like",
             style = "text-align:center; padding: 10px;",
             tags$h3(n_total, style = "margin:0; color:#1b263b;"),
             tags$p("Defense systems found", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding: 10px;",
             tags$h3(n_bact, style = "margin:0; color:#1b263b;"),
             tags$p("Bacteria with data (of 12)", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding: 10px;",
             tags$h3(n_types, style = "margin:0; color:#1b263b;"),
             tags$p("Distinct system types", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding: 10px;",
             tags$h3(n_antidef, style = "margin:0; color:#c9a227;"),
             tags$p("Antidefense (phage-evasion) genes", style = "margin:0; color:#666; font-size:12px;")))
    )
  })

  output$defense_overview_heatmap <- renderPlotly({
    d <- all_defense_systems_data
    if (nrow(d) == 0) {
      return(plotly_empty() %>%
               layout(annotations = list(x = 0.5, y = 0.5, text = "No defense system data available.",
                                          xref = "paper", yref = "paper", showarrow = FALSE)))
    }
    mat <- d %>%
      count(Bacterium, type, name = "n") %>%
      tidyr::pivot_wider(names_from = type, values_from = n, values_fill = 0)

    bact_order <- mat$Bacterium
    type_cols <- setdiff(colnames(mat), "Bacterium")
    z <- as.matrix(mat[, type_cols, drop = FALSE])
    rownames(z) <- bact_order

    plot_ly(
      x = type_cols, y = bact_order, z = z, type = "heatmap",
      colors = colorRamp(c("#f4f6f9", "#c9a227", "#0d1b2a")),
      hovertemplate = "%{y} | %{x}: %{z}<extra></extra>"
    ) %>%
      layout(
        xaxis = list(title = "", tickangle = -45),
        yaxis = list(title = "", automargin = TRUE),
        margin = list(l = 10, b = 120)
      )
  })

  output$defense_overview_bar <- renderPlotly({
    d <- all_defense_systems_data
    if (nrow(d) == 0) {
      return(plotly_empty() %>%
               layout(annotations = list(x = 0.5, y = 0.5, text = "No defense system data available.",
                                          xref = "paper", yref = "paper", showarrow = FALSE)))
    }
    counts <- d %>% count(Bacterium, activity, name = "n")
    plot_ly(counts, x = ~n, y = ~Bacterium, color = ~activity, type = "bar", orientation = "h",
            colors = c("Defense" = "#1b263b", "Antidefense" = "#c9a227")) %>%
      layout(barmode = "stack",
             xaxis = list(title = "Number of systems"),
             yaxis = list(title = "", automargin = TRUE),
             legend = list(orientation = "h", x = 0, y = -0.15))
  })

  output$defense_overview_table <- renderDT({
    d <- all_defense_systems_data
    if (nrow(d) == 0) return(datatable(data.frame(), options = list(dom = "t")))
    if (!identical(input$defense_overview_activity_filter, "all")) {
      d <- d %>% filter(activity == input$defense_overview_activity_filter)
    }
    if (!identical(input$defense_overview_bact_filter, "all")) {
      d <- d %>% filter(Bacterium_filename == input$defense_overview_bact_filter)
    }
    d <- d %>%
      select(Bacterium, Type = type, Subtype = subtype, Activity = activity,
             `Genes in system` = genes_count, Profiles = name_of_profiles_in_sys) %>%
      arrange(Bacterium, Type)
    datatable(d, rownames = FALSE, selection = "none",
              options = list(pageLength = 15, scrollX = TRUE))
  })

  # --- CRISPR ARRAYS & CAS SYSTEMS ---
  output$crispr_missing_note <- renderUI({
    tagList(
      if (length(crispr_bacteria_zero_arrays) > 0) {
        p(style = "color:#1b7a3d; font-size: 12.5px;",
          icon("circle-check"), " CRISPRCasFinder ran clean (no arrays found) for: ",
          paste(crispr_bacteria_zero_arrays, collapse = ", "), ".")
      } else NULL,
      if (length(crispr_bacteria_missing) > 0) {
        p(style = "color:#a33; font-size: 12.5px;",
          icon("triangle-exclamation"), " No CRISPRCasFinder run yet for: ",
          paste(crispr_bacteria_missing, collapse = ", "), ".")
      } else NULL
    )
  })

  output$crispr_summary_stats <- renderUI({
    d <- all_crispr_arrays_data
    if (nrow(d) == 0) {
      return(p(style = "color:#666;", "No CRISPRCasFinder data loaded yet."))
    }
    n_arrays <- nrow(d)
    n_bact <- length(unique(d$Bacterium))
    n_spacers <- sum(as.numeric(d$Spacers_Nb), na.rm = TRUE)
    n_cas_clusters <- if (nrow(all_cas_systems_data) > 0) {
      all_cas_systems_data %>% distinct(Bacterium_filename, Cluster_Start, Cluster_End) %>% nrow()
    } else 0
    fluidRow(
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(n_arrays, style = "margin:0; color:#1b263b;"),
             tags$p("CRISPR arrays found", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(n_spacers, style = "margin:0; color:#c9a227;"),
             tags$p("Total spacers", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(n_cas_clusters, style = "margin:0; color:#5b7c99;"),
             tags$p("Cas gene clusters found", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(paste0(n_bact, " / 12"), style = "margin:0; color:#1b263b;"),
             tags$p("Bacteria with an array", style = "margin:0; color:#666; font-size:12px;")))
    )
  })

  output$crispr_arrays_table <- renderDT({
    d <- all_crispr_arrays_data
    if (nrow(d) == 0) return(datatable(data.frame(), options = list(dom = "t")))
    if (!identical(input$crispr_bact_filter, "all")) {
      d <- d %>% filter(Bacterium_filename == input$crispr_bact_filter)
    }
    d <- d %>%
      select(Bacterium, Contig = Sequence, Start, End, `Length (bp)` = Length,
             Orientation, `Repeat consensus` = DR_Consensus, `Spacers` = Spacers_Nb,
             `Evidence level (1-4)` = Evidence_Level,
             `Repeat conservation %` = Conservation_DRs_pct) %>%
      arrange(Bacterium, Start)
    datatable(d, rownames = FALSE, selection = "none",
              options = list(pageLength = 15, scrollX = TRUE))
  })

  output$cas_systems_table <- renderDT({
    d <- all_cas_systems_data
    if (nrow(d) == 0) return(datatable(data.frame(), options = list(dom = "t")))
    if (!identical(input$cas_bact_filter, "all")) {
      d <- d %>% filter(Bacterium_filename == input$cas_bact_filter)
    }
    d <- d %>%
      select(Bacterium, Contig = Sequence, `Cas type` = Cas_Cluster_Type,
             `Cluster start` = Cluster_Start, `Cluster end` = Cluster_End,
             `Gene` = Gene_Subtype, `Gene start` = Gene_Start, `Gene end` = Gene_End,
             Strand = Gene_Orientation) %>%
      arrange(Bacterium, `Cluster start`, `Gene start`)
    datatable(d, rownames = FALSE, selection = "none",
              options = list(pageLength = 15, scrollX = TRUE))
  })

  # --- BIOSYNTHETIC GENE CLUSTERS (antiSMASH) ---
  output$antismash_missing_note <- renderUI({
    tagList(
      if (length(antismash_bacteria_zero_regions) > 0) {
        p(style = "color:#1b7a3d; font-size: 12.5px;",
          icon("circle-check"), " antiSMASH ran clean (no BGC regions found) for: ",
          paste(antismash_bacteria_zero_regions, collapse = ", "), ".")
      } else NULL,
      if (length(antismash_bacteria_missing) > 0) {
        p(style = "color:#a33; font-size: 12.5px;",
          icon("triangle-exclamation"), " No antiSMASH run yet for: ",
          paste(antismash_bacteria_missing, collapse = ", "), ".")
      } else NULL
    )
  })

  output$antismash_summary_stats <- renderUI({
    d <- all_antismash_data
    if (nrow(d) == 0) {
      return(p(style = "color:#666;", "No antiSMASH data loaded yet."))
    }
    n_regions <- nrow(d)
    n_bact <- length(unique(d$Bacterium))
    n_known <- sum(!is.na(d$Known_Cluster_Accession) & d$Known_Cluster_Accession != "")
    fluidRow(
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(n_regions, style = "margin:0; color:#1b263b;"),
             tags$p("BGC regions found", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(length(antismash_categories), style = "margin:0; color:#c9a227;"),
             tags$p("Distinct product categories", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(n_known, style = "margin:0; color:#5b7c99;"),
             tags$p("Regions with a MIBiG hit", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(paste0(n_bact, " / 12"), style = "margin:0; color:#1b263b;"),
             tags$p("Bacteria with an antiSMASH run", style = "margin:0; color:#666; font-size:12px;")))
    )
  })

  # -- BGC region viewer (community-wide page): mirrors antiSMASH's own
  # region view -- pick a bacterium, pick one of its predicted regions, see
  # the gene-arrow diagram + region/MIBiG details, exactly like antiSMASH's
  # native results page.
  bgc_viewer_regions <- reactive({
    req(input$bgc_viewer_bact)
    f <- file.path("data", paste0(input$bgc_viewer_bact, "_antismash_regions.tsv"))
    if (!file.exists(f)) return(NULL)
    d <- tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(d) || nrow(d) == 0) return(NULL)
    d %>% arrange(Start)
  })

  observeEvent(input$bgc_viewer_bact, {
    d <- bgc_viewer_regions()
    if (is.null(d)) {
      updateSelectInput(session, "bgc_viewer_region", choices = character(0))
      return()
    }
    choices <- setNames(d$Region, paste0(d$Region, " (", d$Predicted_Product, ")"))
    updateSelectInput(session, "bgc_viewer_region", choices = choices)
  }, ignoreInit = FALSE)

  output$bgc_viewer_header <- renderUI({
    d <- bgc_viewer_regions()
    req(d, input$bgc_viewer_region)
    row <- d %>% filter(Region == input$bgc_viewer_region)
    req(nrow(row) > 0)
    row <- row[1, ]
    tags$div(
      style = "margin: 8px 0 4px 0; padding: 8px 10px; background:#f4f6f9; border-radius:4px;",
      tags$b(row$Contig, " -- ", row$Region, " -- ", row$Predicted_Product),
      tags$br(),
      tags$span(style = "color:#666; font-size:12px;",
                "Location: ", format(row$Start, big.mark = ","), " - ", format(row$End, big.mark = ","),
                " nt (", format(row$Length_bp, big.mark = ","), " bp)",
                if (!is.na(row$Category) && row$Category != "") paste0("  |  Category: ", row$Category) else "")
    )
  })

  output$bgc_viewer_plot <- renderPlot({
    req(input$bgc_viewer_bact, input$bgc_viewer_region)
    p <- render_bgc_region_plot(input$bgc_viewer_bact, input$bgc_viewer_region)
    req(p)
    p
  })

  output$bgc_viewer_mibig <- renderUI({
    d <- bgc_viewer_regions()
    req(d, input$bgc_viewer_region)
    row <- d %>% filter(Region == input$bgc_viewer_region)
    req(nrow(row) > 0)
    row <- row[1, ]
    if (is.na(row$Known_Cluster_Accession) || row$Known_Cluster_Accession == "") {
      return(p(style = "color:#888; font-size:12.5px; margin-top:6px;",
                icon("circle-info"), " No similar known cluster found in the MIBiG database."))
    }
    tags$div(
      style = "margin-top:6px; padding:8px 10px; border-left:3px solid #c9a227; background:#fbf8ef;",
      tags$b("Most similar known cluster (MIBiG): "), row$Known_Cluster_Description,
      tags$br(),
      tags$span(style = "color:#666; font-size:12px;",
                "Accession: ", row$Known_Cluster_Accession,
                "  |  Similarity: ", row$Known_Cluster_Similarity_pct, "%")
    )
  })

  # --- AMR & SPECIALTY GENES ---

  # -- Curated AMRFinderPlus results --
  output$amrfinder_missing_note <- renderUI({
    tagList(
      if (length(amrfinder_bacteria_zero_hits) > 0) {
        p(style = "color:#1b7a3d; font-size: 12.5px;",
          icon("circle-check"), " AMRFinderPlus ran clean (zero AMR elements found) for: ",
          paste(amrfinder_bacteria_zero_hits, collapse = ", "), ".")
      } else NULL,
      if (length(amrfinder_bacteria_missing) > 0) {
        p(style = "color:#a33; font-size: 12.5px;",
          icon("triangle-exclamation"), " No AMRFinderPlus run yet for: ",
          paste(amrfinder_bacteria_missing, collapse = ", "),
          ". Those bacteria fall back to the keyword screen below.")
      } else NULL
    )
  })

  output$amrfinder_summary_stats <- renderUI({
    d <- all_amrfinder_data
    if (nrow(d) == 0) {
      return(p(style = "color:#666;", "No AMRFinderPlus data loaded yet."))
    }
    n_total <- nrow(d)
    n_bact <- length(unique(d$Bacterium))
    n_classes <- length(unique(d$Class[!is.na(d$Class)]))
    n_amr_type <- sum(d$Type == "AMR", na.rm = TRUE)
    fluidRow(
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(n_total, style = "margin:0; color:#1b263b;"),
             tags$p("Total elements found", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(n_amr_type, style = "margin:0; color:#c9a227;"),
             tags$p("AMR-type hits", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(n_classes, style = "margin:0; color:#5b7c99;"),
             tags$p("Distinct drug classes", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(paste0(n_bact, " / 12"), style = "margin:0; color:#1b263b;"),
             tags$p("Bacteria with an AMRFinderPlus run", style = "margin:0; color:#666; font-size:12px;")))
    )
  })

  output$amrfinder_table <- renderDT({
    d <- all_amrfinder_data
    if (nrow(d) == 0) return(datatable(data.frame(), options = list(dom = "t")))
    if (!identical(input$amrfinder_bact_filter, "all")) {
      d <- d %>% filter(Bacterium_filename == input$amrfinder_bact_filter)
    }
    d <- d %>%
      select(Bacterium, Gene = `Element symbol`, Name = `Element name`, Type, Class, Subclass,
             `% Identity` = `% Identity to reference`, `% Coverage` = `% Coverage of reference`,
             Method, Contig = `Contig id`, Start, Stop) %>%
      arrange(Bacterium, Class)
    datatable(d, rownames = FALSE, selection = "none",
              options = list(pageLength = 15, scrollX = TRUE))
  })

  # -- Keyword screen (approximate, all 12 bacteria) --
  output$specialty_genes_summary_stats <- renderUI({
    d <- specialty_genes_data
    if (nrow(d) == 0) {
      return(p(style = "color:#666;", "No hits found by this keyword screen."))
    }
    n_total <- nrow(d)
    n_amr <- sum(d$Category == "AMR")
    n_vir <- sum(d$Category == "Virulence")
    n_bact <- length(unique(d$Bacterium))
    fluidRow(
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(n_total, style = "margin:0; color:#1b263b;"),
             tags$p("Total hits", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(n_amr, style = "margin:0; color:#c9a227;"),
             tags$p("AMR-family hits", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(n_vir, style = "margin:0; color:#5b7c99;"),
             tags$p("Virulence-associated hits", style = "margin:0; color:#666; font-size:12px;"))),
      column(width = 3, tags$div(style = "text-align:center; padding:10px;",
             tags$h3(n_bact, style = "margin:0; color:#1b263b;"),
             tags$p("Bacteria with at least one hit (of 12)", style = "margin:0; color:#666; font-size:12px;")))
    )
  })

  output$specialty_genes_heatmap <- renderPlotly({
    d <- specialty_genes_data
    if (nrow(d) == 0) {
      return(plotly_empty() %>%
               layout(annotations = list(x = 0.5, y = 0.5, text = "No hits found by this keyword screen.",
                                          xref = "paper", yref = "paper", showarrow = FALSE)))
    }
    mat <- d %>% count(Bacterium, Subcategory, name = "n") %>%
      tidyr::pivot_wider(names_from = Subcategory, values_from = n, values_fill = 0)
    bact_order <- mat$Bacterium
    cols <- setdiff(colnames(mat), "Bacterium")
    z <- as.matrix(mat[, cols, drop = FALSE])
    rownames(z) <- bact_order
    plot_ly(x = cols, y = bact_order, z = z, type = "heatmap",
            colors = colorRamp(c("#f4f6f9", "#c9a227", "#0d1b2a")),
            hovertemplate = "%{y} | %{x}: %{z}<extra></extra>") %>%
      layout(xaxis = list(title = "", tickangle = -45),
             yaxis = list(title = "", automargin = TRUE),
             margin = list(l = 10, b = 160))
  })

  output$specialty_genes_table <- renderDT({
    d <- specialty_genes_data
    if (nrow(d) == 0) return(datatable(data.frame(), options = list(dom = "t")))
    if (!identical(input$specialty_genes_category_filter, "all")) {
      d <- d %>% filter(Category == input$specialty_genes_category_filter)
    }
    if (!identical(input$specialty_genes_bact_filter, "all")) {
      d <- d %>% filter(Filename == input$specialty_genes_bact_filter)
    }
    d <- d %>%
      select(Bacterium, Category, Subcategory, Gene, Product, `Locus Tag` = Locus_Tag) %>%
      arrange(Bacterium, Category, Subcategory)
    datatable(d, rownames = FALSE, selection = "none",
              options = list(pageLength = 15, scrollX = TRUE))
  })

  output$amr_phenotype_content <- renderUI({
    req(input$amr_phenotype_bact)
    pheno <- load_amr_phenotype(input$amr_phenotype_bact)
    bact_name <- bacteria_names[match(input$amr_phenotype_bact, bacteria_filenames)]
    if (is.null(pheno) || nrow(pheno) == 0) {
      return(tagList(
        p(icon("exclamation-triangle"),
          paste0(" No lab-tested AMR phenotype data on file for ", bact_name, ". To add it:")),
        tags$ol(
          tags$li("Create a TSV file named: ",
                  tags$code(paste0(input$amr_phenotype_bact, "_amr_phenotype.tsv"))),
          tags$li("Include columns such as: antibiotic, method (MIC/disk-diffusion), result (S/I/R), mic_value"),
          tags$li("Place the file in the 'data/' folder and reload the page.")
        )
      ))
    }
    DTOutput("amr_phenotype_table")
  })

  output$amr_phenotype_table <- renderDT({
    pheno <- load_amr_phenotype(input$amr_phenotype_bact)
    req(pheno)
    datatable(pheno, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  # --- COMPARE TWO BACTERIA ---

  compare_bact_name <- function(filename) {
    idx <- match(filename, bacteria_filenames)
    if (is.na(idx)) filename else bacteria_names[idx]
  }

  output$compare_genome_stats_table <- renderDT({
    req(input$compare_bact_a, input$compare_bact_b)
    if (is.null(genome_info_data) || !("bacterium_name_clean" %in% colnames(genome_info_data))) {
      return(datatable(data.frame(Message = "genome_info.tsv not available."), rownames = FALSE, options = list(dom = 't')))
    }
    row_a <- genome_info_data %>% filter(bacterium_name_clean == input$compare_bact_a)
    row_b <- genome_info_data %>% filter(bacterium_name_clean == input$compare_bact_b)

    fields  <- c("size_bp", "gc_content", "num_genes", "assembly_level", "shape", "lifestyle", "type")
    labels  <- c("Genome Size (bp)", "GC Content (%)", "Number of Genes", "Assembly Level", "Shape", "Lifestyle", "Type")
    get_val <- function(row, col) {
      if (nrow(row) > 0 && col %in% colnames(row) && !is.na(row[[col]][1]) && nzchar(as.character(row[[col]][1]))) {
        as.character(row[[col]][1])
      } else "NA"
    }

    df <- data.frame(
      Metric = labels,
      A = vapply(fields, function(f) get_val(row_a, f), character(1)),
      B = vapply(fields, function(f) get_val(row_b, f), character(1)),
      stringsAsFactors = FALSE
    )
    colnames(df) <- c("Metric", compare_bact_name(input$compare_bact_a), compare_bact_name(input$compare_bact_b))
    datatable(df, rownames = FALSE, options = list(dom = 't', pageLength = 20))
  })

  output$compare_cog_chart <- renderPlotly({
    req(input$compare_bact_a, input$compare_bact_b)
    validate(need(input$compare_bact_a != input$compare_bact_b, "Select two different bacteria to compare."))
    a_name <- compare_bact_name(input$compare_bact_a)
    b_name <- compare_bact_name(input$compare_bact_b)
    combined <- bind_rows(
      compute_cog_counts(input$compare_bact_a) %>% mutate(Bacterium = a_name),
      compute_cog_counts(input$compare_bact_b) %>% mutate(Bacterium = b_name)
    )
    validate(need(nrow(combined) > 0, "No COG data available for one or both bacteria."))
    plot_ly(combined, x = ~COG_category, y = ~Count, color = ~Bacterium, type = "bar",
            colors = c("#1b263b", "#c9a227")) %>%
      layout(barmode = "group", xaxis = list(title = "COG Category"), yaxis = list(title = "Genes"),
             legend = list(orientation = "h", x = 0, y = 1.15))
  })

  output$compare_localization_chart <- renderPlotly({
    req(input$compare_bact_a, input$compare_bact_b)
    validate(need(input$compare_bact_a != input$compare_bact_b, "Select two different bacteria to compare."))
    a_name <- compare_bact_name(input$compare_bact_a)
    b_name <- compare_bact_name(input$compare_bact_b)
    combined <- bind_rows(
      compute_localization_counts(input$compare_bact_a) %>% mutate(Bacterium = a_name),
      compute_localization_counts(input$compare_bact_b) %>% mutate(Bacterium = b_name)
    )
    validate(need(nrow(combined) > 0, "No localization data available for one or both bacteria."))
    plot_ly(combined, x = ~Localization, y = ~Count, color = ~Bacterium, type = "bar",
            colors = c("#1b263b", "#c9a227")) %>%
      layout(barmode = "group", xaxis = list(title = ""), yaxis = list(title = "Proteins"),
             legend = list(orientation = "h", x = 0, y = 1.15))
  })

  output$compare_mobileog_chart <- renderPlotly({
    req(input$compare_bact_a, input$compare_bact_b)
    validate(need(input$compare_bact_a != input$compare_bact_b, "Select two different bacteria to compare."))
    a_name <- compare_bact_name(input$compare_bact_a)
    b_name <- compare_bact_name(input$compare_bact_b)
    combined <- bind_rows(
      compute_mobileog_counts(input$compare_bact_a) %>% mutate(Bacterium = a_name),
      compute_mobileog_counts(input$compare_bact_b) %>% mutate(Bacterium = b_name)
    )
    validate(need(nrow(combined) > 0, "No mobileOG data available for one or both bacteria."))
    plot_ly(combined, x = ~Category, y = ~Count, color = ~Bacterium, type = "bar",
            colors = c("#1b263b", "#c9a227")) %>%
      layout(barmode = "group", xaxis = list(title = ""), yaxis = list(title = "Genes"),
             legend = list(orientation = "h", x = 0, y = 1.15))
  })

  # Real sequence-based ortholog overlap: reciprocal best hits between the two
  # genomes from data/ortholog_best_hits.tsv (k-mer containment similarity,
  # computed offline across all 12 proteomes -- see Gene Neighborhood page for
  # the same underlying method and its validation).
  compare_shared_genes_data <- reactive({
    req(input$compare_bact_a, input$compare_bact_b)
    req(input$compare_bact_a != input$compare_bact_b)
    req(!is.null(ortholog_hits), nrow(ortholog_hits) > 0)

    a_total <- all_genes_index %>% filter(Filename == input$compare_bact_a) %>% nrow()
    b_total <- all_genes_index %>% filter(Filename == input$compare_bact_b) %>% nrow()

    info_a <- all_genes_index %>% filter(Filename == input$compare_bact_a) %>%
      select(Locus_Tag, Gene, Product)
    info_b <- all_genes_index %>% filter(Filename == input$compare_bact_b) %>%
      select(Locus_Tag, Gene, Product)

    pairs <- ortholog_hits %>%
      filter(source_bacterium == input$compare_bact_a, target_bacterium == input$compare_bact_b,
             reciprocal == TRUE) %>%
      left_join(info_a, by = c("source_locus_tag" = "Locus_Tag")) %>%
      rename(Gene_A = Gene, Product_A = Product) %>%
      left_join(info_b, by = c("target_locus_tag" = "Locus_Tag")) %>%
      rename(Gene_B = Gene, Product_B = Product)

    list(
      a_total = a_total,
      b_total = b_total,
      n_shared = nrow(pairs),
      table = pairs %>%
        transmute(
          `Gene (A)` = Gene_A, `Product (A)` = Product_A, `Locus Tag (A)` = source_locus_tag,
          `Gene (B)` = Gene_B, `Product (B)` = Product_B, `Locus Tag (B)` = target_locus_tag,
          Score = score
        ) %>%
        arrange(desc(Score))
    )
  })

  output$compare_shared_summary <- renderUI({
    req(input$compare_bact_a, input$compare_bact_b)
    if (input$compare_bact_a == input$compare_bact_b) {
      return(tags$p(style = "color:#888;", "Select two different bacteria to compare."))
    }
    info <- compare_shared_genes_data()
    tags$p(style = "font-size: 14px;",
      sprintf("%d reciprocal-best-hit orthologs found between the two (out of %d genes in %s and %d genes in %s).",
              info$n_shared, info$a_total, compare_bact_name(input$compare_bact_a),
              info$b_total, compare_bact_name(input$compare_bact_b))
    )
  })

  output$compare_shared_genes_table <- renderDT({
    req(input$compare_bact_a, input$compare_bact_b, input$compare_bact_a != input$compare_bact_b)
    info <- compare_shared_genes_data()
    req(info$n_shared > 0)
    datatable(info$table, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  # --- GENE NEIGHBORHOOD / SYNTENY VIEW ---

  # Anchor gene choices depend on which bacterium is selected; repopulate
  # (server-side, since a genome can have a few thousand named genes) whenever
  # the anchor bacterium changes.
  observeEvent(input$neighborhood_bact, {
    req(input$neighborhood_bact)
    choices <- get_gene_choices_for(input$neighborhood_bact)
    updateSelectizeInput(session, "neighborhood_gene", choices = choices, server = TRUE, selected = "")
  }, ignoreInit = FALSE)

  # Resolve the (bacterium, gene symbol) selection to a specific locus tag --
  # the actual anchor key the ortholog table and neighborhood builder use.
  neighborhood_anchor_locus <- reactive({
    req(input$neighborhood_bact, input$neighborhood_gene, nzchar(input$neighborhood_gene))
    match_row <- all_genes_index %>%
      filter(Filename == input$neighborhood_bact, Gene == input$neighborhood_gene) %>%
      slice(1)
    if (nrow(match_row) == 0) return(NULL)
    match_row$Locus_Tag[1]
  })

  output$neighborhood_summary <- renderUI({
    req(input$neighborhood_bact, input$neighborhood_gene, nzchar(input$neighborhood_gene))
    locus_tag <- neighborhood_anchor_locus()
    if (is.null(locus_tag)) {
      return(tags$p(style = "color:#888;", "Gene not found in the selected bacterium."))
    }
    window <- input$neighborhood_window
    if (is.null(window) || is.na(window)) window <- 3
    neigh <- build_neighborhood_data(input$neighborhood_bact, locus_tag, window)
    if (nrow(neigh) == 0) {
      return(tags$p(style = "color:#888;", "No neighborhood data available."))
    }
    n_bact <- length(unique(neigh$BactIndex))
    n_recip <- length(unique(neigh$BactIndex[neigh$IsAnchor | neigh$OrthologReciprocal]))
    tags$p(style = "font-size: 13px; color:#666;",
           sprintf("Ortholog found in %d of 12 community members (%d via reciprocal best hit).",
                   n_bact, n_recip))
  })

  neighborhood_result <- reactive({
    req(input$neighborhood_bact, input$neighborhood_gene, nzchar(input$neighborhood_gene))
    locus_tag <- neighborhood_anchor_locus()
    req(!is.null(locus_tag))
    window <- input$neighborhood_window
    if (is.null(window) || is.na(window)) window <- 3
    render_gene_neighborhood_plot(input$neighborhood_bact, locus_tag, window)
  })

  # gggenes' arrows need real vertical room per facet row, so size the plot
  # to the number of bacteria found rather than using a fixed height.
  output$neighborhood_plot_ui <- renderUI({
    result <- neighborhood_result()
    validate(need(!is.null(result), "No neighborhood data available for this gene."))
    plotOutput("neighborhood_plot", height = paste0(max(320, 95 * result$n_bacteria + 90), "px"))
  })

  output$neighborhood_plot <- renderPlot({
    result <- neighborhood_result()
    validate(need(!is.null(result), "No neighborhood data available for this gene."))
    result$plot
  })

  output$neighborhood_genes_table <- renderDT({
    req(input$neighborhood_bact, input$neighborhood_gene, nzchar(input$neighborhood_gene))
    locus_tag <- neighborhood_anchor_locus()
    req(!is.null(locus_tag))
    window <- input$neighborhood_window
    if (is.null(window) || is.na(window)) window <- 3
    df <- build_neighborhood_display_data(input$neighborhood_bact, locus_tag, window)
    req(!is.null(df))

    display <- df %>%
      mutate(Bacterium = as.character(Bacterium),
             Product = ifelse(is.na(Product) | Product == "", Type, Product)) %>%
      transmute(Bacterium, Gene, Product, `Locus Tag` = Locus_Tag, Strand,
                `Rank vs. anchor` = RelPos, Role = as.character(role)) %>%
      arrange(Bacterium, `Rank vs. anchor`)
    datatable(display, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })

  # --- PAN-GENOME / ORTHOLOG MATRIX ---

  output$pangenome_summary_stats <- renderUI({
    req(!is.null(orthogroup_summary_data))
    counts <- orthogroup_summary_data %>% count(classification)
    get_n <- function(cls) { v <- counts$n[counts$classification == cls]; if (length(v) == 0) 0 else v }
    total_og <- nrow(orthogroup_summary_data)

    stat_box <- function(label, value) {
      column(width = 3,
             div(style = "text-align:center; padding: 10px;",
                 tags$div(style = "font-size: 26px; font-weight: 700; color: var(--navy-800, #1b263b);", value),
                 tags$div(style = "font-size: 12px; color:#666; text-transform: uppercase; letter-spacing: 0.5px;", label)
             )
      )
    }
    fluidRow(
      stat_box("Orthogroups", format(total_og, big.mark = ",")),
      stat_box("Core (12/12)", get_n("core")),
      stat_box("Soft-core (9-11)", get_n("soft-core")),
      stat_box("Shell (2-8)", format(get_n("shell"), big.mark = ","))
    )
  })

  output$pangenome_composition_chart <- renderPlotly({
    req(!is.null(pangenome_composition))
    df <- pangenome_composition %>%
      mutate(BacteriumName = bacteria_names[match(bacterium, bacteria_filenames)],
             classification = factor(classification, levels = c("core", "soft-core", "shell", "unique")))
    plot_ly(df, x = ~n, y = ~BacteriumName, color = ~classification, type = "bar", orientation = "h",
            colors = c("core" = "#1b263b", "soft-core" = "#2c3e57", "shell" = "#c9a227", "unique" = "#a8b4c4")) %>%
      layout(barmode = "stack", xaxis = list(title = "Genes"), yaxis = list(title = ""),
             legend = list(orientation = "h", x = 0, y = 1.1), margin = list(l = 220))
  })

  pangenome_filtered_orthogroups <- reactive({
    req(!is.null(orthogroup_summary_full))
    df <- orthogroup_summary_full
    if (!is.null(input$pangenome_class_filter) && input$pangenome_class_filter != "all") {
      df <- df %>% filter(classification == input$pangenome_class_filter)
    }
    df %>% arrange(desc(n_bacteria), desc(n_genes))
  })

  output$pangenome_orthogroup_table <- renderDT({
    df <- pangenome_filtered_orthogroups() %>%
      transmute(
        `Orthogroup` = orthogroup_id,
        `Representative Gene` = ifelse(is.na(representative_gene), "-", representative_gene),
        `Sample Product` = ifelse(is.na(sample_product), "-", sample_product),
        `# Bacteria` = n_bacteria,
        `# Genes` = n_genes,
        `Classification` = classification
      )
    datatable(df, rownames = FALSE, selection = "single",
              options = list(pageLength = 15, scrollX = TRUE))
  })

  selected_orthogroup_id <- reactiveVal(NULL)

  observeEvent(input$pangenome_orthogroup_table_rows_selected, {
    sel <- input$pangenome_orthogroup_table_rows_selected
    req(sel)
    df <- pangenome_filtered_orthogroups()
    req(nrow(df) >= sel)
    selected_orthogroup_id(df$orthogroup_id[sel])
  })

  output$pangenome_selected_og_title <- renderUI({
    if (is.null(selected_orthogroup_id())) {
      "Gene membership (select an orthogroup above)"
    } else {
      paste("Gene membership:", selected_orthogroup_id())
    }
  })

  output$pangenome_orthogroup_members_table <- renderDT({
    req(!is.null(selected_orthogroup_id()), !is.null(orthogroup_genes_data))
    df <- orthogroup_genes_data %>%
      filter(orthogroup_id == selected_orthogroup_id()) %>%
      mutate(BacteriumName = bacteria_names[match(bacterium, bacteria_filenames)]) %>%
      transmute(
        Bacterium = ifelse(is.na(BacteriumName), bacterium, BacteriumName),
        Gene = ifelse(is.na(gene) | gene == "", "-", gene),
        Product = ifelse(is.na(product) | product == "", "-", product),
        `Locus Tag` = locus_tag
      ) %>%
      arrange(Bacterium)
    datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE, dom = 'tip'))
  })

  orthogroup_tree_result <- reactive({
    req(selected_orthogroup_id())
    build_orthogroup_tree(selected_orthogroup_id())
  })

  output$pangenome_orthogroup_tree_ui <- renderUI({
    if (is.null(selected_orthogroup_id())) {
      return(p(style = "color:#666; font-size: 13px;", "Select an orthogroup above to see its gene tree."))
    }
    result <- orthogroup_tree_result()
    if (is.null(result$plot)) {
      return(p(style = "color:#a33; font-size: 13px;", result$error))
    }
    tagList(
      if (!is.null(result$error)) p(style = "color:#a33; font-size: 12px;", result$error) else NULL,
      girafeOutput("pangenome_orthogroup_tree", height = paste0(max(220, 40 * result$n_members + 60), "px"))
    )
  })

  output$pangenome_orthogroup_tree <- renderGirafe({
    result <- orthogroup_tree_result()
    req(result$plot)
    girafe(
      ggobj = result$plot,
      width_svg = 7,
      height_svg = max(2.5, 0.45 * result$n_members + 0.6),
      options = list(
        opts_zoom(min = 0.5, max = 6),
        opts_hover(css = "cursor:pointer;fill:#c9a227;stroke:#c9a227;stroke-width:2px;", reactive = TRUE),
        opts_hover_inv(css = "opacity:0.35;"),
        opts_tooltip(
          css = "background-color:#1f2937;color:#f9fafb;padding:10px;border-radius:6px;font-size:13px;",
          opacity = 0.95, use_fill = FALSE, use_stroke = FALSE
        ),
        opts_sizing(rescale = TRUE, width = 1)
      )
    )
  })

  # --- OBSERVERS: BUTTONS & CLICKS ---
  
  # Bacteria Selection Observer (resetting)
  lapply(1:12, function(i) {
    shinyjs::onclick(
      id = paste0("bacteria", i, "-img"),
      expr = {
        current_bacterium(list(name = bacteria_names[i], filename = bacteria_filenames[i], index = i))
        updateTabItems(session, "sidebar_menu", selected = "details_page")
        active_content_type(NULL) 
        selected_localization_type(NULL)
        selected_annotation_type(NULL)
        selected_mobileog_category(NULL)
        selected_cog_category(NULL)
        selected_clade(NULL)
      }
    )
  })
  
  # Back Button
  observeEvent(input$back_to_collection, {
    updateTabItems(session, "sidebar_menu", selected = "omm12")
    active_content_type(NULL)
  })
  
  # NEW: Observer for phylogeny plot interactions
  # COMMENTED OUT for now -- Similarity Data & Strain Information feature paused.
  # Restore by uncommenting; depends on current_bacterium_similarity_data() below.
  # observeEvent(input$phylogeny_plot_enhanced_selected, {
  #   tip_clicked <- input$phylogeny_plot_enhanced_selected
  #   similarity_data <- current_bacterium_similarity_data()
  #
  #   if (!is.null(tip_clicked) && tip_clicked %in% current_bacterium()$filename) {
  #     # Check if similarity data exists
  #     if (!is.null(similarity_data) && nrow(similarity_data) > 0) {
  #       # Find the correct column to match tip labels
  #       join_col <- intersect(c("strain_name","label","name"), colnames(similarity_data))[1]
  #
  #       # Filter similarity info for the clicked tip
  #       tip_info <- similarity_data %>% filter(.data[[join_col]] == tip_clicked)
  #
  #       # Show modal popup with the data
  #       showModal(modalDialog(
  #         title = paste0("Strain Information: ", tip_clicked),
  #         DT::datatable(tip_info, options = list(scrollX = TRUE, pageLength = 5)),
  #         easyClose = TRUE,
  #         size = "l"
  #       ))
  #     } else {
  #       showModal(modalDialog(
  #         title = paste0("Strain Information: ", tip_clicked),
  #         "No similarity data available for this strain.",
  #         easyClose = TRUE
  #       ))
  #     }
  #   }
  # })
  
  # Plotly click observers
  observeEvent(event_data("plotly_click", source = "cog_barplot_clicks"), {
    ed <- event_data("plotly_click", source = "cog_barplot_clicks")
    if (!is.null(ed) && "x" %in% names(ed)) {
      selected_letter <- ed$x[1]
      selected_cog_category(trimws(selected_letter))
      legend <- cog_legend()
      if (!is.null(legend)) {
        desc <- legend %>%
          filter(Letter == selected_letter) %>%
          pull(Description) %>%
          head(1)
        display_text <- paste0(selected_letter, ": ", desc)
        shinyjs::runjs(paste0('$("#cog-selected-type").text("', display_text, '");'))
      }
    } else {
      selected_cog_category(NULL)
      shinyjs::runjs('$("#cog-selected-type").text("(Select a COG category from the bar chart)");')
    }
  })
  
  observeEvent(event_data("plotly_click", source = "localization_cog_heatmap_clicks"), {
    ed <- event_data("plotly_click", source = "localization_cog_heatmap_clicks")
    if (!is.null(ed) && "y" %in% names(ed)) {
      selected_label <- trimws(as.character(ed$y[1]))
      if (nzchar(selected_label)) {
        selected_localization_type(selected_label)
        shinyjs::runjs(paste0('$("#localization-selected-type").text("', selected_label, '");'))
      }
    }
  })

  observeEvent(event_data("plotly_click", source = "localization_treemap_clicks"), {
    ed <- event_data("plotly_click", source = "localization_treemap_clicks")
    if (!is.null(ed) && "x" %in% names(ed)) {
      selected_label <- ed$x[1]
      if (tolower(selected_label) != 'total' && selected_label != "") {
        selected_localization_type(trimws(selected_label))
        shinyjs::runjs(paste0('$("#localization-selected-type").text("', trimws(selected_label), '");'))
      } else {
        selected_localization_type(NULL)
        shinyjs::runjs('$("#localization-selected-type").text("(Select a location from the chart)");')
      }
    } else {
      selected_localization_type(NULL)
      shinyjs::runjs('$("#localization-selected-type").text("(Select a location from the chart)");')
    }
  })
  
  observeEvent(event_data("plotly_click", source = "annotation_pie_chart_clicks"), {
    ed <- event_data("plotly_click", source = "annotation_pie_chart_clicks")
    if (!is.null(ed) && "pointNumber" %in% names(ed)) {
      selected_index <- ed$pointNumber + 1
      stats <- annotation_stats()
      if (selected_index <= nrow(stats)) {
        selected_label <- stats$type[selected_index]
        selected_annotation_type(selected_label)
        shinyjs::runjs(paste0('$("#annotation-selected-type").text("', selected_label, '");'))
      } else {
        selected_annotation_type(NULL)
        shinyjs::runjs('$("#annotation-selected-type").text("Select a type from the chart...");')
      }
    } else {
      selected_annotation_type(NULL)
      shinyjs::runjs('$("#annotation-selected-type").text("Select a type from the chart...");')
    }
  })
  
  observeEvent(event_data("plotly_click", source = "mobileog_donut_clicks"), {
    ed <- event_data("plotly_click", source = "mobileog_donut_clicks")
    if (!is.null(ed) && "pointNumber" %in% names(ed)) {
      selected_index <- ed$pointNumber + 1
      category_counts <- mobileog_category_counts()
      if (selected_index <= nrow(category_counts)) {
        selected_label <- category_counts$Category[selected_index]
        selected_mobileog_category(trimws(selected_label))
        shinyjs::runjs(paste0('$("#mobileog-selected-type").text("', trimws(selected_label), '");'))
      } else {
        selected_mobileog_category(NULL)
        shinyjs::runjs('$("#mobileog-selected-type").text("Select a category from the chart...");')
      }
    } else {
      selected_mobileog_category(NULL)
      shinyjs::runjs('$("#mobileog-selected-type").text("Select a category from the chart...");')
    }
  })
  
  observeEvent(event_data("plotly_click", source = "defense_donut_clicks"), {
    ed <- event_data("plotly_click", source = "defense_donut_clicks")
    if (!is.null(ed) && "pointNumber" %in% names(ed)) {
      idx <- ed$pointNumber + 1
      counts <- defense_type_counts()
      if (idx <= nrow(counts)) {
        lbl <- counts$Type[idx]
        selected_defense_type(trimws(lbl))
        shinyjs::runjs(paste0('$("#defense-selected-type").text("', trimws(lbl), '");'))
      } else {
        selected_defense_type(NULL)
        shinyjs::runjs('$("#defense-selected-type").text("Select a system from the chart...");')
      }
    } else {
      selected_defense_type(NULL)
      shinyjs::runjs('$("#defense-selected-type").text("Select a system from the chart...");')
    }
  })
  
  selected_kegg_node <- reactiveVal(NULL)
  kegg_node_names    <- reactiveVal(NULL) 
  
  observeEvent(event_data("plotly_click", source = "kegg_sankey_clicks"), {
    ed <- event_data("plotly_click", source = "kegg_sankey_clicks")
    req(ed, "pointNumber" %in% names(ed))
    names_vec <- kegg_node_names()
    req(names_vec)
    idx <- ed$pointNumber[1] + 1          # plotly is 0-indexed
    if (idx >= 1 && idx <= length(names_vec)) {
      selected_kegg_node(names_vec[idx])
    }
  })
  # in the nav observer
  observeEvent(input$nav_metabolic_model, {
    cat(">>> NAV metabolic_model clicked; bacterium =",
        if (!is.null(current_bacterium())) current_bacterium()$filename else "NONE", "\n")
    active_content_type("metabolic_model")
    send_scroll_message("metabolic_model_box")
  })
  
    # in current_gem
  current_gem <- reactive({
    req(current_bacterium())
    fn <- current_bacterium()$filename
    path <- file.path("data", paste0(fn, "_model.xml"))
    cat(">>> GEM looking for:", path, "| exists =", file.exists(path), "\n")
    g <- load_gem_model(fn)
    cat(">>> GEM loaded; is.null =", is.null(g),
        if (!is.null(g)) paste("| reactions =", nrow(g$reactions),
                               "| subsystems =", nrow(g$subsystems)) else "", "\n")
    g
  })
  
  #observeEvent(list(current_gem(), active_content_type()), {
   # cat(">>> selectize observer fired; active =", active_content_type(), "\n")
    #req(active_content_type() == "metabolic_model")
    #gem <- current_gem()
    #req(gem, nrow(gem$subsystems) > 0)
    #ch <- gem_subsystem_choices()
    #req(length(ch) > 0)
    #updateSelectizeInput(session, "gem_subsystem",
     #                    choices = ch, selected = ch[1], server = TRUE)
    #cat(">>> updateSelectizeInput called\n")
  #}, ignoreInit = FALSE)
  # Side Panel Navigation Observers
  send_scroll_message <- function(target_id) {
    session$sendCustomMessage(type = "scrollTo", message = list(id = target_id))
  }
  
  observeEvent(input$nav_general_info, { 
    active_content_type(NULL)
    send_scroll_message("General_Information_Section") 
  })
  observeEvent(input$nav_genomic_info, { 
    active_content_type(NULL)
    send_scroll_message("Genomic_Information_Section") 
  })
  observeEvent(input$nav_proteomic_info, { 
    active_content_type(NULL)
    send_scroll_message("Proteomic_Information_Section") 
  })
  observeEvent(input$nav_metabolic_info, { 
    active_content_type("metabolic")
    send_scroll_message("Metabolic_Information_Section") 
  })
  observeEvent(input$nav_phylogeny, { 
    active_content_type("phylogeny") 
    send_scroll_message("phylogeny_content_ui_box") 
  })
  observeEvent(input$nav_genome_map, { 
    active_content_type("genome") 
    send_scroll_message("genome_plot_ui_box") 
  })
  
  observeEvent(input$nav_cog, { 
    active_content_type("cog") 
    send_scroll_message("cog_content_box") 
  })
  observeEvent(input$nav_kegg, { 
    active_content_type("kegg") 
    send_scroll_message("kegg_content_box") 
  })
  
    observeEvent(input$nav_mobileog, { 
    active_content_type("mobileog") 
    send_scroll_message("mobileog_content_box") 
  })
  observeEvent(input$nav_localization, { 
    active_content_type("localization") 
    send_scroll_message("localization_content_box") 
  })
  observeEvent(input$nav_defense, {
    active_content_type("defense")
    send_scroll_message("defense_content_box")
  })
  observeEvent(input$nav_crispr, {
    active_content_type("crispr")
    send_scroll_message("crispr_content_box")
  })
  observeEvent(input$nav_bgc, {
    active_content_type("bgc")
    send_scroll_message("bgc_content_box")
  })
  observeEvent(input$nav_specialty_genes, {
    active_content_type("specialty_genes")
    send_scroll_message("specialty_genes_content_box")
  })
  observeEvent(input$defense_activity, {
    selected_defense_type(NULL)
    shinyjs::runjs('$("#defense-selected-type").text("Select a system from the chart...");')
  })
  
  observeEvent(input$search_gene, {
    req(active_content_type() == "genome")
    
    bacterium <- current_bacterium()
    req(bacterium)
    
    slider_input_id <- paste0("region_", bacterium$filename)
    searched_gene <- input$search_gene
    if (is.null(searched_gene) || searched_gene == "") return()
    
    genes_df <- current_bacterium_genes()
    req(nrow(genes_df) > 0)
    
    matched_genes <- genes_df %>%
      filter(
        contains_ci(gene, searched_gene) |
          contains_ci(description, searched_gene)
      )
    
    if (nrow(matched_genes) == 0) return()
    
    target_start <- matched_genes$start[1]
    target_end   <- matched_genes$end[1]
    
    buffer <- 5000
    new_start <- max(1, target_start - buffer)
    new_end   <- min(max(genes_df$end), target_end + buffer)
    
    updateSliderInput(
      session,
      slider_input_id,
      value = c(new_start, new_end)
    )
  })
  
  observeEvent(current_bacterium(), {
    req(active_content_type() == "genome")
    
    genes_df <- current_bacterium_genes()
    req(genes_df)
    
    suggestions <- unique(c(genes_df$gene, genes_df$description))
    
    updateSelectizeInput(
      session,
      "search_gene",
      choices = suggestions,
      server = FALSE
    )
  })
  
  # Copy data/{bacterium}{suffix} into a download; if the file is missing, send a
# short text note instead of failing silently.
copy_data_file <- function(base_name, suffix, dest) {
  src <- file.path("data", paste0(base_name, suffix))
  if (file.exists(src)) {
    file.copy(src, dest, overwrite = TRUE)
  } else {
    writeLines(paste0("This file is not available yet for ", base_name,
                      " (expected data/", base_name, suffix, ")."), dest)
  }
}

get_base_filename <- reactive({
    req(current_bacterium())
    return(current_bacterium()$filename)
  })
  # --- UI RENDERERS (PARENT) ---
  
  output$bacteria_details_ui <- renderUI({
    req(current_bacterium())
    bacterium <- current_bacterium()
    
    tagList(
      actionButton("back_to_collection", "Back to OMM12 Collection", class = "btn-primary"),
      h2(bacterium$name, class = "details-page-header"),
      
      tags$div(id = "General_Information_Section",
               box(
                 title = "General Information", width = 12, solidHeader = TRUE, status = "info",
                 fluidRow(
                   column(width = 8, get_general_info_html(bacterium$filename)),
                   column(width = 4,
                          div(style = "text-align: center;",
                              # "?v=<file time>" makes browsers reload the picture whenever it is regenerated
                              tags$img(src = paste0("images/", bacterium$filename, ".png?v=",
                                                    as.integer(file.mtime(file.path("www", "images", paste0(bacterium$filename, ".png"))))),
                                       alt = paste("Taxonomic lineage of", bacterium$name),
                                       style = "width: 100%; max-width: 300px; aspect-ratio: 1000 / 1420; object-fit: contain; border: 1px solid #ddd; padding: 5px;"),
                              taxonomy_source_note(bacterium$filename)))
                 )
               ),
               box(
                 title = "Data Provenance", width = 12, solidHeader = TRUE, status = "primary",
                 collapsible = TRUE, collapsed = TRUE,
                 p(style = "color:#666; font-size: 12px;",
                   "Which tool produced each analysis on this page, and when it was last generated."),
                 get_provenance_html(bacterium$filename)
               )
      ),
      tags$hr(),
      
      if (isTRUE(active_content_type() %in% c("phylogeny", "genome", "cog",
                                              "kegg", "mobileog", "localization",
                                              "defense", "crispr", "bgc", "specialty_genes"))) {
        tags$div(id = "Genomic_Information_Section",
                 tags$h3("Genomic Information", class = "section-divider"))
      } else { NULL },
      
      # conditional sections
      if (isTRUE(active_content_type() == "phylogeny")) {
        tagList(
          tags$div(id = "phylogeny_content_ui_box", 
                   tags$h4("Phylogeny"), 
                   uiOutput("phylogeny_content_ui")
          ), 
          tags$hr()
        )
      } else { NULL },
      
      if (isTRUE(active_content_type() == "genome")) {
        tagList(
          tags$div(id = "genome_plot_ui_box",
                   tags$h4("Genome Map"),
                   uiOutput("genome_plot_ui")
          ),
          tags$hr()
        )
      } else { NULL },
      if (isTRUE(active_content_type() == "kegg")) {
        tagList(
          tags$div(id = "kegg_content_box", 
                   tags$h4("KEGG Pathways"), 
                   uiOutput("kegg_content")
          ), 
          tags$hr()
        )
      } else { NULL },
      
      if (isTRUE(active_content_type() == "annotation")) {
        tagList(
          tags$div(id = "annotation_content_box", 
                   tags$h4("Prokka Annotation"), 
                   uiOutput("annotation_content")
          ), 
          tags$hr()
        )
      } else { NULL },
      
      if (isTRUE(active_content_type() == "cog")) {
        tagList(
          tags$div(id = "cog_content_box", 
                   tags$h4("COG"), 
                   uiOutput("cog_content")
          ), 
          tags$hr()
        )
      } else { NULL },
      
      if (isTRUE(active_content_type() == "mobileog")) {
        tagList(
          tags$div(id = "mobileog_content_box", 
                   tags$h4("Bacterial Mobile Genetic Elements"), 
                   uiOutput("mobileog_content")
          ), 
          tags$hr()
        )
      } else { NULL },
      
      if (isTRUE(active_content_type() == "localization")) {
        tagList(
          tags$div(id = "localization_content_box", 
                   tags$h4("Localization"), 
                   uiOutput("localization_content")
          ), 
          tags$hr()
        )
      } else { NULL },
      if (isTRUE(active_content_type() == "defense")) {
        tagList(
          tags$div(id = "defense_content_box",
                   tags$h4("Defense Systems"),
                   uiOutput("defense_content")
          ),
          tags$hr()
        )
      } else { NULL },
      if (isTRUE(active_content_type() == "crispr")) {
        tagList(
          tags$div(id = "crispr_content_box",
                   tags$h4("CRISPR & Cas Systems"),
                   uiOutput("crispr_bact_content")
          ),
          tags$hr()
        )
      } else { NULL },
      if (isTRUE(active_content_type() == "bgc")) {
        tagList(
          tags$div(id = "bgc_content_box",
                   tags$h4("Biosynthetic Gene Clusters"),
                   uiOutput("bgc_bact_content")
          ),
          tags$hr()
        )
      } else { NULL },
      if (isTRUE(active_content_type() == "specialty_genes")) {
        tagList(
          tags$div(id = "specialty_genes_content_box",
                   tags$h4("AMR & Specialty Genes"),
                   uiOutput("specialty_genes_bact_content")
          ),
          tags$hr()
        )
      } else { NULL },

      #tags$div(id = "Proteomic_Information_Section", 
       #        tags$h3("Proteomic Information", style = "margin-top: 30px; border-bottom: 2px solid #00008B; padding-bottom: 5px; color: #00008B;")
      #), 
      if (isTRUE(active_content_type() == "metabolic_model")) {
        tags$div(id = "Metabolic_Information_Section",
                 tags$h3("Metabolic Information", class = "section-divider"))
      } else { NULL },
      
      if (isTRUE(active_content_type() == "metabolic_model")) {
        tagList(
          tags$div(id = "metabolic_model_box",
                   tags$h4("Genome-Scale Metabolic Model"),
                   uiOutput("gem_content")
          ),
          tags$hr()
        )
      } else { NULL }
    )
  })
  
  output$phylogeny_content_ui <- renderUI({
    req(current_bacterium(), isTRUE(active_content_type() == "phylogeny"))

    available_trees <- phylogeny_available_trees(current_bacterium()$filename)

    tree_box <- if (length(available_trees) == 0) {
      box(
        title = "Phylogenetic Tree - Interactive View",
        width = 12, solidHeader = TRUE, status = "warning",
        p(icon("exclamation-triangle"),
          " No phylogenetic tree has been built yet for ", tags$b(current_bacterium()$name), "."),
        p(style = "color:#666; font-size: 13px;",
          "Whole-genome and 16S rRNA gene trees (Type (Strain) Genome Server) are currently available ",
          "for the reference strain Acutalibacter muris KB18; trees for the other members will be added. ",
          "The relationships of all 12 members are shown on the ", tags$b("OMM12 Relatedness Tree"), " page.")
      )
    } else {
      box(
        title = "Phylogenetic Tree - Interactive View",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        fluidRow(
          column(12,
                 tags$h5("The tree builder service from Type (Strain) Genome Server, calculates both genome and 16S rRNA gene GBDP trees including branch support14. Phylogenies are inferred using FastME 2.1.4 with a BioNJ starting tree and Subtree Pruning and Regrafting postprocessing55,56; the search settings for original and bootstrapped distance matrices are identical."
                 )),
          column(
            width = 12,
            div(
              style = "padding: 10px; background-color: #e8edf3; border-left: 4px solid #c9a227; border-radius: 4px; margin-bottom: 15px;",
              tags$strong(icon("info-circle"), " Interactive Features:"),
              tags$ul(
                style = "margin-top: 8px; margin-bottom: 5px;",
                tags$li(tags$strong("Zoom:"), " Use mouse wheel to zoom in/out (0.5x - 10x)"),
                tags$li(tags$strong("Pan:"), " Click and drag to move around the tree"),
                tags$li(tags$strong("Tooltips:"), " Hover over nodes and labels to see detailed information"),
                tags$li(tags$strong("Highlight:"), " Click on internal nodes to highlight entire clades"),
                tags$li(tags$strong("Export:"), " Use the toolbar button to save as PNG")
              )
            )
          ),
          if (length(available_trees) > 1) {
            column(width = 6,
                   radioButtons("phylo_tree_kind", "Tree source",
                                choices = c("Whole-genome (GBDP) tree" = "genome",
                                            "16S rRNA gene tree" = "16s"),
                                selected = "genome", inline = TRUE))
          } else NULL,
          column(width = if (length(available_trees) > 1) 6 else 12,
                 radioButtons("phylo_branch_mode", "Branch style",
                              choices = c("Phylogram (branch length = genetic distance)" = "phylogram",
                                          "Cladogram (equal spacing, topology only)" = "cladogram"),
                              selected = "phylogram", inline = TRUE))
        ),
        girafeOutput("phylogeny_plot_enhanced", height = "800px")
      )
    }

    tagList(
      tree_box

      # Conditional similarity data table -- COMMENTED OUT for now.
      # Similarity Data & Strain Information feature paused; restore by
      # uncommenting this block (and re-adding the trailing comma after tree_box).
      # ,
      # if (!is.null(current_bacterium_similarity_data())) {
      #   box(
      #     title = HTML("<i class='fa fa-table'></i> Similarity Data & Strain Information"),
      #     width = 12,
      #     solidHeader = TRUE,
      #     status = "info",
      #     p("This table shows similarity metrics and additional information for each strain in the phylogenetic tree."),
      #     DTOutput("similarity_table")
      #   )
      # } else {
      #   box(
      #     title = "Similarity Data",
      #     width = 12,
      #     solidHeader = TRUE,
      #     status = "warning",
      #     p(icon("exclamation-triangle"),
      #       " No similarity data file found for this bacterium. To add similarity data:"),
      #     tags$ol(
      #       tags$li("Create a CSV file named: ",
      #               tags$code(paste0(current_bacterium()$filename, "_similarity.csv"))),
      #       tags$li("Include columns: strain_name, similarity_to_KB1, and any other metrics"),
      #       tags$li("Place the file in the 'data/' folder"),
      #       tags$li("Reload the page to see the data in tooltips and this table")
      #     )
      #   )
      # }
    )
  })
  
  output$phylogeny_plot_enhanced <- renderGirafe({
    req(current_bacterium())
    
    similarity_data <- current_bacterium_similarity_data()
    tree_kind <- if (!is.null(input$phylo_tree_kind)) input$phylo_tree_kind else "genome"
    branch_mode <- if (!is.null(input$phylo_branch_mode)) input$phylo_branch_mode else "phylogram"
    p <- get_phylogeny_plot_enhanced(current_bacterium()$filename, similarity_data,
                                      tree_kind = tree_kind, branch_mode = branch_mode)
    req(p)
    
    # Highlight selected clade if any
    if (!is.null(selected_clade())) {
      p <- p + geom_hilight(node = selected_clade(), fill = "gold", alpha = 0.3)
    }
    
    girafe(
      ggobj = p,
      width_svg = 14,
      height_svg = 10,
      options = list(
        opts_zoom(min = 0.5, max = 10),
        opts_hover(
          css = "cursor:pointer;fill:#c9a227;stroke:#c9a227;stroke-width:2px;",
          reactive = TRUE
        ),
        opts_hover_inv(css = "opacity:0.3;"),
        opts_tooltip(
          css = "background-color:#1f2937;color:#f9fafb;padding:12px;border-radius:6px;font-size:13px;box-shadow:0 4px 6px rgba(0,0,0,0.3);font-family:Arial, sans-serif;line-height:1.5;",
          opacity = 0.95,
          use_fill = FALSE,
          use_stroke = FALSE,
          delay_mouseover = 200,
          delay_mouseout = 500
        ),
        opts_selection(
          type = "single",
          css = "fill:#dc2626;stroke:#dc2626;stroke-width:3px;"
        ),
        opts_toolbar(
          position = "topright",
          saveaspng = TRUE,
          pngname = paste0("phylo_tree_", current_bacterium()$filename)
        ),
        opts_sizing(rescale = TRUE, width = 1)
      )
    )
  })
  
  output$similarity_table <- renderDT({
    similarity_data <- current_bacterium_similarity_data()
    req(similarity_data)
    
    datatable(
      similarity_data,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        searching = TRUE,
        ordering = TRUE,
        autoWidth = TRUE,
        columnDefs = list(
          list(className = 'dt-center', targets = '_all')
        )
      ),
      class = 'cell-border stripe hover',
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: left; color: #374151; font-size: 14px; font-weight: bold; padding: 10px;',
        'Strain Similarity and Metadata Information'
      ),
      filter = 'top'
    )
  })
  
  output$defense_content <- renderUI({
    req(active_content_type() == "defense")
    
    # Which activities actually exist in this strain's file (for an info note)
    sys <- current_defense_systems()
    avail <- if (!is.null(sys) && "activity" %in% colnames(sys)) {
      sort(unique(trimws(sys$activity)))
    } else character(0)
    
    box(
      title = paste("Defense Systems for", current_bacterium()$name),
      width = 12, solidHeader = TRUE, status = "warning",
      fluidRow(
        column(12,
               tags$h5("Defense systems predicted with DefenseFinder. Choose an activity class below: 'Defense' systems protect the host (e.g. restriction-modification, CBASS); 'Antidefense' systems counteract host defenses (e.g. anti-CRISPR).")
        )
      ),
      fluidRow(
        column(12,
               radioButtons(
                 "defense_activity",
                 label = "Activity class:",
                 choices = c("Defense", "Antidefense"),
                 selected = if ("Defense" %in% avail) "Defense"
                 else if (length(avail) > 0) avail[1]
                 else "Defense",
                 inline = TRUE
               )
        )
      ),
      fluidRow(
        column(width = 6,
               plotlyOutput("defense_donut_chart", height = "400px")
        ),
        column(width = 6,
               h4("Defense System Types"),
               DTOutput("defense_summary_table")
        )
      ),
      tags$hr(style = "border-top: 1px solid #ccc;"),
      fluidRow(
        column(width = 12,
               h4(HTML(paste0("Genes in System: ",
                              tags$span(id = "defense-selected-type",
                                        "Select a system from the chart...")))),
               DTOutput("defense_genes_table")
        )
      )
    )
  })
  
  output$defense_donut_chart <- renderPlotly({
    counts <- defense_type_counts()
    if (is.null(counts) || nrow(counts) == 0 || counts$Type[1] == "No data") {
      return(plotly_empty() %>%
               layout(
                 title = "",
                 annotations = list(
                   x = 0.5, y = 0.5,
                   text = paste0("No ", input$defense_activity, " systems found in this genome."),
                   xref = "paper", yref = "paper", showarrow = FALSE,
                   font = list(size = 16, color = "#1b263b")
                 )
               ))
    }
    plot_ly(counts, labels = ~Type, values = ~Count, type = "pie",
            hole = 0.6,
            marker = list(line = list(color = "#ffffff", width = 1)),
            textinfo = "label+value",
            insidetextorientation = "radial",
            hoverinfo = "text",
            text = ~paste(Type, ":", Count),
            source = "defense_donut_clicks") %>%
      layout(title = paste(input$defense_activity, "System Distribution"),
             showlegend = TRUE) %>%
      plotly::event_register("plotly_click")
  })
  
  output$defense_summary_table <- renderDT({
    counts <- defense_type_counts()
    if (is.null(counts) || nrow(counts) == 0 || counts$Type[1] == "No data") {
      return(datatable(
        data.frame(Message = paste0("No ", input$defense_activity, " systems available.")),
        options = list(dom = "t"), rownames = FALSE))
    }
    datatable(counts, options = list(dom = "t", pageLength = 20), rownames = FALSE,
              colnames = c("System Type", "Count"))
  })
  
  output$defense_genes_table <- renderDT({
    datatable(
      filtered_defense_data_reactive(),
      options = list(dom = "tip", pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$crispr_bact_content <- renderUI({
    req(active_content_type() == "crispr", current_bacterium())
    bf <- current_bacterium()$filename
    arrays <- load_crispr_arrays(bf)
    cas <- load_cas_systems(bf)

    tagList(
      box(
        title = paste("CRISPR arrays for", current_bacterium()$name),
        width = 12, solidHeader = TRUE, status = "warning",
        p(style = "color:#666; font-size: 12.5px;",
          "See the community-wide ", tags$b("CRISPR & Cas Systems"), " page for the evidence-level ",
          "and orphan-array method note."),
        if (!is.null(arrays) && nrow(arrays) > 0) {
          DTOutput("crispr_bact_arrays_table")
        } else if (!is.null(arrays)) {
          p(style = "color:#1b7a3d; font-size: 13px;",
            icon("circle-check"), " CRISPRCasFinder ran for this bacterium and found no arrays.")
        } else {
          p(style = "color:#666; font-size: 13px;",
            "CRISPRCasFinder results are not yet available for this strain.")
        }
      ),
      box(
        title = paste("Cas gene clusters for", current_bacterium()$name),
        width = 12, solidHeader = TRUE, status = "warning",
        if (!is.null(cas) && nrow(cas) > 0) {
          DTOutput("crispr_bact_cas_table")
        } else if (!is.null(cas)) {
          p(style = "color:#1b7a3d; font-size: 13px;",
            icon("circle-check"), " CRISPRCasFinder ran for this bacterium and found no Cas gene cluster ",
            "(any CRISPR array above is an orphan array).")
        } else {
          p(style = "color:#666; font-size: 13px;",
            "CRISPRCasFinder results are not yet available for this strain.")
        }
      )
    )
  })

  output$crispr_bact_arrays_table <- renderDT({
    req(current_bacterium())
    d <- load_crispr_arrays(current_bacterium()$filename)
    req(d, nrow(d) > 0)
    d <- d %>%
      select(Contig = Sequence, Start, End, `Length (bp)` = Length, Orientation,
             `Repeat consensus` = DR_Consensus, Spacers = Spacers_Nb,
             `Evidence level (1-4)` = Evidence_Level,
             `Repeat conservation %` = Conservation_DRs_pct) %>%
      arrange(Start)
    datatable(d, rownames = FALSE, selection = "none",
              options = list(pageLength = 10, scrollX = TRUE))
  })

  output$crispr_bact_cas_table <- renderDT({
    req(current_bacterium())
    d <- load_cas_systems(current_bacterium()$filename)
    req(d, nrow(d) > 0)
    d <- d %>%
      select(Contig = Sequence, `Cas type` = Cas_Cluster_Type,
             `Cluster start` = Cluster_Start, `Cluster end` = Cluster_End,
             Gene = Gene_Subtype, `Gene start` = Gene_Start, `Gene end` = Gene_End,
             Strand = Gene_Orientation) %>%
      arrange(`Cluster start`, `Gene start`)
    datatable(d, rownames = FALSE, selection = "none",
              options = list(pageLength = 10, scrollX = TRUE))
  })

  output$bgc_bact_content <- renderUI({
    req(active_content_type() == "bgc", current_bacterium())
    bf <- current_bacterium()$filename
    f <- file.path("data", paste0(bf, "_antismash_regions.tsv"))
    regions <- if (file.exists(f)) {
      tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
    } else NULL

    tagList(
      box(
        title = paste("Predicted BGC regions for", current_bacterium()$name),
        width = 12, solidHeader = TRUE, status = "warning",
        p(style = "color:#666; font-size: 12.5px;",
          "See the community-wide ", tags$b("Biosynthetic Gene Clusters"), " page for the antiSMASH/",
          "MIBiG method note."),
        if (!is.null(regions) && nrow(regions) > 0) {
          tagList(p(style = "color:#666; font-size: 12.5px;",
                    "Click a row to see that region's gene-arrow diagram below."),
                  DTOutput("bgc_bact_table"))
        } else if (!is.null(regions)) {
          p(style = "color:#1b7a3d; font-size: 13px;",
            icon("circle-check"), " antiSMASH ran for this bacterium and found no BGC regions.")
        } else {
          p(style = "color:#666; font-size: 13px;",
            "antiSMASH results are not yet available for this strain.")
        }
      ),
      if (!is.null(regions) && nrow(regions) > 0) {
        box(
          title = "Gene-arrow diagram for selected region",
          width = 12, solidHeader = TRUE, status = "warning",
          p(style = "color:#666; font-size: 12px;",
            "Colored by antiSMASH's own gene classification: dark red = core biosynthetic genes, ",
            "pink = additional biosynthetic genes, blue = transport-related, green = regulatory, ",
            "grey = other/uninvolved genes in the region."),
          shinycssloaders::withSpinner(plotOutput("bgc_region_plot", height = "260px"),
                                        type = 6, color = "#1b263b")
        )
      } else NULL
    )
  })

  output$bgc_bact_table <- renderDT({
    req(current_bacterium())
    f <- file.path("data", paste0(current_bacterium()$filename, "_antismash_regions.tsv"))
    req(file.exists(f))
    d <- tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
    req(d, nrow(d) > 0)
    d <- d %>%
      select(Region, Contig, Start, End, `Length (bp)` = Length_bp,
             `Predicted product` = Predicted_Product, Category,
             `Known cluster (MIBiG)` = Known_Cluster_Accession,
             `Known cluster description` = Known_Cluster_Description,
             `Similarity %` = Known_Cluster_Similarity_pct) %>%
      arrange(Start)
    datatable(d, rownames = FALSE, selection = "single",
              options = list(pageLength = 10, scrollX = TRUE))
  })

  output$bgc_region_plot <- renderPlot({
    req(current_bacterium())
    f <- file.path("data", paste0(current_bacterium()$filename, "_antismash_regions.tsv"))
    req(file.exists(f))
    d <- tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
    req(d, nrow(d) > 0)
    d <- d %>% arrange(Start)
    row_idx <- input$bgc_bact_table_rows_selected
    region_id <- if (!is.null(row_idx) && length(row_idx) > 0) d$Region[row_idx[1]] else d$Region[1]
    p <- render_bgc_region_plot(current_bacterium()$filename, region_id)
    req(p)
    p
  })

  output$specialty_genes_bact_content <- renderUI({
    req(active_content_type() == "specialty_genes", current_bacterium())
    bf <- current_bacterium()$filename
    hits <- if (nrow(specialty_genes_data) > 0) {
      specialty_genes_data %>% filter(Filename == bf) %>%
        select(Category, Subcategory, Gene, Product, `Locus Tag` = Locus_Tag) %>%
        arrange(Category, Subcategory)
    } else data.frame()
    pheno <- load_amr_phenotype(bf)
    amrf <- load_amrfinder(bf)

    tagList(
      box(
        title = paste("Curated AMR calls (AMRFinderPlus) for", current_bacterium()$name),
        width = 12, solidHeader = TRUE, status = "success",
        if (!is.null(amrf) && nrow(amrf) > 0) {
          DTOutput("specialty_genes_bact_amrfinder_table")
        } else if (!is.null(amrf)) {
          p(style = "color:#1b7a3d; font-size: 13px;",
            icon("circle-check"), " AMRFinderPlus ran for this bacterium and found zero AMR elements.")
        } else {
          p(style = "color:#666; font-size: 13px;",
            "AMRFinderPlus results are not yet available for this strain; the keyword screen below is shown instead.")
        }
      ),
      box(
        title = paste("Keyword screen for", current_bacterium()$name),
        width = 12, solidHeader = TRUE, status = "warning",
        p(style = "color:#666; font-size: 12.5px;",
          "Keyword screen against Prokka's own gene/product annotations for antimicrobial-resistance ",
          "and virulence-associated gene families — not a curated CARD/ResFinder/VFDB alignment call. ",
          "See the community-wide ", tags$b("AMR & Specialty Genes"), " page for the full method note."),
        if (nrow(hits) == 0) {
          p(style = "color:#666;", "No keyword hits for this bacterium.")
        } else {
          DTOutput("specialty_genes_bact_table")
        }
      ),
      box(
        title = "AMR phenotype (lab-tested)",
        width = 12, solidHeader = TRUE, status = "warning",
        if (is.null(pheno) || nrow(pheno) == 0) {
          tagList(
            p(icon("exclamation-triangle"), " No lab-tested AMR phenotype data on file for this bacterium. To add it:"),
            tags$ol(
              tags$li("Create a TSV file named: ", tags$code(paste0(bf, "_amr_phenotype.tsv"))),
              tags$li("Include columns such as: antibiotic, method (MIC/disk-diffusion), result (S/I/R), mic_value"),
              tags$li("Place the file in the 'data/' folder and reload the page.")
            )
          )
        } else {
          DTOutput("specialty_genes_bact_pheno_table")
        }
      )
    )
  })

  output$specialty_genes_bact_table <- renderDT({
    req(current_bacterium())
    bf <- current_bacterium()$filename
    hits <- if (nrow(specialty_genes_data) > 0) {
      specialty_genes_data %>% filter(Filename == bf) %>%
        select(Category, Subcategory, Gene, Product, `Locus Tag` = Locus_Tag) %>%
        arrange(Category, Subcategory)
    } else data.frame()
    datatable(hits, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
  })

  output$specialty_genes_bact_amrfinder_table <- renderDT({
    req(current_bacterium())
    amrf <- load_amrfinder(current_bacterium()$filename)
    req(amrf)
    d <- amrf %>%
      select(Gene = `Element symbol`, Name = `Element name`, Type, Class, Subclass,
             `% Identity` = `% Identity to reference`, `% Coverage` = `% Coverage of reference`,
             Method, Contig = `Contig id`, Start, Stop) %>%
      arrange(Class)
    datatable(d, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
  })

  output$specialty_genes_bact_pheno_table <- renderDT({
    req(current_bacterium())
    pheno <- load_amr_phenotype(current_bacterium()$filename)
    req(pheno)
    datatable(pheno, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
  })


  output$gem_table <- renderDT({
    df <- gem_filtered()
    validate(need(!is.null(df) && nrow(df) > 0, "No reactions in this pathway."))
    datatable(
      df, rownames = FALSE, extensions = "Buttons",
      options = list(
        pageLength = 10, scrollX = TRUE,
        dom = "Bfrtip", buttons = c("copy", "csv", "excel"),
        columnDefs = list(list(width = "35%", targets = 2))
      )
    )
  })
  
  current_bacterium_genes <- reactive({
    req(active_content_type() == "genome", current_bacterium())
    bacterium_filename <- current_bacterium()$filename
    
    gff_data <- load_prokka_gff(bacterium_filename)
    
    if (!is.null(gff_data) && nrow(gff_data) > 0) {
      genes_df <- gff_data %>%
        filter(!is.na(type), type %in% c("CDS", "tRNA", "rRNA", "tmRNA", "misc_RNA")) %>%
        filter(!is.na(start), !is.na(end), !is.na(strand)) %>%
        mutate(
          gene = ifelse(is.na(gene) | gene == "", locus_tag, gene),
          gene = ifelse(is.na(gene) | gene == "", paste0(type, "_", row_number()), gene),
          description = ifelse(is.na(product) | product == "", type, product),
          strand = ifelse(strand == "-", "-", "+"),
          COG_category = ifelse(is.na(COG_category) | COG_category == "", "Unknown", COG_category)
        ) %>%
        select(gene, start, end, strand, description, type,COG_category)
      
      if (nrow(genes_df) > 0) {
        return(genes_df)
      }
    }
    
    file_path <- file.path("data", paste0(bacterium_filename, "_gene_df.tsv"))
    if (file.exists(file_path)) {
      genes_df <- read_tsv(file_path, show_col_types = FALSE) %>%
        filter(!is.na(start), !is.na(end), !is.na(strand)) %>%
        mutate(
          strand = ifelse(tolower(strand) %in% c("reverse", "-"), "-", "+"),
          type = "CDS"
        )
      
      if (!"type" %in% colnames(genes_df)) {
        genes_df <- genes_df %>% mutate(type = "CDS")
      }
      
      return(genes_df)
    } else {
      return(NULL)
    }
  })
  
  output$genome_plot_ui <- renderUI({
    req(active_content_type() == "genome")
    bacterium <- current_bacterium()
    genes_df <- current_bacterium_genes()
    
    if (is.null(genes_df) || nrow(genes_df) == 0) {
      return(box(width = 12, status = "warning", solidHeader = TRUE,
                 paste("No gene data found for", bacterium$name)))
    }
    
    current_genome_min <- min(genes_df$start)
    current_genome_max <- max(genes_df$end)
    
    box(
      title = "Genome Viewer (Linear / Circular)", width = 12, solidHeader = TRUE, status = "info",

      # --- Control toolbar: region, search, color-by, all in one aligned panel ---
      div(
        style = paste(
          "background-color:#f4f6f9; border:1px solid #e2e8f0; border-radius:6px;",
          "padding:16px 18px; margin-bottom:16px;"
        ),
        fluidRow(
          column(width = 7,
                 sliderInput(paste0("region_", bacterium$filename),
                             "Genomic region",
                             min = current_genome_min, max = current_genome_max,
                             value = c(current_genome_min, min(current_genome_min + 50000, current_genome_max)),
                             step = 1000, width = "100%")
          ),
          column(width = 3,
                 tags$label("Search gene", style = "font-weight: 600; display:block; margin-bottom: 6px;"),
                 uiOutput("genome_search_ui")
          ),
          column(width = 2,
                 radioButtons("gene_map_color_choice",
                              "Color genes by",
                              choices = c("Default (Yellow)", "COG Category"),
                              selected = "Default (Yellow)"
                 )
          )
        )
      ),

      tabsetPanel(
        tabPanel("Linear", plotlyOutput("gene_plot", height = "250px"),
                 fluidRow(
                   style = "margin-top: 12px;",
                   column(3, downloadButton("download_prokka_tsv", "Download Prokka TSV", class = "btn-sm")),
                   column(3, downloadButton("download_prokka_gff", "Download Prokka GFF", class = "btn-sm")),
                   column(3, downloadButton("download_prokka_gbk", "Download Prokka GBK", class = "btn-sm")),
                   column(3, downloadButton("download_cog_gff", "Download COG table (TSV)", class = "btn-sm"))
                 )
        ),
        tabPanel("Circular",
                 p(style = "color:#666; font-size: 12px; margin-top: 12px;",
                   "Interactive genome map (CGView.js) — scroll to zoom, drag to pan, click features for details."),
                 tags$div(
                   id = "cgview_container",
                   style = "width: 100%; height: 600px; border: 1px solid #e2e8f0; border-radius: 6px; margin-top: 8px;"
                 )
        )
      )
    )
  })

  # (Re)draw the interactive circular genome map whenever the genome tab is
  # opened or the selected bacterium changes.
  observeEvent(list(active_content_type(), current_bacterium()), {
    req(active_content_type() == "genome", current_bacterium())
    gbk_url <- paste0("data/", current_bacterium()$filename, "_prokka.gbk")
    session$sendCustomMessage("loadCGView", list(gbkUrl = gbk_url))
  })
  
  output$gene_plot <- renderPlotly({
    req(active_content_type() == "genome", current_bacterium())
    bacterium <- current_bacterium()
    slider_input_id <- paste0("region_", bacterium$filename)
    filter_input_id <- paste0("feature_filter_", bacterium$filename)
    searched_gene <- input$search_gene

    req(input[[slider_input_id]])

    genes_df <- current_bacterium_genes()
    req(!is.null(genes_df))
    
    if (!is.null(input[[filter_input_id]]) && length(input[[filter_input_id]]) > 0) {
      genes_df <- genes_df %>% filter(type %in% input[[filter_input_id]])
    }
    
    if (nrow(genes_df) == 0) {
      return(ggplotly(ggplot() +
                        annotate(x = 0.5, y = 0.5,
                                 label = "No features selected",
                                 size = 5) +
                        theme_void()))
    }
    
    region_start <- max(1, input[[slider_input_id]][1])
    region_end <- min(max(genes_df$end), input[[slider_input_id]][2])
    
    filtered_genes <- genes_df %>%
      filter(start <= region_end, end >= region_start)
    
    if (nrow(filtered_genes) == 0) {
      return(ggplotly(ggplot() +
                        annotate(x = 0.5, y = 0.5,
                                 label = "No genes in selected region.") +
                        theme_void()))
    }
    
    wrap_width <- 50000
    filtered_genes <- filtered_genes %>%
      mutate(row = floor((start - min(start)) / wrap_width), y = 0)
    
    arrow_df <- filtered_genes %>%
      rowwise() %>%
      do({
        gene <- .
        len <- gene$end - gene$start
        arrow_head <- min(300, len * 0.3)
        body_len <- len - arrow_head
        x0 <- gene$start
        x1 <- x0 + body_len
        x2 <- gene$end
        
        if (gene$strand == "+") {
          xs <- c(x0, x1, x1, x2, x1, x1, x0)
        } else {
          xs <- c(x2, x1, x1, x0, x1, x1, x2)
        }
        
        data.frame(
          gene = gene$gene,
          row = gene$row,
          x = xs,
          y = c(-0.2, -0.2, -0.4, 0, 0.4, 0.2, 0.2),
          strand = gene$strand,
          description = gene$description,
          type = gene$type,
          COG_category = gene$COG_category,
          text = paste0(
            "<b>Gene:</b> ", gene$gene,
            "<br><b>Type:</b> ", gene$type,
            "<br><b>Product:</b> ", gene$description,
            "<br><b>Start:</b> ", gene$start,
            "<br><b>End:</b> ", gene$end,
            "<br><b>Strand:</b> ", gene$strand,
            "<br><b>COG:</b> ", gene$COG_category
          )
        )
      }) %>%
      ungroup()
    
    # Define base colors for non-CDS features
    base_type_colors <- c(
      "tRNA" = "#FD8D3C",
      "rRNA" = "#74C476",
      "tmRNA" = "#9E9AC8",
      "misc_RNA" = "#F768A1"
    )
    
    # 🌟 CONDITIONAL COLORING LOGIC 🌟
    if (input$gene_map_color_choice == "COG Category") {
      
      # Define Group colors (used as values for the scale)
      group_colors_map <- c(
        "INFORMATION STORAGE AND PROCESSING" = "#4C72B0",
        "CELLULAR PROCESSES AND SIGNALING" 	= "#DD8452",
        "METABOLISM" 	= "#55A868",
        "POORLY CHARACTERIZED" 	= "#8172B2",
        "UNKNOWN"="grey"
      )
      
      # 🟢 CRITICAL STEP: Join the COG Group name to the data
      legend_df <- cog_legend()     # Assumed table with Group and Letter
      
      arrow_df <- arrow_df %>%
        # Join the COG Group name (Group) based on the COG Letter (COG_category)
        left_join(legend_df %>% select(Letter, Group), 
                  by = c("COG_category" = "Letter")) %>%
        
        # If CDS but Group is NA (e.g., 'None' or 'Unknown' COG), assign 'Unknown COG'
        mutate(COG_Group = ifelse(is.na(Group) & type == "CDS", "UNKNOWN", Group)) %>%
        
        # Define fill_group: Use COG Group Name for CDS, use Type Name for others
        mutate(fill_group = ifelse(type == "CDS", COG_Group, type))
      
      # 🟢 Define the combined color map (Group Name/Type Name -> Color)
      combined_colors <- c(group_colors_map, 
                           "Unknown" = "grey", 
                           base_type_colors)
      
      fill_colors <- combined_colors
      fill_name <- "COG Group / Feature Type"
      
      # Since fill_group now contains the correct label (Group Name), the labels argument is simpler
      final_labels <- names(fill_colors) 
      
    } else {
      # 2. Default Coloring Setup (Yellow for CDS, base colors for others)
      default_colors <- c("CDS" = "yellow", base_type_colors)
      
      # Ensure fill_group is defined for the else case too
      arrow_df <- arrow_df %>%
        mutate(fill_group = type)
      
      fill_colors <- default_colors
      fill_name <- "Feature Type"
      final_labels <- names(default_colors)
    }
    
    if (is.null(searched_gene) || searched_gene == "") {
      searched_gene <- NA
    }
    
    arrow_df <- arrow_df %>%
      mutate(
        highlight = case_when(
          !is.na(searched_gene) &
            (
              contains_ci(gene, searched_gene) |
                contains_ci(description, searched_gene)
            ) ~ TRUE,
          TRUE ~ FALSE 
        )
      ) %>%
      mutate(
        display_fill = ifelse(highlight, "HIGHLIGHT_GENE", fill_group)
      )
    
    
    line_df <- filtered_genes %>%
      group_by(row) %>%
      summarise(x_start = min(start), x_end = max(end), .groups = "drop") %>%
      mutate(y = 0)
    
    highlight_color <- "#FF0000"  # bright red
    
    fill_colors <- c(fill_colors, "HIGHLIGHT_GENE" = highlight_color)
    
    final_labels <- c(final_labels, "HIGHLIGHT_GENE" = "Highlighted Gene")
    
    p <- ggplot() +
      geom_segment(data = line_df,
                   aes(x = x_start, xend = x_end, y = y, yend = y),
                   color = "black", linewidth = 0.2) +
      geom_polygon(
        data = arrow_df,
        aes(x = x, y = y, group = gene, fill = display_fill, text = text),
        color = ifelse(arrow_df$highlight, "black", "black"),
        linewidth = ifelse(arrow_df$highlight, 1.2, 0.3)
      ) +
      # Final colors and labels are passed
      scale_fill_manual(values = fill_colors, name = fill_name, labels = final_labels) +
      facet_wrap(~ row, scales = "free_x", ncol = 1) +
      theme_minimal() +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom"
      )
    
    ggplotly(p, tooltip = "text") %>%
      config(scrollZoom = TRUE) %>%
      layout(legend = list(orientation = "h", x = 0.5, xanchor = "center"))
  })
  
  output$genome_search_ui <- renderUI({
    req(active_content_type() == "genome")
    bacterium <- current_bacterium()
    req(bacterium)
    genes_df <- current_bacterium_genes()
    gene_list <- sort(unique(c(genes_df$gene, genes_df$description)))
    
    selectizeInput(
      inputId = "search_gene",
      label = NULL,
      # leading "" = empty default: start with no gene selected unless one was
      # passed in (e.g. from Cross-Genome Search); otherwise selectize would
      # auto-select the first gene and jump the map to it
      choices = c("", gene_list),
      selected = if (is.null(pending_gene_search())) "" else pending_gene_search(),
      multiple = FALSE,
      options = list(
        placeholder = "Search gene…",
        highlight = TRUE,
        maxOptions = 50
      ),
      width = "100%"
    )
  })
  
  output$annotation_content <- renderUI({
    req(current_bacterium(), isTRUE(active_content_type() == "annotation"))
    
    tagList(
      fluidRow(
        box(
          title = "Feature Annotation Summary", width = 5, solidHeader = TRUE, status = "primary",
          plotlyOutput("annotation_pie_chart", height = "300px")
        ),
        box(
          title = HTML("Annotation Stats <span id='annotation-selected-type' style='color: #FF8C00; font-size: 14px; margin-left: 10px;'>Select a type from the chart...</span>"), 
          width = 7, solidHeader = TRUE, status = "info",
          DTOutput("annotation_details_table")
        )
      ),
      tags$hr(),
      box(
        title = "Prokka Annotation Raw Summary", width = 12, solidHeader = TRUE, status = "default",
        p("This table displays the summary statistics from the Prokka run (tsv file)."),
        DTOutput("prokka_summary_table")
      )
    )
  })
  
  output$annotation_pie_chart <- renderPlotly({
    data <- annotation_stats()
    req(data, nrow(data) > 0, all(c("type", "Count") %in% colnames(data)))
    
    if (data$type[1] == "No data") {
      return(
        ggplotly(
          ggplot() + annotate(x = 0, y = 0, label = "No GFF data available.") + theme_void()
        )
      )
    }
    
    plot_ly(
      data, 
      labels = ~type, 
      values = ~Count, 
      type = 'pie', 
      source = "annotation_pie_chart_clicks",
      marker = list(colors = RColorBrewer::brewer.pal(min(nrow(data), 9), "Set1"),
                    line = list(color = '#FFFFFF', width = 1)),
      textinfo = 'percent',
      hoverinfo = 'text',
      text = ~paste('Type:', type, '<br>Count:', Count)
    ) %>%
      layout(title = 'Annotated Feature Types',
             showlegend = TRUE,
             margin = list(t = 50, b = 50))
  })
  
  output$annotation_details_table <- renderDT({
    data <- filtered_annotation_data_reactive()
    req(data)
    
    if ("Message" %in% colnames(data)) {
      return(datatable(data, options = list(dom = 't', ordering = FALSE, pageLength = 1)))
    }
    
    datatable(
      data,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        searching = TRUE
      )
    )
  })
  
  output$prokka_summary_table <- renderDT({
    data <- current_prokka_tsv()
    req(data, nrow(data) > 0)
    
    datatable(
      data,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        searching = TRUE
      )
    )
  })
  
  
  output$cog_content <- renderUI({
    req(active_content_type() == "cog")
    eggnog <- load_eggnog_annotation(current_bacterium()$filename)
    tagList(
      box(
        title = paste("COG Category Distribution for", current_bacterium()$name),
        width = 12, solidHeader = TRUE, status = "info",
        fluidRow(
          column(width = 12,
                 plotlyOutput("cog_barplot", height = "500px")
          )
        ),
        tags$hr(style = "border-top: 1px solid #ccc;"),
        fluidRow(
          column(
            width = 12,
            h4(HTML(paste0("Genes in COG Category: ", tags$span(id="cog-selected-type", "(Select a COG category from the bar chart)")))),
            DTOutput("cog_genes_table")
          )
        )
      ),
      # eggNOG-mapper only run so far for the reference strain -- this box
      # simply doesn't render for the other 11 bacteria (eggnog is NULL).
      if (!is.null(eggnog) && nrow(eggnog) > 0) {
        box(
          title = paste("eggNOG-mapper Functional Annotation for", current_bacterium()$name,
                         "(reference strain only)"),
          width = 12, solidHeader = TRUE, status = "primary", collapsible = TRUE, collapsed = TRUE,
          p(style = "color:#666; font-size: 12.5px;",
            "Supplemental annotation from ", tags$a(href = "https://github.com/eggnogdb/eggnog-mapper",
                                                       target = "_blank", "eggNOG-mapper"),
            " — orthologous groups, Pfam domains, and a preferred gene name, alongside the COG category ",
            "already used above. Currently only run for this reference strain; richer/complementary to, ",
            "not a replacement for, Prokka's own annotation used everywhere else in this app."),
          DTOutput("eggnog_table")
        )
      } else NULL
    )
  })

  output$eggnog_table <- renderDT({
    req(current_bacterium())
    d <- load_eggnog_annotation(current_bacterium()$filename)
    req(d, nrow(d) > 0)
    d <- d %>%
      select(`Locus Tag` = Locus_Tag, `COG category` = COG_category, Description,
             `Preferred name` = Preferred_Name, `Orthologous groups` = OGs, Pfam = PFAMs,
             Score, `E-value` = Evalue) %>%
      arrange(`Locus Tag`)
    datatable(d, rownames = FALSE, selection = "none",
              options = list(pageLength = 10, scrollX = TRUE))
  })
  output$kegg_content <- renderUI({
    req(active_content_type() == "kegg")
    df <- current_bacterium_kegg_data()
    
    if (is.null(df) || nrow(df) == 0) {
      return(box(width = 12, status = "warning", solidHeader = TRUE,
                 "KEGG pathway annotation is not yet available for this strain."))
    }
    
    # keep only complete hierarchy rows, trim stray whitespace/CR
    df$Main_Pathway <- trimws(df$Main_Pathway)
    df$Sub_pathway  <- trimws(df$Sub_pathway)
    df$Pathway      <- trimws(df$Pathway)
    df <- df[complete.cases(df[, c("Main_Pathway","Sub_pathway","Pathway")]), ]
    
    main_levels <- sort(unique(df$Main_Pathway))
    
    # Build one collapsible box per Main_Pathway
    main_panels <- lapply(seq_along(main_levels), function(mi) {
      main <- main_levels[mi]
      df_main <- df[df$Main_Pathway == main, ]
      sub_levels <- sort(unique(df_main$Sub_pathway))
      
      # Inside each main, a collapsible sub-box per Sub_pathway
      sub_panels <- lapply(seq_along(sub_levels), function(si) {
        sub <- sub_levels[si]
        df_sub <- df_main[df_main$Sub_pathway == sub, ]
        n_genes <- length(unique(df_sub$locus_tag))
        
        # unique, safe output id for this sub-table
        tbl_id <- paste0("kegg_tbl_", mi, "_", si)
        
        box(
          title = tags$span(
            style = "color:#1b263b; font-weight:600;",
            sub,
            tags$span(
              style = paste(
                "margin-left:8px; padding:2px 8px; border-radius:10px;",
                "background-color:#e8edf3; color:#2c3e57;",
                "font-weight:500; font-size:0.85em;"
              ),
              paste0(n_genes, " genes")
            )
          ),
          width = 12, collapsible = TRUE, collapsed = TRUE,
          solidHeader = FALSE, status = "primary",
          DTOutput(tbl_id)
        )
      })

      box(
        title = tags$span(
          style = "color:#ffffff; font-weight:600;",
          main,
          tags$span(
            style = paste(
              "margin-left:8px; padding:2px 10px; border-radius:10px;",
              "background-color:#c9a227; color:#0d1b2a;",
              "font-weight:600; font-size:0.85em;"
            ),
            paste0(length(sub_levels), " sub-pathways, ",
                   length(unique(df_main$locus_tag)), " genes")
          )
        ),
        width = 12, collapsible = TRUE, collapsed = TRUE,
        solidHeader = TRUE, status = "info",
        sub_panels
      )
    })
    
    tagList(
      box(
        title = paste("KEGG Pathways for", current_bacterium()$name),
        width = 12, solidHeader = TRUE, status = "info",
        p("Expand a main pathway, then a sub-pathway, to see the genes involved."),
        main_panels
      )
    )
  })
  # Render all KEGG accordion sub-tables
  
  # Render all KEGG accordion sub-tables
  observe({
    req(active_content_type() == "kegg")
    df <- current_bacterium_kegg_data()
    req(df)
    
    df$Main_Pathway <- trimws(df$Main_Pathway)
    df$Sub_pathway  <- trimws(df$Sub_pathway)
    df$Pathway      <- trimws(df$Pathway)
    df <- df[complete.cases(df[, c("Main_Pathway","Sub_pathway","Pathway")]), ]
    
    main_levels <- sort(unique(df$Main_Pathway))
    
    for (mi in seq_along(main_levels)) {
      df_main <- df[df$Main_Pathway == main_levels[mi], ]
      sub_levels <- sort(unique(df_main$Sub_pathway))
      
      for (si in seq_along(sub_levels)) {
        # capture the SUBSET DATA itself, not just indices
        df_sub <- df_main[df_main$Sub_pathway == sub_levels[si], ]
        tbl_id <- paste0("kegg_tbl_", mi, "_", si)
        
        local({
          df_local  <- df_sub      # freeze this iteration's data
          id_local  <- tbl_id
          
          output[[id_local]] <- renderDT({
            df_disp <- df_local %>%
              transmute(
                `Locus Tag`   = locus_tag,
                `Gene`        = ifelse(is.na(Gene_Symbol) | Gene_Symbol == "",
                                       locus_tag, Gene_Symbol),
                `Description` = Description,
                `KO`          = KO,
                `KO Module`   = KO_module,
                `Module`      = Pathway_1,
                `Pathway`     = Pathway
              ) %>%
              distinct()
            
            datatable(
              df_disp,
              rownames = FALSE,
              options = list(pageLength = 5, scrollX = TRUE, dom = "tip")
            )
          })
        })
      }
    }
  })
  
  output$kegg_table <- renderDT({
    df <- kegg_filter()
    validate(need(!is.null(df) && nrow(df) > 0,
                  "Click a node in the Sankey diagram to see the genes in that pathway."))
    
    datatable(
      df,
      rownames = FALSE,
      extensions = "Buttons",
      caption = htmltools::tags$caption(
        style = "caption-side: top; text-align: left; font-weight: bold;",
        paste0("Genes in: ", selected_kegg_node())
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel")
      )
    )
  })
  output$cog_barplot <- renderPlotly({
    counts <- cog_category_counts()
    if (is.null(counts) || nrow(counts) == 0 || sum(counts$Count) == 0) {
      return(plotly_empty() %>%
               layout(
                 title = 'No COG Data Available',
                 annotations = list(x = 0.5, y = 0.5, text = "COG TSV file not loaded or is empty.",
                                    xref = "paper", yref = "paper", showarrow = FALSE, font = list(size = 20, color = "red"))
               ) %>% plotly::event_register("plotly_click")
      )
    }
    
    present_categories <- as.character(counts$COG_category)
    counts$COG_category <- factor(counts$COG_category, levels = present_categories, exclude = NULL)
    full_colors <- cog_colors()
    filtered_colors <- full_colors[present_categories]
    counts$Color_Hex <- filtered_colors[counts$COG_category]
    
    p <- plot_ly(
      data = counts, x = ~COG_category, y = ~Count, type = 'bar',
      color = I(counts$Color_Hex),
      hovertemplate = paste("Category: %{x} (%{text})<br>", "Count: %{y}<extra></extra>"),
      text = ~Description,
      source = "cog_barplot_clicks"
    ) %>%
      layout(title = 'COG Category Counts', xaxis = list(title = "COG Category"),
             yaxis = list(title = "Number of Genes"), showlegend = FALSE) %>%
      plotly::event_register("plotly_click")
    return(p)
  })
  
  output$cog_genes_table <- renderDT({
    datatable(
      filtered_cog_data_reactive(),
      options = list(dom = 'tip', pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  output$localization_content <- renderUI({
    req(active_content_type() == "localization")
    loc <- load_localization_csv(current_bacterium()$filename)
    if (is.null(loc) || nrow(loc) == 0) {
      return(box(
        title = paste("Subcellular Localization Predictions for", current_bacterium()$name),
        width = 12, solidHeader = TRUE, status = "warning",
        p(icon("circle-info"), " Subcellular localization (DeepLocPro) has not yet been predicted for this strain.")
      ))
    }
    box(
      title = paste("Subcellular Localization Predictions for", current_bacterium()$name),
      width = 12, solidHeader = TRUE, status = "primary",
      fluidRow(
        column(12,
               tags$h5("DeepLocPro predicts the subcellular localization of prokaryotic proteins. It can differentiate between 6 different localizations: Cytoplasm, Cytoplasmic membrane, Periplasm, Outer membrane, Cell wall and surface, and extracellular space"                       )
        ),
        column(width = 12,
               h4("Localization Distribution"),
               plotlyOutput("localization_barplot", height = "400px")
        )
        # Cell Compartment Overview (schematic diagram) — disabled for now, not
        # informative enough as-is. Server-side code (build_localization_cell_diagram(),
        # output$localization_cell_diagram) is left in place to revisit later.
        # column(width = 6,
        #        h4("Cell Compartment Overview"),
        #        p(style = "color:#666; font-size: 12px;",
        #          "Schematic overview only — outer membrane/periplasm and cell wall/surface aren't ",
        #          "both present in the same cell; which applies depends on this bacterium's Gram status."),
        #        uiOutput("localization_cell_diagram")
        # )
      ),
      tags$hr(style = "border-top: 1px solid #ccc;"),
      fluidRow(
        column(
          width = 12,
          h4("Localization by Functional Category (COG)"),
          p(style = "color:#666; font-size: 13px;",
            "Which COG functional groups dominate each predicted compartment — darker cells mean more proteins."),
          shinycssloaders::withSpinner(plotlyOutput("localization_cog_heatmap", height = "400px"),
                                        type = 6, color = "#1b263b")
        )
      ),
      tags$hr(style = "border-top: 1px solid #ccc;"),
      fluidRow(
        column(
          width = 12,
          h4("Protein Projection by Predicted Localization"),
          p(style = "color:#666; font-size: 13px;",
            "Each point is one protein, projected from its 6-class DeepLocPro probability scores onto 2 components. ",
            "Note: built from the model's output probabilities, not sequence embeddings, so clusters mostly restate ",
            "the predicted class rather than validating it independently."),
          shinycssloaders::withSpinner(plotlyOutput("localization_pca_plot", height = "450px"),
                                        type = 6, color = "#1b263b")
        )
      ),
      tags$hr(style = "border-top: 1px solid #ccc;"),
      fluidRow(
        column(
          width = 12,
          h4(HTML(paste0("Genes at Location: ", tags$span(id="localization-selected-type", "(Select a location from the chart)")))),
          DTOutput("localization_genes_table")
        )
      )
    )
  })
  
  output$localization_barplot <- renderPlotly({
    counts <- localization_word_counts()
    if (nrow(counts) <= 1 && counts$word[1] == "No Data") {
      return(plotly_empty() %>% layout(title = 'No Localization Data Available'))
    }
    
    # Using a bar chart as an interactive proxy for a treemap/wordcloud
    p <- plot_ly(
      data = counts, type = "bar", x = ~word, y = ~freq, color = ~word,
      source = "localization_treemap_clicks"
    ) %>%
      layout(title = "Subcellular Localization Counts", xaxis = list(title = "Localization Category"),
             yaxis = list(title = "Count"), showlegend = FALSE) %>%
      plotly::event_register("plotly_click")
    return(p)
  })
  
  output$localization_genes_table <- renderDT({
    datatable(
      filtered_localization_data_reactive(),
      options = list(dom = 'tip', pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  output$mobileog_content <- renderUI({
    req(active_content_type() == "mobileog")
    box(
      title = paste("MobileOG Category Distribution for", current_bacterium()$name),
      width = 12, solidHeader = TRUE, status = "warning",
      fluidRow(
        column(12,
               tags$h5("Mobile genetic elements (MGEs) are DNA segments—such as plasmids, phages, integrative, transposable, and conjugative elements—that can move between genomes. We created mobileOG-db, a curated catalog of protein families that define the core life-cycle functions of these MGEs. The results below are derived from MGEs and their corresponding mobileOG classifications.
                       mobileOG-db can be accessed at ", tags$a(href = "https://mobileogdb.flsi.cloud.vt.edu/", "mobileogdb.flsi.cloud.vt.edu", target = "_blank"))
        ),
        column(12,
               downloadButton("download_mobileog_csv", "Download mobile elements", class = "btn-primary", style = "margin-bottom: 10px;")
        ),
        tags$hr(style = "margin-top: 5px; margin-bottom: 15px;")
      ),
      fluidRow(
        column(width = 6,
               plotlyOutput("mobileog_donut_chart", height = "400px")
        ),
        column(width = 6,
               h4("Top MobileOG Categories"),
               DTOutput("mobileog_summary_table")
        )
      ),
      tags$hr(style = "border-top: 1px solid #ccc;"),
      fluidRow(
        column(
          width = 12,
          h4(HTML(paste0("Genes in MobileOG Category: ", tags$span(id="mobileog-selected-type", "(Select a category from the chart)")))),
          DTOutput("mobileog_genes_table")
        )
      )
    )
  })
  
  output$mobileog_donut_chart <- renderPlotly({
    counts <- mobileog_category_counts()
    if (is.null(counts) || nrow(counts) == 0 || sum(counts$Count) == 0) {
      return(plotly_empty() %>%
               layout(
                 title = 'No MobileOG Data Available',
                 annotations = list(x = 0.5, y = 0.5, text = "MobileOG CSV file not loaded or is empty.",
                                    xref = "paper", yref = "paper", showarrow = FALSE, font = list(size = 20, color = "red"))
               ) %>% plotly::event_register("plotly_click")
      )
    }
    
    plot_ly(counts, labels = ~Category, values = ~Count, type = 'pie',
            hole = 0.6,
            marker = list(line = list(color = '#FFFFFF', width = 1)),
            textinfo = 'label+percent',
            insidetextorientation = 'radial',
            hoverinfo = 'text',
            text = ~paste(Category, ':', Count),
            source = "mobileog_donut_clicks") %>%
      layout(title = 'MobileOG Category Distribution',
             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             showlegend = TRUE) %>%
      plotly::event_register("plotly_click")
  })
  
  output$mobileog_summary_table <- renderDT({
    counts <- mobileog_category_counts()
    if (is.null(counts) || nrow(counts) == 0 || sum(counts$Count) == 0) {
      return(datatable(data.frame(Message = "No MobileOG category summary loaded."), rownames = FALSE))
    }
    
    datatable(
      counts,
      options = list(dom = 't', pageLength = 10, autoWidth = TRUE),
      rownames = FALSE,
      colnames = c("Category", "Count")
    )
  })
  
  output$mobileog_genes_table <- renderDT({
    datatable(
      filtered_mobileog_data_reactive(),
      options = list(dom = 'tip', pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  
  output$kegg_sankey <- renderPlotly({
    req(active_content_type() == "kegg")
    df <- current_bacterium_kegg_data()
    req(df)
    
    df <- df[complete.cases(df[, c("Main_Pathway","Sub_pathway","Pathway")]), ]
    df <- unique(df)
    df_s <- df[, c("Main_Pathway","Sub_pathway","Pathway")]
    
    build_links <- function(data, c1, c2) {
      aggregate(list(value = data[[c1]]),
                by = list(source = data[[c1]], target = data[[c2]]), FUN = length)
    }
    links <- rbind(build_links(df_s, "Main_Pathway", "Sub_pathway"),
                   build_links(df_s, "Sub_pathway", "Pathway"))
    
    nodes <- data.frame(name = unique(c(links$source, links$target)),
                        stringsAsFactors = FALSE)
    links$source_id <- match(links$source, nodes$name) - 1
    links$target_id <- match(links$target, nodes$name) - 1
    
    kegg_node_names(nodes$name)          # <-- KEY: expose names to the observer
    
    plot_ly(
      type = "sankey",
      orientation = "h",
      node = list(
        label = nodes$name,
        pad = 15, thickness = 25,
        line = list(color = "white", width = 0.5)
      ),
      link = list(
        source = links$source_id,
        target = links$target_id,
        value  = links$value
      ),
      source = "kegg_sankey_clicks"
    ) %>%
      layout(title = paste("KEGG Pathways for", current_bacterium()$name),
             font = list(size = 12)) %>%
      event_register("plotly_click")
  })
  
  
  
  # 1. Download Prokka TSV
  output$download_prokka_tsv <- downloadHandler(
    filename = function() {
      # e.g., acutalibacter_muris_kb18.tsv
      paste0(get_base_filename(), "_prokka.tsv")
    },
    content = function(file) {
      copy_data_file(get_base_filename(), "_prokka.tsv", file)
    }
  )
  
  output$download_circular_image <- downloadHandler(
    filename = function() {
      # e.g., acutalibacter_muris_kb18.tsv
      paste0(get_base_filename(), "_circular.png")
    },
    content = function(file) {
      copy_data_file(get_base_filename(), "_circular.png", file)
    }
  )
  
  
  # 2. Download Prokka GFF
  output$download_prokka_gff <- downloadHandler(
    filename = function() {
      # e.g., acutalibacter_muris_kb18.gff
      paste0(get_base_filename(), "_prokka.gff")
    },
    content = function(file) {
      copy_data_file(get_base_filename(), "_prokka.gff", file)
    }
  )
  
  # 3. Download Prokka GBK (GenBank)
  output$download_prokka_gbk <- downloadHandler(
    filename = function() {
      # e.g., acutalibacter_muris_kb18.gbk
      paste0(get_base_filename(), "_prokka.gbk")
    },
    content = function(file) {
      copy_data_file(get_base_filename(), "_prokka.gbk", file)
    }
  )
  
  # 4. Download COG table (one row per gene: locus tag, COG category, description)
  # (output id kept as download_cog_gff so the existing button keeps working;
  #  no *_cog.gff files exist -- the COG data lives in *_cog.tsv)
  output$download_cog_gff <- downloadHandler(
    filename = function() paste0(get_base_filename(), "_cog.tsv"),
    content = function(file) copy_data_file(get_base_filename(), "_cog.tsv", file)
  )
  
  
  
  # 1. Download mobile_og TSV
  output$download_mobileog_csv <- downloadHandler(
    filename = function() {
      # e.g., acutalibacter_muris_kb18.tsv
      paste0(get_base_filename(), "_mobileog.csv")
    },
    content = function(file) {
      copy_data_file(get_base_filename(), "_mobileog.csv", file)
    }
  )

  # Genome-scale metabolic model (gapseq SBML)
  output$download_gem_sbml <- downloadHandler(
    filename = function() paste0(get_base_filename(), "_model.xml"),
    content = function(file) copy_data_file(get_base_filename(), "_model.xml", file)
  )

  # 5. Download raw genome FASTA (nucleotide sequence)
  output$download_genome_data <- downloadHandler(
    filename = function() {
      # e.g., acutalibacter_muris_kb18_genome.fasta
      paste0(get_base_filename(), "_genome.fasta")
    },
    content = function(file) {
      copy_data_file(get_base_filename(), "_genome.fasta", file)
    }
  )

}



shinyApp(ui = ui, server = server)
