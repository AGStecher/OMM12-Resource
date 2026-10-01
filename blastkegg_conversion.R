# -----------------------------
# Convert KEGG-module-like text
# -----------------------------

input_file  <- "data/acutalibacter_muris_kb18_16_kegg.txt"     # your raw text file
output_file <- "data/acutalibacter_muris_kb18_16_kegg_converted_output.tsv"

# Read file
lines <- readLines(input_file)

# Storage vectors
heading <- ""
ko <- ""
results <- data.frame(GeneID = character(),
                      KO = character(),
                      Heading = character(),
                      stringsAsFactors = FALSE)

for (ln in lines) {
  
  # Detect section heading (starts with Mxxxxx)
  if (grepl("^M\\d{5}", ln)) {
    heading <- trimws(ln)
  }
  
  # Detect KO term (starts with Kxxxxx)
  else if (grepl("^K\\d{5}", ln)) {
    ko <- trimws(ln)
  }
  
  # Detect gene IDs (start with G…)
  else if (grepl("^G[A-Za-z0-9_\\.]+", ln)) {
    gene_ids <- unlist(strsplit(ln, ","))
    gene_ids <- trimws(gene_ids)
    
    for (g in gene_ids) {
      results <- rbind(
        results,
        data.frame(GeneID = g, KO = ko, Heading = heading, stringsAsFactors = FALSE)
      )
    }
  }
}

# Write output
write.table(results, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Conversion complete! File saved as:", output_file, "\n")
