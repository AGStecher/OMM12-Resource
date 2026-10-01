library(readr)
library(dplyr)

get_genome_info_html <- function(bacterium_name_clean) {
  
  # Path to the data file. It is assumed to be in the 'data' folder.
  data_file_path <- "data/genome_info.tsv"
  
  # Check if the data file exists
  if (!file.exists(data_file_path)) {
    return(
      paste0(
        "<h2>Genome Information Not Found</h2>",
        "<p>The required data file (", data_file_path, ") could not be found.</p>"
      )
    )
  }
  
  # Read the data and filter for the specific bacterium
  genome_info_df <- read_tsv(data_file_path, show_col_types = FALSE) %>%
    filter(bacterium_name_clean == bacterium_name_clean) # The second bacterium_name_clean is the one from the function argument
  
  # Check if data exists for the current bacterium
  if (nrow(genome_info_df) == 0) {
    return(
      paste0(
        "<h2>Genome Information Not Found</h2>",
        "<p>No specific genome data is available for ", gsub("_", " ", bacterium_name_clean), " in the data file.</p>"
      )
    )
  }
  
  data <- as.list(genome_info_df[1, ])
  
  # --- Build the HTML content using tags and variables from the data frame ---
  html_content <- paste0(
    "<h2>Genome Information</h2>",
    "<div class='genome-info-container'>",
    "<div class='genome-figure'>",
    "<img src='", data$figure_path, "' alt='Genome Map' style='width: 100%; border: 1px solid #ddd;'>",
    "<p class='figure-caption'>A schematic map of the genome for ", data$display_name, ".</p>",
    "</div>",
    "<div class='genome-details'>",
    "<h3>Key Statistics:</h3>",
    "<ul>",
    "<li><strong>Genome Size:</strong> ", data$size_mbp, " Mbp</li>",
    "<li><strong>GC Content:</strong> ", data$gc_content, " %</li>",
    "<li><strong>Number of Genes:</strong> ", data$num_genes, "</li>",
    "</ul>",
    "<h3>Description:</h3>",
    "<p>", data$description, "</p>",
    "<h3>External Resources:</h3>",
    "<ul>",
    "<li><a href='", data$ncbi_link, "' target='_blank'>View on NCBI Genome</a></li>",
    "<li><a href='", data$biocyc_link, "' target='_blank'>View on BioCyc</a></li>",
    "</ul>",
    "</div>",
    "</div>"
  )
  
  return(html_content)
}