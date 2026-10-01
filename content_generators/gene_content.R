# This function will generate the specific HTML content for Gene Info
get_gene_info_html <- function(bacterium_clean_name) {
  # You can make this more dynamic, e.g., load gene lists from a file
  
  html_content <- paste0(
    "<h2>Gene Information</h2>",
    "<p>This page lists key genes and their functions for <strong>", gsub("_", " ", bacterium_clean_name), "</strong>.</p>",
    "<ul>",
    "<li>Gene A: Function related to metabolism.</li>",
    "<li>Gene B: Involved in stress response.</li>",
    "<li>Gene C: Virulence factor.</li>",
    "</ul>",
    "<p>For a complete list, please refer to the gene database link.</p>"
  )
  
  return(html_content)
}