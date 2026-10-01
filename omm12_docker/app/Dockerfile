# OMM12 Genomic Resource -- Shiny app image
# Modeled on daisybio/namco's app/Dockerfile, trimmed down: no conda/picrust2/
# usearch/FastTree needed since this app runs no bioinformatics tools live --
# it only reads precomputed annotation files and renders them.
#
# Build context is the OMM12_website folder (see docker-compose.yml);
# .dockerignore in that folder decides what gets copied in.
#
# R version: set this to your local R.version.string (run it in RStudio).
# The versioned rocker/r-ver tags pin CRAN to a dated Posit Package Manager
# snapshot, so package versions stay fixed and install as fast binaries.
ARG R_VERSION=4.4.2
FROM rocker/r-ver:${R_VERSION}

# System libraries needed at build/run time by the R packages below
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libcairo2-dev \
    libxt-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libglpk-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# CRAN packages. No repos= argument on purpose: this uses the dated snapshot
# rocker configured for this R version (see note above).
RUN R -e "install.packages(c( \
    'shiny', 'shinydashboard', 'shinydashboardPlus', 'shinyjs', \
    'ggplot2', 'plotly', 'dplyr', 'readr', 'ape', 'DT', 'stringr', \
    'tidyr', 'wordcloud2', 'ggiraph', 'shinySearchbar', 'shinyWidgets', \
    'networkD3', 'xml2', 'shinycssloaders', 'gggenes', 'BiocManager' \
  ))"

# Bioconductor packages (not on CRAN)
RUN R -e "BiocManager::install(c('ggtree', 'ggtreeExtra'), update = FALSE, ask = FALSE)"

# Fail the build right here if any package the app loads is missing
RUN R -e "pkgs <- c('shiny','shinydashboard','shinydashboardPlus','shinyjs','ggplot2','plotly','dplyr','readr','ape','ggtree','DT','stringr','tidyr','wordcloud2','ggiraph','ggtreeExtra','shinySearchbar','shinyWidgets','networkD3','xml2','shinycssloaders','gggenes'); miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; if (length(miss)) stop('Missing R packages: ', paste(miss, collapse = ', '))"

COPY app.R /srv/shiny-server/app.R
COPY data /srv/shiny-server/data
COPY www /srv/shiny-server/www

WORKDIR /srv/shiny-server

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', host = '0.0.0.0', port = 3838)"]
