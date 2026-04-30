#!/usr/bin/env Rscript

# Install R packages required for downstream analysis steps.
#
# Required by:
#   Step 06: IGHV and VJ usage analysis
#   Step 07: pairwise repertoire correlation analysis
#   Step 08: CDR3 clonotype network analysis
#   Step 09: V-region mutation analysis

cran_mirror <- "https://cloud.r-project.org"

required_packages <- c(
  "stringdist",
  "igraph",
  "dplyr",
  "tidyr",
  "data.table",
  "stringr"
)

installed_packages <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!pkg %in% installed_packages) {
    message("Installing R package: ", pkg)
    install.packages(pkg, repos = cran_mirror)
  } else {
    message("R package already installed: ", pkg)
  }
}

message("All required R packages are installed.")
