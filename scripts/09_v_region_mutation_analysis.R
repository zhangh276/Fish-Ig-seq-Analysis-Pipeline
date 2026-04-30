#!/usr/bin/env Rscript

# Step 09: V-region somatic hypermutation analysis
#
# This script summarizes V-region mutation burden based on mutation tables
# generated from VDJRegion clonotypes in Step 05.
#
# Main read-weighted outputs:
#   cloneCount
#   allReads
#   unmutatedReads
#   unmutatedReadsFraction
#   mutatedReads
#   mutatedReadsFraction
#   meanMutationsPerRead
#
# Input:
#   metadata.csv
#   results/05_vdjregion_mutation_tsv/<sample_id>_mutation_IGH.tsv
#
# Output:
#   results/09_v_region_mutation/v_region_mutation_stat.tsv

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

metadata_file <- ifelse(length(args) >= 1, args[1], "metadata.csv")
mutation_dir  <- ifelse(length(args) >= 2, args[2], "results/05_vdjregion_mutation_tsv")
out_dir       <- ifelse(length(args) >= 3, args[3], "results/09_v_region_mutation")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

to_numeric_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", "", x)
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  return(x)
}

# ------------------------------------------------------------
# Read metadata
# ------------------------------------------------------------

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

sample <- read.csv(
  file = metadata_file,
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"sample_id" %in% colnames(sample)) {
  stop("The metadata file must contain a column named 'sample_id'.")
}

sample_ID <- sample$sample_id

# ------------------------------------------------------------
# Initialize output table
# ------------------------------------------------------------

clone_stat <- sample

clone_stat$cloneCount <- 0
clone_stat$allReads <- 0
clone_stat$unmutatedReads <- 0
clone_stat$unmutatedReadsFraction <- 0
clone_stat$mutatedReads <- 0
clone_stat$mutatedReadsFraction <- 0
clone_stat$meanMutationsPerRead <- 0

# ------------------------------------------------------------
# Required columns
# ------------------------------------------------------------

length_cols <- c(
  "nLengthFR1",
  "nLengthCDR1",
  "nLengthFR2",
  "nLengthCDR2",
  "nLengthFR3"
)

mutation_cols <- c(
  "nMutationsCountFR1Substitutions",
  "nMutationsCountCDR1Substitutions",
  "nMutationsCountFR2Substitutions",
  "nMutationsCountCDR2Substitutions",
  "nMutationsCountFR3Substitutions"
)

required_cols <- c(
  "cloneId",
  "readCount",
  "readFraction",
  length_cols,
  mutation_cols
)

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------

for (i in seq_along(sample_ID)) {
  sid <- sample_ID[i]

  message("Processing sample: ", sid)

  mutation_file <- file.path(mutation_dir, paste0(sid, "_mutation_IGH.tsv"))

  if (!file.exists(mutation_file)) {
    mutation_file_alt <- file.path(mutation_dir, paste0(sid, "_mutation.tsv"))

    if (file.exists(mutation_file_alt)) {
      mutation_file <- mutation_file_alt
    } else {
      stop(
        "Mutation file not found for sample ",
        sid,
        ". Tried:\n",
        file.path(mutation_dir, paste0(sid, "_mutation_IGH.tsv")),
        "\n",
        mutation_file_alt
      )
    }
  }

  full_length_mutation <- read.table(
    file = mutation_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  missing_cols <- required_cols[!required_cols %in% colnames(full_length_mutation)]

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ",
      mutation_file,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  full_length_mutation_sub <- full_length_mutation[, required_cols]

  # Remove rows containing region_not_covered, consistent with the reference code.
  full_length_mutation_sub <- full_length_mutation_sub %>%
    filter(!if_any(everything(), ~ str_detect(as.character(.), "region_not_covered")))

  clone_stat$cloneCount[i] <- nrow(full_length_mutation_sub)

  if (nrow(full_length_mutation_sub) == 0) {
    message("  No covered V-region clonotypes retained for sample: ", sid)
    next
  }

  # Convert readCount and mutation-related columns to numeric.
  full_length_mutation_sub$readCount <- to_numeric_safe(full_length_mutation_sub$readCount)

  full_length_mutation_sub[, length_cols] <- lapply(
    full_length_mutation_sub[, length_cols],
    to_numeric_safe
  )

  full_length_mutation_sub[, mutation_cols] <- lapply(
    full_length_mutation_sub[, mutation_cols],
    to_numeric_safe
  )

  # Total analysed V-region length and total substitution mutations.
  full_length_mutation_sub$nLength <- rowSums(
    full_length_mutation_sub[, length_cols],
    na.rm = TRUE
  )

  full_length_mutation_sub$nMutations <- rowSums(
    full_length_mutation_sub[, mutation_cols],
    na.rm = TRUE
  )

  # Read-weighted SHM burden.
  full_length_mutation_sub$SHM_reads <- (
    full_length_mutation_sub$nMutations *
      full_length_mutation_sub$readCount
  )

  all_reads <- sum(full_length_mutation_sub$readCount, na.rm = TRUE)

  unmutated_reads <- sum(
    full_length_mutation_sub$readCount[full_length_mutation_sub$nMutations == 0],
    na.rm = TRUE
  )

  mutated_reads <- sum(
    full_length_mutation_sub$readCount[full_length_mutation_sub$nMutations != 0],
    na.rm = TRUE
  )

  clone_stat$allReads[i] <- all_reads
  clone_stat$unmutatedReads[i] <- unmutated_reads
  clone_stat$mutatedReads[i] <- mutated_reads

  if (all_reads > 0) {
    clone_stat$unmutatedReadsFraction[i] <- unmutated_reads / all_reads
    clone_stat$mutatedReadsFraction[i] <- mutated_reads / all_reads
    clone_stat$meanMutationsPerRead[i] <- sum(
      full_length_mutation_sub$SHM_reads,
      na.rm = TRUE
    ) / all_reads
  } else {
    clone_stat$unmutatedReadsFraction[i] <- 0
    clone_stat$mutatedReadsFraction[i] <- 0
    clone_stat$meanMutationsPerRead[i] <- 0
  }

  # Optional: write the filtered mutation table for each sample.
  write.table(
    full_length_mutation_sub,
    file = file.path(out_dir, paste0(sid, "_V_region_mutation_filtered.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message("  cloneCount: ", clone_stat$cloneCount[i])
  message("  allReads: ", clone_stat$allReads[i])
  message("  unmutatedReadsFraction: ", round(clone_stat$unmutatedReadsFraction[i], 4))
  message("  meanMutationsPerRead: ", round(clone_stat$meanMutationsPerRead[i], 4))
}

# ------------------------------------------------------------
# Write summary output
# ------------------------------------------------------------

write.table(
  clone_stat,
  file = file.path(out_dir, "v_region_mutation_stat.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Step 09 finished: V-region mutation analysis completed.")
message("Output directory: ", out_dir)
