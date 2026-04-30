#!/usr/bin/env Rscript

# Step 06: IGHV and VJ usage analysis
#
# This script calculates IGHV usage and VJ combination usage
# under two weighting schemes:
#   1. V3J_weight: each CDR3 clonotype is counted once
#   2. reads_weight: each CDR3 clonotype is weighted by readCount
#
# Required input:
#   metadata.csv
#   results/03_cdr3_tsv/<sample_id>_CDR3_clones_IGH.tsv
#   example_data/reference/v-genes.IGH.fasta
#   example_data/reference/j-genes.IGH.fasta

args <- commandArgs(trailingOnly = TRUE)

metadata_file <- ifelse(length(args) >= 1, args[1], "metadata.csv")
cdr3_dir      <- ifelse(length(args) >= 2, args[2], "results/03_cdr3_tsv")
ref_dir       <- ifelse(length(args) >= 3, args[3], "example_data/reference")
out_dir       <- ifelse(length(args) >= 4, args[4], "results/06_ighv_vj_usage")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

v_fasta <- file.path(ref_dir, "v-genes.IGH.fasta")
j_fasta <- file.path(ref_dir, "j-genes.IGH.fasta")

# ------------------------------------------------------------
# Functions
# ------------------------------------------------------------

remove_allele <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^(.*?)\\*.*$", "\\1", x)
  x <- trimws(x)
  x[x == ""] <- NA
  return(x)
}

normalize_mixcr_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)

  # Keep the first hit if multiple hits exist.
  x <- sapply(strsplit(x, ",", fixed = TRUE), `[`, 1)
  x <- trimws(x)

  # Keep only the part before the first "*".
  x <- remove_allele(x)

  return(x)
}

extract_gene_from_fasta <- function(fasta_file, gene_prefix) {
  if (!file.exists(fasta_file)) {
    stop("Reference FASTA file not found: ", fasta_file)
  }

  lines <- readLines(fasta_file, warn = FALSE)
  headers <- lines[grepl("^>", lines)]

  if (length(headers) == 0) {
    stop("No FASTA headers found in: ", fasta_file)
  }

  headers <- sub("^>", "", headers)
  gene_names <- c()

  for (h in headers) {
    h <- trimws(h)

    fields <- unlist(strsplit(h, "[|[:space:]]+"))
    fields <- trimws(fields)
    fields <- fields[fields != ""]

    gene_field <- fields[grepl(gene_prefix, fields)]

    if (length(gene_field) > 0) {
      gene_name <- gene_field[1]
    } else {
      gene_name <- fields[1]
    }

    gene_name <- remove_allele(gene_name)

    if (!is.na(gene_name) && gene_name != "") {
      gene_names <- c(gene_names, gene_name)
    }
  }

  gene_names <- unique(gene_names)
  gene_names <- gene_names[!is.na(gene_names) & gene_names != ""]
  gene_names <- sort(gene_names)

  return(gene_names)
}

normalize_freq <- function(x) {
  x[is.na(x)] <- 0
  total <- sum(x)

  if (total == 0) {
    return(x)
  } else {
    return(x / total)
  }
}

# ------------------------------------------------------------
# Read metadata and reference genes
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

all_V_segment <- data.frame(
  V_segment = extract_gene_from_fasta(v_fasta, "IGHV"),
  stringsAsFactors = FALSE
)

all_J_segment <- data.frame(
  J_segment = extract_gene_from_fasta(j_fasta, "IGHJ"),
  stringsAsFactors = FALSE
)

all_VJ_combination <- expand.grid(
  V_segment = all_V_segment$V_segment,
  J_segment = all_J_segment$J_segment,
  stringsAsFactors = FALSE
)

all_VJ_combination$VJ_combination <- paste(
  all_VJ_combination$V_segment,
  all_VJ_combination$J_segment,
  sep = "_"
)

all_VJ_combination <- data.frame(
  VJ_combination = sort(all_VJ_combination$VJ_combination),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Initialize output matrices
# ------------------------------------------------------------

IGHV_V3J_weight_matrix <- all_V_segment
IGHV_reads_weight_matrix <- all_V_segment

VJ_V3J_weight_matrix <- all_VJ_combination
VJ_reads_weight_matrix <- all_VJ_combination

clone_stat <- sample
clone_stat$V3J <- 0

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------

for (i in seq_along(sample_ID)) {
  sid <- sample_ID[i]
  message("Processing sample: ", sid)

  cdr3_file <- file.path(cdr3_dir, paste0(sid, "_CDR3_clones_IGH.tsv"))

  if (!file.exists(cdr3_file)) {
    stop("CDR3 clonotype file not found: ", cdr3_file)
  }

  CDR3_clone <- read.table(
    file = cdr3_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  required_cols <- c("readCount", "allVHitsWithScore", "allJHitsWithScore")
  missing_cols <- required_cols[!required_cols %in% colnames(CDR3_clone)]

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ",
      cdr3_file,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  CDR3_clone_sub <- CDR3_clone[, required_cols]

  clone_stat$V3J[i] <- nrow(CDR3_clone_sub)

  CDR3_clone_sub$V_segment <- normalize_mixcr_gene(
    CDR3_clone_sub$allVHitsWithScore
  )

  CDR3_clone_sub$J_segment <- normalize_mixcr_gene(
    CDR3_clone_sub$allJHitsWithScore
  )

  CDR3_clone_sub$VJ_combination <- paste(
    CDR3_clone_sub$V_segment,
    CDR3_clone_sub$J_segment,
    sep = "_"
  )

  CDR3_clone_sub$readCount <- as.numeric(
    gsub(",", "", CDR3_clone_sub$readCount)
  )

  CDR3_clone_sub$readCount[is.na(CDR3_clone_sub$readCount)] <- 0

  # ==========================================================
  # 1. IGHV usage under V3J_weight
  # ==========================================================

  V3J_V_freq <- aggregate(
    rep(1, nrow(CDR3_clone_sub)),
    by = list(V_segment = CDR3_clone_sub$V_segment),
    FUN = sum
  )

  colnames(V3J_V_freq) <- c("V_segment", "V_segment_freq")

  tmp <- merge(
    all_V_segment,
    V3J_V_freq,
    by = "V_segment",
    all.x = TRUE,
    sort = FALSE
  )

  tmp <- tmp[match(all_V_segment$V_segment, tmp$V_segment), ]
  tmp$V_segment_freq <- normalize_freq(tmp$V_segment_freq)

  IGHV_V3J_weight_matrix[[paste0(sid, "_IG")]] <- tmp$V_segment_freq

  # ==========================================================
  # 2. IGHV usage under reads_weight
  # ==========================================================

  reads_V_freq <- aggregate(
    CDR3_clone_sub$readCount,
    by = list(V_segment = CDR3_clone_sub$V_segment),
    FUN = sum
  )

  colnames(reads_V_freq) <- c("V_segment", "V_segment_freq")

  tmp <- merge(
    all_V_segment,
    reads_V_freq,
    by = "V_segment",
    all.x = TRUE,
    sort = FALSE
  )

  tmp <- tmp[match(all_V_segment$V_segment, tmp$V_segment), ]
  tmp$V_segment_freq <- normalize_freq(tmp$V_segment_freq)

  IGHV_reads_weight_matrix[[paste0(sid, "_IG")]] <- tmp$V_segment_freq

  # ==========================================================
  # 3. VJ usage under V3J_weight
  # ==========================================================

  V3J_VJ_freq <- aggregate(
    rep(1, nrow(CDR3_clone_sub)),
    by = list(VJ_combination = CDR3_clone_sub$VJ_combination),
    FUN = sum
  )

  colnames(V3J_VJ_freq) <- c("VJ_combination", "VJ_combination_freq")

  tmp <- merge(
    all_VJ_combination,
    V3J_VJ_freq,
    by = "VJ_combination",
    all.x = TRUE,
    sort = FALSE
  )

  tmp <- tmp[match(all_VJ_combination$VJ_combination, tmp$VJ_combination), ]
  tmp$VJ_combination_freq <- normalize_freq(tmp$VJ_combination_freq)

  VJ_V3J_weight_matrix[[paste0(sid, "_IG")]] <- tmp$VJ_combination_freq

  # ==========================================================
  # 4. VJ usage under reads_weight
  # ==========================================================

  reads_VJ_freq <- aggregate(
    CDR3_clone_sub$readCount,
    by = list(VJ_combination = CDR3_clone_sub$VJ_combination),
    FUN = sum
  )

  colnames(reads_VJ_freq) <- c("VJ_combination", "VJ_combination_freq")

  tmp <- merge(
    all_VJ_combination,
    reads_VJ_freq,
    by = "VJ_combination",
    all.x = TRUE,
    sort = FALSE
  )

  tmp <- tmp[match(all_VJ_combination$VJ_combination, tmp$VJ_combination), ]
  tmp$VJ_combination_freq <- normalize_freq(tmp$VJ_combination_freq)

  VJ_reads_weight_matrix[[paste0(sid, "_IG")]] <- tmp$VJ_combination_freq
}

# ------------------------------------------------------------
# Write output files
# ------------------------------------------------------------

write.table(
  clone_stat,
  file = file.path(out_dir, "clone_stat.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  IGHV_V3J_weight_matrix,
  file = file.path(out_dir, "IGHV_V3J_weight_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  IGHV_reads_weight_matrix,
  file = file.path(out_dir, "IGHV_reads_weight_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  VJ_V3J_weight_matrix,
  file = file.path(out_dir, "VJ_V3J_weight_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  VJ_reads_weight_matrix,
  file = file.path(out_dir, "VJ_reads_weight_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Step 06 finished: IGHV and VJ usage analysis completed.")
message("Output directory: ", out_dir)
