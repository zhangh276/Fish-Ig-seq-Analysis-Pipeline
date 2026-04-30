#!/usr/bin/env Rscript

# Step 08: CDR3 clonotype network analysis
#
# This script constructs CDR3 clonotype networks using Levenshtein distance.
#
# Input:
#   metadata.csv
#   results/03_cdr3_tsv/<sample_id>_CDR3_clones_IGH.tsv
#
# Output:
#   results/08_cdr3_network/<sample_id>_CDR3_network_nodes.tsv
#   results/08_cdr3_network/<sample_id>_CDR3_network_edges.tsv
#   results/08_cdr3_network/<sample_id>_network_LD1_CDR3_nSeq.pdf
#   results/08_cdr3_network/clone_network_stat.tsv
#
# Network definition:
#   Nodes: unique CDR3 clonotypes, defined by V gene + nSeqCDR3 + J gene
#   Edges: pairs of clonotypes with Levenshtein distance = 1 between nSeqCDR3 sequences

suppressPackageStartupMessages({
  library(stringdist)
  library(igraph)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

metadata_file <- ifelse(length(args) >= 1, args[1], "metadata.csv")
cdr3_dir      <- ifelse(length(args) >= 2, args[2], "results/03_cdr3_tsv")
out_dir       <- ifelse(length(args) >= 3, args[3], "results/08_cdr3_network")
max_nodes     <- ifelse(length(args) >= 4, as.integer(args[4]), 3000)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

normalize_mixcr_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)

  # Keep the first hit if multiple hits exist.
  x <- sapply(strsplit(x, ",", fixed = TRUE), `[`, 1)
  x <- trimws(x)

  # Remove allele and score information after the first "*".
  # Example: IGHV1-1*01(123.4) -> IGHV1-1
  x <- sub("^(.*?)\\*.*$", "\\1", x)

  x <- trimws(x)
  x[x == ""] <- NA
  return(x)
}

clean_read_count <- function(x) {
  x <- as.numeric(gsub(",", "", x))
  x[is.na(x)] <- 0
  return(x)
}

build_ld1_edges <- function(seq_ids, cdr3_seqs) {
  n <- length(cdr3_seqs)

  if (n < 2) {
    return(data.frame(query = character(), target = character(), LD = numeric()))
  }

  edges <- list()
  edge_index <- 1

  for (i in 1:(n - 1)) {
    s1 <- cdr3_seqs[i]
    s2_vec <- cdr3_seqs[(i + 1):n]

    # Levenshtein distance = 1 is only possible when sequence length differs by <= 1.
    len_diff <- abs(nchar(s1) - nchar(s2_vec))
    candidate_idx <- which(len_diff <= 1)

    if (length(candidate_idx) == 0) {
      next
    }

    j_indices <- (i + 1):n
    j_candidates <- j_indices[candidate_idx]
    s2_candidates <- cdr3_seqs[j_candidates]

    d <- stringdist(s1, s2_candidates, method = "lv")

    hit <- which(d == 1)

    if (length(hit) > 0) {
      edges[[edge_index]] <- data.frame(
        query = seq_ids[i],
        target = seq_ids[j_candidates[hit]],
        LD = d[hit],
        stringsAsFactors = FALSE
      )
      edge_index <- edge_index + 1
    }
  }

  if (length(edges) == 0) {
    return(data.frame(query = character(), target = character(), LD = numeric()))
  }

  edge_df <- do.call(rbind, edges)
  return(edge_df)
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

sample$V3J <- 0
sample$lcc <- 0
sample$edges <- 0
sample$network_nodes <- 0

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

  CDR3_clone_IGM <- read.table(
    file = cdr3_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  required_cols <- c(
    "nSeqCDR3",
    "readCount",
    "allVHitsWithScore",
    "allJHitsWithScore"
  )

  missing_cols <- required_cols[!required_cols %in% colnames(CDR3_clone_IGM)]

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ",
      cdr3_file,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  CDR3_clone_sub_IGM <- CDR3_clone_IGM[, required_cols]

  CDR3_clone_sub_IGM$allVHitsWithScore <- normalize_mixcr_gene(
    CDR3_clone_sub_IGM$allVHitsWithScore
  )

  CDR3_clone_sub_IGM$allJHitsWithScore <- normalize_mixcr_gene(
    CDR3_clone_sub_IGM$allJHitsWithScore
  )

  CDR3_clone_sub_IGM$readCount <- clean_read_count(
    CDR3_clone_sub_IGM$readCount
  )

  CDR3_clone_sub_IGM <- CDR3_clone_sub_IGM %>%
    filter(
      !is.na(nSeqCDR3),
      nSeqCDR3 != "",
      !is.na(allVHitsWithScore),
      !is.na(allJHitsWithScore)
    )

  CDR3_clone_sub_IGM$V3J <- paste(
    CDR3_clone_sub_IGM$allVHitsWithScore,
    CDR3_clone_sub_IGM$nSeqCDR3,
    CDR3_clone_sub_IGM$allJHitsWithScore,
    sep = "_"
  )

  # Merge duplicated V3J clonotypes if present.
  CDR3_clone_sub_IGM <- CDR3_clone_sub_IGM %>%
    group_by(V3J, nSeqCDR3, allVHitsWithScore, allJHitsWithScore) %>%
    summarise(readCount = sum(readCount), .groups = "drop")

  sample$V3J[i] <- nrow(CDR3_clone_sub_IGM)

  message("  Number of V3J clonotypes: ", nrow(CDR3_clone_sub_IGM))

  if (nrow(CDR3_clone_sub_IGM) > max_nodes) {
    message(
      "  Number of clonotypes exceeds max_nodes = ",
      max_nodes,
      ". Keeping the top ",
      max_nodes,
      " clonotypes by readCount for network analysis."
    )

    CDR3_clone_sub_IGM <- CDR3_clone_sub_IGM %>%
      arrange(desc(readCount)) %>%
      slice_head(n = max_nodes)
  }

  CDR3_clone_sub_IGM$isotype <- "IGM"

  merged_df <- CDR3_clone_sub_IGM[, c("V3J", "nSeqCDR3", "readCount", "isotype")]
  merged_df$seq_ID <- paste0("seq_", seq_len(nrow(merged_df)))

  merged_df <- merged_df[, c("seq_ID", "nSeqCDR3", "V3J", "readCount", "isotype")]

  node_file <- file.path(out_dir, paste0(sid, "_CDR3_network_nodes.tsv"))

  write.table(
    merged_df,
    file = node_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  message("  Node table written to: ", node_file)

  # ----------------------------------------------------------
  # Compute Levenshtein-distance-1 edges
  # ----------------------------------------------------------

  dist_LD_unique <- build_ld1_edges(
    seq_ids = merged_df$seq_ID,
    cdr3_seqs = merged_df$nSeqCDR3
  )

  edge_file <- file.path(out_dir, paste0(sid, "_CDR3_network_edges.tsv"))

  write.table(
    dist_LD_unique,
    file = edge_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  message("  Edge table written to: ", edge_file)

  sample$edges[i] <- nrow(dist_LD_unique)

  # ----------------------------------------------------------
  # Build graph and calculate LCC
  # ----------------------------------------------------------

  if (nrow(dist_LD_unique) > 0) {
    g <- graph_from_data_frame(
      dist_LD_unique[, c("query", "target")],
      vertices = merged_df,
      directed = FALSE
    )
  } else {
    g <- make_empty_graph(n = nrow(merged_df), directed = FALSE)
    V(g)$name <- merged_df$seq_ID
    V(g)$nSeqCDR3 <- merged_df$nSeqCDR3
    V(g)$V3J <- merged_df$V3J
    V(g)$readCount <- merged_df$readCount
    V(g)$isotype <- merged_df$isotype
  }

  sample$network_nodes[i] <- vcount(g)

  comp <- components(g)

  if (length(comp$csize) > 0) {
    lcc_id <- which.max(comp$csize)
    g_lcc <- induced_subgraph(g, vids = which(comp$membership == lcc_id))
    sample$lcc[i] <- vcount(g_lcc)
  } else {
    sample$lcc[i] <- 0
  }

  message("  LCC size: ", sample$lcc[i])

  # ----------------------------------------------------------
  # Plot network
  # ----------------------------------------------------------

  V(g)$label <- NA
  V(g)$size <- log2(as.numeric(V(g)$readCount) + 1) / 1.5

  V(g)$color <- ifelse(
    V(g)$isotype == "IGM",
    "#F77F7F",
    "black"
  )

  pdf_file <- file.path(out_dir, paste0(sid, "_network_LD1_CDR3_nSeq.pdf"))

  pdf(file = pdf_file, width = 3, height = 3)

  if (vcount(g) > 0) {
    layout <- layout_with_fr(g, niter = 2000, grid = "nogrid")

    plot(
      g,
      layout = layout,
      vertex.frame.color = "transparent",
      edge.color = "black",
      vertex.label = NA
    )
  }

  dev.off()

  message("  Network plot written to: ", pdf_file)
}

# ------------------------------------------------------------
# Write summary statistics
# ------------------------------------------------------------

stat_file <- file.path(out_dir, "clone_network_stat.tsv")

write.table(
  sample,
  file = stat_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Step 08 finished: CDR3 clonotype network analysis completed.")
message("Output directory: ", out_dir)
