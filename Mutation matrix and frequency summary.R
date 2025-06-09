###################################################
#        Mutation matrix and frequency summary    #
#             Version 1.0 (June 2025)             #
#        Written by Adam A. Capoferri, PhD        #               
###################################################

# This script can be used two fold. For mutation frequency you will want to run this without a collapsed.fasta file. The output is number/frequency for A>G, G>A, A<>G, etc mutations, total number of mutations, total A, C, T, G changes, and the Transition/Transversion ratio. The other file that is generated is a data matrix for each nucleotide frequency for all sequences relative to the reference sequence (ensure it is the first sequence in the FASTA file). If there is no change "." other the nucleotide change is shown as well as any gaps. For uncollapsed sequences this may not be a completely necessary table. However, it is for collapsed sequences that will then be used to generate highlighter plots. This is because the associated R script uses this data matrix based csv file as the input. 

# As always, ensure all packages are installed first.

library(Biostrings)
library(tidyverse)

# ---- Load aligned sequences ----
fasta_file <- "201403_alignment with master_collapsed.fasta"  # <-- Change to your file, needs to be non-collapsed
seqs <- readDNAStringSet(fasta_file)
seq_mat <- as.matrix(seqs)

ref <- seq_mat[1, ]
query_mat <- seq_mat[-1, , drop = FALSE]
mutation_types <- c("A", "T", "C", "G")

# ---- Build mutation matrix (match = ".", mismatch = base) ----
comparison <- apply(query_mat, 1, function(row) {
  ifelse(row == ref, ".", as.character(row))
})
comparison <- t(comparison)
rownames(comparison) <- names(seqs)[-1]
colnames(comparison) <- seq_len(ncol(query_mat))

# ---- Mutation counts per position ----
mut_counts <- lapply(seq_len(ncol(query_mat)), function(j) {
  ref_base <- as.character(ref[j])
  col <- as.character(query_mat[, j])
  muts <- col[col != ref_base]
  table(factor(muts[muts %in% mutation_types], levels = mutation_types))
})
mut_counts_matrix <- do.call(rbind, mut_counts)
mut_counts_matrix <- t(mut_counts_matrix)  # Rows = A/T/C/G, Cols = positions
rownames(mut_counts_matrix) <- paste0("count_", mutation_types)
colnames(mut_counts_matrix) <- seq_len(ncol(query_mat))

# ---- Total mutation count per position ----
total_mut <- colSums(mut_counts_matrix)

# ---- Substitution summary ----
subs_to_count <- c(
  "A>G", "G>A", "C>T", "T>C",  # Transitions
  "A>C", "C>A", "A>T", "T>A",  # Transversions
  "C>G", "G>C", "G>T", "T>G",  # Transversions
  "A<>G", "C<>T", "A<>C", "A<>T", "C<>G", "G<>T"  # Symmetric
)
sub_counts <- setNames(rep(0, length(subs_to_count)), subs_to_count)
base_freqs <- setNames(rep(0, 4), paste0("freq_", mutation_types))
base_totals <- setNames(rep(0, 4), paste0("total_", mutation_types))

# Ti/Tv groups
transitions <- c("A>G", "G>A", "C>T", "T>C")
transversions <- c("A>C", "C>A", "A>T", "T>A", "C>G", "G>C", "G>T", "T>G")

ti_count <- 0
tv_count <- 0

# Tally substitutions
for (i in 2:nrow(seq_mat)) {
  for (j in seq_along(ref)) {
    ref_base <- as.character(ref[j])
    query_base <- as.character(seq_mat[i, j])
    
    if (ref_base != query_base &&
        ref_base %in% mutation_types &&
        query_base %in% mutation_types) {
      
      pair <- paste0(ref_base, ">", query_base)
      reverse <- paste0(query_base, ">", ref_base)
      symmetric <- paste0(sort(c(ref_base, query_base)), collapse = "<>")
      
      if (pair %in% names(sub_counts)) sub_counts[pair] <- sub_counts[pair] + 1
      if (reverse %in% names(sub_counts)) sub_counts[reverse] <- sub_counts[reverse] + 1
      if (symmetric %in% names(sub_counts)) sub_counts[symmetric] <- sub_counts[symmetric] + 1
      
      if (pair %in% transitions) ti_count <- ti_count + 1
      if (pair %in% transversions) tv_count <- tv_count + 1
    }
    
    if (query_base %in% mutation_types) {
      base_freqs[paste0("freq_", query_base)] <- base_freqs[paste0("freq_", query_base)] + 1
      base_totals[paste0("total_", query_base)] <- base_totals[paste0("total_", query_base)] + 1
    }
  }
}

# ---- Substitution frequencies (as percentages) ----
# Improved substitution frequency calculation
sub_freqs <- setNames(rep(NA, length(subs_to_count)), paste0("freq_", subs_to_count))

for (sub in names(sub_counts)) {
  from_base <- strsplit(sub, "[><]")[[1]][1]
  total_from <- base_freqs[paste0("freq_", from_base)]
  
  # Only compute if total_from > 0, and use 4 decimal places (0.0001%)
  sub_freqs[paste0("freq_", sub)] <- if (total_from > 0) {
    signif(100 * sub_counts[sub] / total_from, 4)
  } else {
    NA
  }
}

# ---- Build mutation matrix ----
mutation_matrix <- as.data.frame(comparison)
mutation_matrix <- rbind(
  reference = as.list(as.character(ref)),
  mutation_matrix,
  count_A = as.list(mut_counts_matrix["count_A", ]),
  count_T = as.list(mut_counts_matrix["count_T", ]),
  count_C = as.list(mut_counts_matrix["count_C", ]),
  count_G = as.list(mut_counts_matrix["count_G", ]),
  total_mut = as.list(total_mut)
)

# ---- Substitution summary table ----
sub_summary <- rbind(counts = sub_counts, frequencies = sub_freqs)
sub_summary <- as.data.frame(t(sub_summary))

# ---- Add totals and Ti/Tv info to summary ----
summary_row <- data.frame(
  counts = c(
    total_mutations = sum(sub_counts),
    total_A = base_totals["total_A"],
    total_C = base_totals["total_C"],
    total_G = base_totals["total_G"],
    total_T = base_totals["total_T"],
    transitions = ti_count,
    transversions = tv_count,
    Ti_Tv_ratio = round(if (tv_count > 0) ti_count / tv_count else NA, 3)
  ),
  frequencies = c(rep(NA, 8))
)

sub_summary_extended <- rbind(sub_summary, summary_row)

# ---- Write to disk ----
#Change the names of the files as needed.

write.csv(mutation_matrix, "201403_alignment with master_collapsed_mutation_matrix.csv", row.names = TRUE)
write.csv(sub_summary_extended, "201403_alignment with master_collapsed_summary_extended.csv", row.names = TRUE)

cat("✔ Written: mutation_collapsed_matrix.csv and substitution_summary.csv\n")
### End ###
