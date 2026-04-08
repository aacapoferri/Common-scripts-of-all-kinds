###################################################
#        Mutation matrix and frequency summary    #
#             Version 2.4 (April 2026)            #
#        Correct directional mutation counting    #
#         Written by Adam A. Capoferri, PhD       #
###################################################

# load libraries
library(Biostrings)
library(tidyverse)

# input fasta (reference sequence must be first sequence)
fasta_file <- "file_name.fasta" 

# read the alignment
seqs <- readDNAStringSet(fasta_file)
seq_mat <- as.matrix(seqs)

# separate the reference from the query sequences (allows pairwise comparison)
ref <- as.character(seq_mat[1, ])
query_mat <- seq_mat[-1, , drop = FALSE]

mutation_types <- c("A", "C", "G", "T")

# build mutation (SNP) matrix where ('.' = MATCH)
mutation_matrix <- apply(query_mat, 1, function(row) {
  ifelse(row == ref, ".", as.character(row))
})
mutation_matrix <- t(mutation_matrix)
rownames(mutation_matrix) <- names(seqs)[-1]
colnames(mutation_matrix) <- seq_len(ncol(query_mat))


# determine the total base composition (including MATCHES)

# expand reference to query matrix shape
ref_mat <- matrix(
  ref,
  nrow = nrow(mutation_matrix),
  ncol = length(ref),
  byrow = TRUE
)

# resolve '.' to reference base
resolved_mat <- mutation_matrix
resolved_mat[resolved_mat == "."] <- ref_mat[resolved_mat == "."]

# count A/C/G/T (exclude gaps and Ns)
base_composition <- table(
  factor(
    resolved_mat[resolved_mat %in% mutation_types],
    levels = mutation_types
  )
)

base_composition_df <- tibble(
  metric = paste0("total_", names(base_composition)),
  count  = as.integer(base_composition)
)


# count direction mutations (once per event)
# for sanity check, A>G should not be equal to G>A in most cases
directional_counts <- matrix(
  0,
  nrow = length(mutation_types),
  ncol = length(mutation_types),
  dimnames = list(
    from = mutation_types,
    to   = mutation_types
  )
)

for (j in seq_len(ncol(query_mat))) {
  ref_base <- ref[j]
  query_bases <- query_mat[, j]
  
  valid <- query_bases != ref_base &
    query_bases %in% mutation_types &
    ref_base %in% mutation_types
  
  if (any(valid)) {
    for (qb in query_bases[valid]) {
      directional_counts[ref_base, qb] <-
        directional_counts[ref_base, qb] + 1
    }
  }
}

# convert directional counts to data frame (an extra precaution)
directional_df <- expand.grid(
  from = mutation_types,
  to   = mutation_types,
  stringsAsFactors = FALSE
) %>%
  mutate(count = as.vector(directional_counts)) %>%
  filter(from != to & count > 0) %>%
  mutate(mutation = paste0(from, ">", to))

# defining the symmetric mutation classes
# this defines the A<>G, C<>T, A<>C, A<>T, C<>G, G<>T as well as A>G, G>A, C>T, T>C, A>C, C>A, A>T, T>A, C>G, G>C, G>T, T>G
symmetric_df <- tibble(
  class = c("A<>G", "C<>T", "A<>C", "A<>T", "C<>G", "G<>T"),
  count = c(
    sum(directional_counts["A","G"], directional_counts["G","A"]),
    sum(directional_counts["C","T"], directional_counts["T","C"]),
    sum(directional_counts["A","C"], directional_counts["C","A"]),
    sum(directional_counts["A","T"], directional_counts["T","A"]),
    sum(directional_counts["C","G"], directional_counts["G","C"]),
    sum(directional_counts["G","T"], directional_counts["T","G"])
  )
)

# if interested in deaminase specific counts
# primarily if looking for APOBEC3 driven mutations

apobec_GR  <- 0   # GR  > AR  (permissive APOBEC3G/F)
apobec_NGR <- 0   # All possible APOBEC3 NGR=GR

# taking the time to explain this section in case other mutations are added. If the motif is NTA this would be j-1 (ref_prev), j (ref_mid), j+1 position (ref_next). qb == ref_mid indicates which base is undergoing the mutation, in this case the ref_mid [j]. If ambiguous base (e.g. R=G/A) is at a position, the possibilities should be listed as the "%in%" section. Ensure the language is consistent--most errors come from this.

for (j in 2:(ncol(query_mat)-1)) {
  
  ref_prev <- ref[j-1]
  ref_mid  <- ref[j]
  ref_next <- ref[j+1]
  
  if (!(ref_prev %in% mutation_types &&
        ref_mid  %in% mutation_types &&
        ref_next %in% mutation_types)) next
  
  query_mid <- query_mat[, j]
  
  valid <- query_mid %in% mutation_types
  if (!any(valid)) next
  
  for (qb in query_mid[valid]) {
    
    if (qb == ref_mid) next

    # --- GR > AR ---
    if (ref_mid == "G" && ref_next %in% c("A","G")) {
      if (qb == "A") {
        apobec_GR <- apobec_GR + 1
      }
    }
    
    # ---NGR > NAG/ NAA ---
    if (ref_prev %in% c("A","C","T","G") &&
        ref_mid == "G" &&
        ref_next %in% c("G","A")) {
      if(qb =="A") {
        apobec_NGR <- apobec_NGR + 1
      }
    }
  }
}

# these should all match the above section
deaminase_df <- tibble(
  mutation = c(
    "GR>AR (permissive APOBEC3G/F)",
    "NGR>NAR"
  ),
  count = c(
    apobec_GR,
    apobec_NGR
  )
)


# determining the number of deaminase mutations per sequence
seq_names <- rownames(query_mat)

context_table <- data.frame(
  sequence = seq_names,
  total_mut = 0,
  GR_to_AR = 0,
  NGR_to_NAR = 0,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(query_mat))) {
  
  seq <- query_mat[i, ]
  
  for (j in 2:(length(ref)-1)) {
    
    ref_prev <- ref[j-1]
    ref_mid  <- ref[j]
    ref_next <- ref[j+1]
    qb <- seq[j]
    
    if (!(qb %in% mutation_types &&
          ref_prev %in% mutation_types &&
          ref_mid %in% mutation_types &&
          ref_next %in% mutation_types)) next
    
    if (qb == ref_mid) next
    
    context_table$total_mut[i] <- context_table$total_mut[i] + 1
    
    # GR > AR
    if (ref_mid == "G" && ref_next %in% c("A","G") && qb == "A") {
      context_table$GR_to_AR[i] <- context_table$GR_to_AR[i] + 1
    }
    
    # NGR>NAR
    if (ref_prev %in% c("A","C","T","G") && ref_mid=="G" && ref_next %in% c("G","A") && qb == "A") {
      context_table$NGR_to_NAR[i] <- context_table$NGR_to_NAR[i] + 1
    }
    
  }
}
###### ~ the deaminase tangent is over ######

# count the total number of mutations
total_mutations <- sum(directional_counts)

# determine the transition vs transversion ratio (Ti/Tv)
transitions <- symmetric_df$count[symmetric_df$class %in% c("A<>G", "C<>T")]
transversions <- total_mutations - sum(transitions)
ti_tv_ratio <- ifelse(transversions > 0,
                      sum(transitions) / transversions,
                      NA)


# total mutation counts
base_counts <- colSums(directional_counts)

base_counts_df <- tibble(
  metric = paste0("count", names(base_counts)),
  count  = as.integer(base_counts)
)

# building the summary table
# make sure all of the elements you wish to include are present
sub_summary_extended <- bind_rows(
  directional_df %>% select(mutation, count),
  symmetric_df   %>% transmute(mutation = class, count),
  deaminase_df,
  tibble(mutation = "Total_mutations", count = total_mutations),
  tibble(mutation = "Ti/Tv", count = ti_tv_ratio),
  base_composition_df,
  base_counts_df
)


# writing the output csv files, change file names as needed

# this file will list each sequence (ref seq not shown), where '.' indicate a MATCH and alternative will show as A/C/T/G/-
write.csv(
  mutation_matrix,
  "file_name_matrix.csv",
  row.names = TRUE
)

# this file will have the extended table with the directional substitution counts, symmetric counts, deaminase, total mutations, Ti/Tv, base composition, etc
write.csv(
  sub_summary_extended,
  "file_name__matrix_extended.csv",
  row.names = FALSE
)

# this file will list each sequence and the number of total mutations and GR>AR and NGR>NAR
write.csv(
  context_table,
  "file_name_deaminase_context_per_sequence.csv",
  row.names = FALSE
)

### END ###
