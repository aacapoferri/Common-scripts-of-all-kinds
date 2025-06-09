###################################################
#              Collapsing sequences               #
#             Version 1.0 (June 2025)             #
#           Written by Adam A. Capoferri, PhD     #               
###################################################

# Collapse .fasta sequences where you are able to change the name and isolate 1 to n where the number of sequences collapsed is transferred to "()". So PID A isolate (10) indicates 10 sequences were collapsed into one.

# As always, ensure all packages are installed first.

library(Biostrings)

# ---- Load FASTA ----
fasta_file <- "PID X_alignment without master.fasta" #switch with .fasta you want
sequences <- readDNAStringSet(fasta_file)

# ---- Collapse identical sequences ----
seq_strings <- as.character(sequences)
seq_df <- data.frame(name = names(sequences), seq = seq_strings, stringsAsFactors = FALSE)
collapsed <- split(seq_df$name, seq_df$seq)

# ---- Create isolate names with unique number and collapsed count ----
isolate_ids <- seq_along(collapsed)
counts <- lengths(collapsed)
isolate_labels <- paste0("201403 isolate ", isolate_ids, " (", counts, ")")

# ---- Create collapsed DNAStringSet ----
unique_seqs <- DNAStringSet(names(collapsed))
names(unique_seqs) <- isolate_labels

# ---- Write collapsed FASTA ----
writeXStringSet(unique_seqs, "PID X_alignment without master_collapsed.fasta")

# ---- Create and write mapping CSV ----
mapping_df <- data.frame(
  original_name = unlist(collapsed),
  isolate = rep(isolate_labels, times = counts),
  stringsAsFactors = FALSE
)

write.csv(mapping_df, "PID X_alignment without master_collapsed.csv", row.names = FALSE)
### End ###
