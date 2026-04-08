###################################################
#              Collapsing sequences               #
#             Version 1.0 (June 2025)             #
#           Written by Adam A. Capoferri, PhD     #
#    Contact information: adam.capoferri@nih.gov  #
###################################################

# Collapse .fasta sequences where you are able to change the name and isolate 1 to n where the number of sequences collapsed is transferred to "()". So PID A isolate (10) indicates 10 sequences were collapsed into one.

##### PART I: Version for nucleotide based alignments
# Load library
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
isolate_labels <- paste0("PID X isolate ", isolate_ids, " (", counts, ")")

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







##### PART II: Version for amino acid based alignments

# load library
library(Biostrings)

# Function to collapse sequences
collapse_sequences <- function(input_fasta, output_fasta_collapsed) {
  # Read the input FASTA file
  alignment <- readAAStringSet(input_fasta)
  
  # Convert alignment to character vector
  sequences <- as.character(alignment)
  
  # Create a table to count occurrences of each unique sequence
  sequence_counts <- table(sequences)
  
  # Prepare vectors for collapsed names and sequences
  collapsed_names <- character(length(sequence_counts))
  collapsed_sequences <- character(length(sequence_counts))
  
  # Collapse sequences and modify sequence names
  i <- 1
  for (seq in names(sequence_counts)) {
    # Extract an original name for the sequence and append the count
    original_name <- names(alignment)[which(sequences == seq)[1]]
    collapsed_names[i] <- paste0(original_name, "-", sequence_counts[seq])
    collapsed_sequences[i] <- seq
    i <- i + 1
  }
  
  # Create the output alignment with collapsed sequences
  collapsed_alignment <- AAStringSet(collapsed_sequences)
  names(collapsed_alignment) <- collapsed_names
  
  # Write the collapsed alignment to an output file
  writeXStringSet(collapsed_alignment, filepath = output_fasta_collapsed, format = "fasta")
}

# Specify the input and output file paths
input_fasta <- "file_name.fasta"
output_fasta_collapsed <- "file_name_collapsed.fasta"

# Run the function to collapse sequences
collapse_sequences(input_fasta, output_fasta_collapsed)


### End ###
###########
