Collapsed identical amino acid sequences

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
