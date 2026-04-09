# Converting a tsv file to fasta
# fairly simple, and focuses on amino acid here


# Load required library
library(Biostrings)

# Read the TSV file
tsv_data <- read.delim("file_name.tsv", header = TRUE, stringsAsFactors = FALSE)

# Confirm the expected columns are present
if (!all(c("name", "sequence") %in% names(tsv_data))) {
  stop("TSV must have columns named 'name' and 'sequence'")
}

# Use AAStringSet for amino acid sequences
fasta <- AAStringSet(tsv_data$sequence)
names(fasta) <- tsv_data$name

# Write to FASTA
writeXStringSet(fasta, "file_name.fasta")

### END ###
