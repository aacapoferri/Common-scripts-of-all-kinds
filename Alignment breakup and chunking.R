# Breaking up a large alignment

# This is a bit niche, but if you have an alignment of many sequences and need to break it up into small chunks...maybe to run through another program, this will separate into multiple files. Basically, if you have 1000 sequences but a tool can only run 500 at a time, it will generate 2 chunks.


# load library
library(Biostrings) 

# Load the full FASTA file
all_seqs <- readDNAStringSet("file_name_sequences.fasta")

# Set chunk size
chunk_size <- 499
n_chunks <- ceiling(length(all_seqs) / chunk_size)  # Calculate number of chunks

# Get the base name from the FASTA file
base_name <- tools::file_path_sans_ext(basename("file_name_sequences_chunksize.fasta"))

# Split into chunks and write to separate files
for (i in seq_len(n_chunks)) {
  start <- (i - 1) * chunk_size + 1
  end <- min(i * chunk_size, length(all_seqs))
  chunk <- all_seqs[start:end]
  
  # Create the filename using the base name and chunk number
  output_filename <- paste0(base_name, "_chunk_", i, ".fasta")
  
  # Write the chunk to the file with the new name
  writeXStringSet(chunk, filepath = output_filename)
}

### END ###
