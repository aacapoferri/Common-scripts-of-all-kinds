####################################################
# Mutation freq, logo plot, substitution frequency #
#                 Version 1.0                      #
#       Author: Adam A. Capoferri, PhD             #
#      Contact: adam.capoferri@nih.gov             #
####################################################

# This was initially designed for analyzing HIV-1 V3 loop mutations. It is currently written to input all fasta (amino acid) files from a given folder, loop, and export figures as pdf back into the same folder. The first sequence in the fasta file must be your "Reference". 
# Figure is 3 tiered:
  # Top: Mutation Frequency vs. Reference. Compared to the Reference, what is the frequency at each position that there is a mutation (blue bars).
  # Middle: A Logo plot of the Reference sequence.
  # Bottom: The y-axis is scaled to 100 proportion. But at each given position where there is a mutation relative to the Reference, what are the amino acid changes present with their frequency. If at position 1, 0.1% of sequences do not match the Reference, 0.5 are alanine, 0.3 are glycine, 0.2 are valine

# Load libraries
library(Biostrings)
library(ggplot2)
library(ggseqlogo)
library(gridExtra)
library(grid)
library(cowplot)

# input; select file_path
fasta_files <- list.files(path = "file_path", pattern = "\\.fasta$", full.names = TRUE)

# defining amino acids into class and color
aa_chars <- c("A", "V", "I", "L", "M", "F", "W", "P", "G",  # Nonpolar
              "S", "T", "Y", "C", "N", "Q",                 # Polar
              "D", "E",                                     # Acidic (negatively charged)
              "R", "K", "H")                                # Basic (positively charged)

aa_groups <- c(rep("Nonpolar", 9), rep("Polar", 6), rep("Acidic", 2), rep("Basic", 3))

aa_colors <- c(rep("#000000", 9),  # Nonpolar, black
               rep("#2ca02c", 6),  # Polar, green
               rep("#1f77b4", 2),  # Acidic, blue
               rep("#d62728", 3))  # Basic, red

custom_col_scheme <- make_col_scheme(
  chars = aa_chars,
  groups = aa_groups,
  cols = aa_colors
)

# processing each fasta file
for (fasta_file in fasta_files) {
  
  output_pdf <- paste0(sub(".fasta", "", basename(fasta_file)), "_mutfreq_seqlogo_plot_final.pdf") # change output file name as needed
  
  alignment <- readAAStringSet(fasta_file)
  if (length(alignment) < 2) next
  
  alignment_matrix <- do.call(rbind, lapply(as.character(alignment), function(x) strsplit(x, "")[[1]]))
  
  ref_seq <- alignment_matrix[1, ]
  positions <- 1:ncol(alignment_matrix)
  
  # mutation frequency plot
  mutation_freq <- colSums(alignment_matrix != matrix(rep(ref_seq, each = nrow(alignment_matrix)), ncol = ncol(alignment_matrix))) / nrow(alignment_matrix)
  mutation_df <- data.frame(Position = positions, MutationFreq = mutation_freq)
  
  mutation_plot <- ggplot(mutation_df, aes(x = Position - 0.25, y = MutationFreq)) +
    geom_col(fill = "steelblue", width = 0.5) +
    scale_x_continuous(limits = c(0.5, length(positions) + 0.5), breaks = positions) +
    scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.02), expand = c(0, 0)) +
    labs(y = "Mutation Frequency vs Reference (%)", x = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title.y = element_text(size = 11, color = "black", margin = margin(r = 30)),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      plot.margin = unit(c(1, 0.6, 0, 0.6), "lines")
    )
  
  # reference logo plot
  ref_logo_input <- as.character(alignment[1]) # here is defining the 1st sequence in fasta as reference
  
  ref_logo_plot <- ggseqlogo(ref_logo_input, seq_type = "aa", method = "prob",
                             col_scheme = custom_col_scheme) +
    labs(y = "Reference") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_text(size = 11, color = "black", angle = 0, hjust = 0.5, margin = margin(r = 15)),
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black"),
      plot.margin = unit(c(1.8, 0.6, 0.4, 0.6), "lines")
    )
  
  # substitution freqency relative to total mutant plot
  mutation_only_matrix <- alignment_matrix[-1, , drop = FALSE]
  for (i in 1:nrow(mutation_only_matrix)) {
    mutation_only_matrix[i, mutation_only_matrix[i, ] == ref_seq] <- "-"
  }
  
  mutation_only_input <- apply(mutation_only_matrix, 1, paste0, collapse = "")
  
  seq_logo_plot <- ggseqlogo(mutation_only_input, seq_type = "aa", method = "prob",
                             col_scheme = custom_col_scheme) +
    labs(y = "Substitution Frequency\n(relative to total mutant)") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 11, color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title.y = element_text(size = 11, color = "black", margin = margin(r = 15)),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      plot.margin = unit(c(0.2, 0.6, 0.1, 0.6), "lines")
    )
  
  # combining the plots
  aligned_plots <- plot_grid(mutation_plot, ref_logo_plot, seq_logo_plot,
                             ncol = 1, align = "v", axis = "lr",
                             rel_heights = c(0.3, 0.3, 0.3))  # Adjust height ratios
  
  # save plots as pdf
  pdf(output_pdf, width = 12, height = 8)
  print(aligned_plots)
  dev.off()
  
  cat("Plot saved to", output_pdf, "\n")
}

### END ###
