###################################################
#         Highlighter plot with annotation        #
#                  Version 2.0.                   #
#         Written by Adam A. Capoferri, PhD       #
###################################################

# ---- LOAD LIBRARIES ----
library(tidyverse)
library(dplyr)
library(tibble)
library(stringr)

# ---- INPUT FILES ----
mutation_csv <- "name.csv" #change the name of the csv file
output_plot  <- "name.pdf" #change the name of the pdf file

# ---- DEFINE GENOME FEATURES ----
# change the start and end positions before running
features <- tribble( 
  ~label,   ~start, ~end,  ~feature_tier,
  "5' UTR",       1, 503,        3,
  "gag",       504, 1315,        3,
  "L",        2270, 2314,        3,
  "V",        1316, 1642,        2,
  "CDR3",     1643, 1684,        2,
  "J",        1685, 1717,        2,
  "C",        1718, 2269,        2,
  "TCRb",     1316, 2269,        3,
  "V",        2315, 2647,        2,
  "CDR3",     2648, 2695,        2,
  "J",        2696, 2728,        2,
  "C",        2729, 3136,        2,
  "TCRa",     2315, 3136,        3,
  "3' UTR",   3137, 3238,        3
)

# ---- READ MUTATION MATRIX ----
mutation_df <- read.csv(
  mutation_csv,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ---- REMOVE REFERENCE ROW ----
# the first sequence will be considered the reference/consensus
mutation_df <- mutation_df[!rownames(mutation_df) %in% c("reference", "Reference"), ]

# ---- IDENTIFY GENOMIC POSITION COLUMNS (EXCLUDE COUNTS) ----
genome_positions <- colnames(mutation_df) %>%
  str_subset("^[0-9]+$") %>%
  as.numeric()

stopifnot(length(genome_positions) > 0)

# ---- STORE SEQUENCE IDS ----
seq_ids <- rownames(mutation_df)

# ---- COUNT MUTATIONS PER SEQUENCE (FOR ORDERING) ----
mutation_counts <- mutation_df %>%
  rownames_to_column("seq_id") %>%
  pivot_longer(
    cols = matches("^[0-9]+$"),
    names_to = "position",
    values_to = "base"
  ) %>%
  filter(!is.na(base) & base != ".") %>%
  count(seq_id, name = "n_mutations")

# ---- ORDER SEQUENCES: LEAST → MOST MUTATIONS ----
seq_order <- tibble(seq_id = seq_ids) %>%
  left_join(mutation_counts, by = "seq_id") %>%
  mutate(n_mutations = replace_na(n_mutations, 0)) %>%
  arrange(n_mutations) %>%
  mutate(seq_index = rev(row_number()))

# ---- UPDATE SEQUENCE LINES ----
sequence_lines <- seq_order %>%
  transmute(
    seq_id,
    seq_index,
    y_start = seq_index,
    y_end   = seq_index
  )

# ---- FILTER MUTATIONS TO KEPT SEQUENCES ----
mutation_long <- mutation_long %>%
  filter(seq_id %in% seq_order$seq_id)

# ---- LONG-FORM MUTATION DATA ----
mutation_long <- mutation_df %>%
  rownames_to_column("seq_id") %>%
  pivot_longer(
    cols = matches("^[0-9]+$"),
    names_to = "position",
    values_to = "base"
  ) %>%
  mutate(position = as.numeric(position)) %>%
  filter(!is.na(base) & base != ".") %>%
  left_join(seq_order, by = "seq_id")

# ---- PRECOMPUTE X-JITTER FOR MUTATION TICKS ----
set.seed(123)
mutation_long <- mutation_long %>%
  mutate(jitter_x = position + runif(n(), -0.1, 0.1))

# ---- FEATURE Y POSITIONS (TOP OF PLOT) ----
feature_y_base <- max(sequence_lines$seq_index) + 1
features <- features %>%
  mutate(y = feature_y_base + feature_tier * 0.9) #incresaes the tier spacing
set.seed(123)

features <- features %>%
  mutate(
    label_y = if_else(
      feature_tier == 2,
      y - 1,   # drop below rectangle
      y
    ),
    label_x = (start + end) / 3 +
      if_else(feature_tier == 2,
              runif(n(), -3, 3),  # x-jitter for tier 2 only
              0)
  )

# ---- BUILD PLOT ----
p <- ggplot() +
  
  # Sequence baselines
  geom_segment(
    data = sequence_lines,
    aes(
      x = min(genome_positions),
      xend = max(genome_positions),
      y = seq_index,
      yend = seq_index
    ),
    color = "gray80",
    linewidth = 0.4
  ) +
  
  # Mutation ticks
  geom_segment(
    data = mutation_long,
    aes(
      x = jitter_x,
      xend = jitter_x,
      y = seq_index - 0.3,
      yend = seq_index + 0.3,
      color = base
    ),
    linewidth = 0.6
  ) +
  
  # Feature boxes
  geom_rect(
    data = features,
    aes(
      xmin = start,
      xmax = end,
      ymin = y - 0.35,
      ymax = y + 0.35
    ),
    fill = "grey85",
    color = "black"
  ) +
  
  # Feature labels (centered)
  geom_text(
    data = features,
    aes(
      x = (start + end) / 2,
      y = label_y,
      label = label
    ),
    size = 3,
    hjust = 0.5,
    vjust = 0.5
  ) +
  
  # Axes and scales
  scale_x_continuous(
    limits = c(min(genome_positions) - 1,
               max(genome_positions) + 1),
    breaks = pretty(genome_positions),
    name   = "Genomic Position (bp)"
  ) +
  
  scale_y_continuous(
    breaks = sequence_lines$seq_index,
    labels = sequence_lines$seq_id,
    expand = expansion(add = c(1, 3))
  ) +
  
  scale_color_manual(
    values = c(
      A = "red1",
      T = "green3",
      C = "blue2",
      G = "orange1",
      `-` = "gray20"
    )
  ) +
  guides(color = "none") +
  labs(
    x = "Genomic Position (bp)",
    y = NULL,
    title = ""
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# ---- SAVE PLOT ----
ggsave(
  filename = output_plot,
  plot     = p,
  width    = 12,
  height   = 8,
  units    = "in"
)
# End
