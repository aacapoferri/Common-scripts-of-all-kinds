###################################################
#         Highlighter plot with annotationn       #
#             Version 1.0 (June 2025)             #
#         Written by Adam A. Capoferri, PhD       #               
###################################################

# This is a generic highlighter plot script that can be used. The csv file that is used is generated from the "Nucleotide datamatrix and mutation frequency.R script. Here the first sequence (reference) is not shown but any collapsed sequences that have no differences to the reference will be the first line. The order of the sequences go from the #most collapsed to least. Nucleotides for A, C, T, G, (and gaps) are displayed. The mutation relative to the reference is shown in the highlighter plot. 

#When it comes to the annotations features, this can become a little tricky. The current set up is for TCR transgene proviruses. The order can be played around with a tiny bit but things do start to get wonky. The ~feature_tier is key so that the features are flying off somewhere. This is currently set so that the CDR3 feature is set below its annotation map so that it is visible. 

#There is a jitter to help ensure that mutations don't overlap, however if they are too close together the fine-tuning of jitter can't properly separate them without sending other elements off base. It is cumbersome, but some hand editing in Illustrator or similar might be necessary in some cases.

# As always, ensure all packages are installed first.

library(tidyverse)

# ---- MODIFY HERE: Input files ----
mutation_csv <- "PID X_alignment with master_collapsed_mutation_matrix.csv"
output_plot <- "PID X_alignment with master_collapsed_mutation_matrix_highlighter_plot.png"

# ---- DEFINE FEATURES WITH FIXED VISUAL STACKING ORDER ----
features <- tribble(
  ~label,     ~start, ~end,   ~feature_tier,
  "5' UTR",      1, 1316,        3,
  "Psi",        18,  504,        2,
  "gag",       505, 1316,        2,
  "L",        2271, 2315,        3,
  "V",        1317, 1643,        2,
  "CDR3",     1644, 1685,        2,
  "J",        1686, 1718,        2,
  "C",        1719, 2270,        2,
  "TCRb",     1317, 2270,        3,
  "V",        2316, 2648,        2,
  "CDR3",     2649, 2696,        2,
  "J",        2697, 2729,        2,
  "C",        2730, 3137,        2,
  "TCRa",     2316, 3137,        3,
  "3' UTR",   3137, 3275,        3
)

# ---- READ AND PREPARE MUTATION MATRIX ----
mutation_df <- read.csv(mutation_csv, row.names = 1, check.names = FALSE)
ref <- as.character(mutation_df["reference", ])
mutation_df <- mutation_df[!(rownames(mutation_df) %in% c("reference", "count_A", "count_T", "count_C", "count_G", "total_mut")), ]

# ---- ORDER SEQUENCES BY COLLAPSED COUNT (Descending) ----
seq_order <- rownames(mutation_df) %>%
  tibble(seq_id = .) %>%
  mutate(count = as.numeric(str_extract(seq_id, "(?<=\\()\\d+(?=\\))"))) %>%
  arrange(desc(count)) %>%
  mutate(seq_index = rev(row_number()))

# ---- CONVERT TO LONG FORMAT ----
mutation_long <- mutation_df %>%
  rownames_to_column("seq_id") %>%
  pivot_longer(-seq_id, names_to = "position", values_to = "base") %>%
  mutate(position = as.integer(position)) %>%
  filter(!is.na(base) & base != ".") %>%
  left_join(seq_order, by = "seq_id")

# ---- FORMAT MUTATION LABELS ----
ref_vec <- as.character(ref)
names(ref_vec) <- colnames(mutation_df)

mutation_long <- mutation_long %>%
  mutate(
    type = ifelse(base %in% c("A", "T", "C", "G"), "mutation", "gap_or_N"),
    seq_id = factor(seq_id, levels = seq_order$seq_id),
    base_label = ifelse(base %in% c("A", "T", "C", "G"), base, "-"),
    ref_base = ref_vec[as.character(position)],
    label = ifelse(base_label == "-", "-", paste0(ref_base, ">", base_label))
  )

n_seq <- length(seq_order$seq_id)

# ---- CALCULATE Y POSITION BASED ON FIXED VISUAL TIERS ----
features <- features %>%
  mutate(
    ymin = n_seq + 1.5 + feature_tier,
    ymax = n_seq + 2.5 + feature_tier,
    ytext = n_seq + 2 + feature_tier
  )

# ---- JITTER MUTATION POSITIONS ----
set.seed(42)
mutation_long <- mutation_long %>%
  mutate(jitter_x = position + runif(n(), -0.2, 0.2))

# ---- ALL SEQUENCE LINES ----
sequence_lines <- seq_order %>% select(seq_id, seq_index)

# ---- PLOT ----
p <- ggplot() +
  geom_segment(data = sequence_lines,
               aes(x = min(as.integer(names(ref_vec))),
                   xend = max(as.integer(names(ref_vec))),
                   y = seq_index, yend = seq_index),
               color = "gray80", linewidth = 0.4) +
  geom_segment(data = mutation_long,
               aes(x = jitter_x, xend = jitter_x,
                   y = seq_index - 0.3, yend = seq_index + 0.3,
                   color = base_label),
               linewidth = 0.6) +
  geom_text(data = mutation_long,
            aes(x = jitter_x, y = seq_index, label = label),
            color = "black", size = 2.5,
            position = position_nudge(y = 0.5), vjust = 0) +
  scale_y_continuous(breaks = sequence_lines$seq_index,
                  labels = sequence_lines$seq_id,
                  expand = expansion(add = 1)) +
  scale_color_manual(values = c(
    A = "blue", T = "green3", C = "orange", G = "red", `-` = "gray30"
  )) +
  guides(color = "none") +
  labs(x = "Position (bp)", y = NULL, title = "PID 201403 Proviral structures") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# ---- ADD FEATURE BOXES ----
if (nrow(features) > 0) {
  p <- p +
    geom_rect(data = features,
              aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE,
              fill = "white", color = "black", alpha = 0.4) +
    geom_text(
      data = features,
      aes(
        x = (start + end)/2,
        y = ytext,
        label = label,
        hjust = 0.5,
        vjust=ifelse(label=="CDR3",3.5,0.5)
      ),
      inherit.aes = FALSE,
      size = 3.5
    )
}

# ---- SAVE PLOT ----
ggsave(output_plot, p, width = 14, height = 6 + max(features$feature_tier), dpi = 300)
cat("\u2714 Saved plot:", output_plot, "\n")
### End ###
