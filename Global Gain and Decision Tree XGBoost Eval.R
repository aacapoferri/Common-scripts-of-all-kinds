####################################################
# Global Gain and Decision Tree XGBoost Evaluation #
#                 Version 1.0                      #
#       Author: Adam A. Capoferri, PhD             #
#      Contact: adam.capoferri@nih.gov             #
####################################################

# This script is meant to evaluate a machine-learning prediciton model. It was initially created to for the CRF01AE tropism prediction but should be applicable for any XGBoost visualization. 

# Publication-ready XGBoost visualization (tidymodels)
# Outputs:
# fig_xgb_importance.png
# fig_xgb_tree.png
# fig_xgb_combined.png
# feature importance table (includes Gain, Cover, and Frequency)


# Configuration
rdata_path        <- "filename_xgb_fit.RData" #change name as needed
out_prefix        <- "fig_xgb"
top_n_importance  <- 15
plot_width_px     <- 1600
plot_height_px    <- 1200


# Load libraries
suppressPackageStartupMessages({
  library(tidymodels)
  library(xgboost)
  library(dplyr)
  library(ggplot2)
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)
  library(magick)
})

# Load RData file and helpers
if (!file.exists(rdata_path)) {
  stop("RData file not found: ", rdata_path)
}

loaded_objs <- load(rdata_path)

find_by_class <- function(class_name) {
  for (nm in loaded_objs) {
    obj <- get(nm, envir = .GlobalEnv)
    if (inherits(obj, class_name)) return(obj)
  }
  NULL
}

wf  <- find_by_class("workflow")
rec <- find_by_class("recipe")

if (!is.null(wf) && is.null(rec)) {
  rec <- tryCatch(workflows::extract_recipe(wf), error = function(e) NULL)
}

train_df <- NULL
for (nm in loaded_objs) {
  obj <- get(nm, envir = .GlobalEnv)
  if (is.data.frame(obj)) {
    train_df <- obj
    break
  }
}


# Extract xgboost model
xgb_fit <- NULL

if (!is.null(wf)) {
  xgb_fit <- tryCatch(
    workflows::extract_fit_engine(wf),
    error = function(e) NULL
  )
}

if (is.null(xgb_fit)) {
  for (nm in loaded_objs) {
    obj <- get(nm, envir = .GlobalEnv)
    if (inherits(obj, "xgb.Booster")) {
      xgb_fit <- obj
      break
    }
  }
}

if (is.null(xgb_fit)) {
  stop("No xgboost model found in RData file.")
}


# Recover feature names
feat_names <- xgb_fit$feature_names

if (is.null(feat_names) && !is.null(rec) && !is.null(train_df)) {
  prep_rec <- tryCatch(recipes::prep(rec), error = function(e) NULL)
  if (!is.null(prep_rec)) {
    X <- tryCatch(
      recipes::bake(prep_rec, new_data = train_df, composition = "matrix"),
      error = function(e) NULL
    )
    if (!is.null(X)) feat_names <- colnames(X)
  }
}


# Feature importance ("Panel A")
imp <- xgboost::xgb.importance(
  model = xgb_fit,
  feature_names = feat_names
)

if (nrow(imp) == 0) stop("No importance values returned.")

imp_tbl <- imp %>%
  as_tibble() %>%
  arrange(desc(Gain)) %>%
  mutate(Feature = factor(Feature, levels = rev(Feature)))


# Export feature importance table as CSV 
# Included this here so that one may make figure outside of R
write.csv(
  imp_tbl %>% arrange(desc(Gain)),
  paste0(out_prefix, "_importance.csv"),
  row.names = FALSE
)

top_n <- min(top_n_importance, nrow(imp_tbl))

p_imp <- ggplot(
  imp_tbl %>% slice_head(n = top_n),
  aes(x = Feature, y = Gain)
) +
  geom_col(fill = "#367BF5") +
  coord_flip() +
  labs(
    title = "Global feature importance (Gain)",
    x = NULL,
    y = "Gain"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# Saving figure image, change file name as needed
imp_png <- paste0(out_prefix, "_importance.png")
ggsave(imp_png, p_imp, width = 7, height = 5, dpi = 300)


# Representative tree ("Panel B") For this, Tree zero will be generated
trees_dt <- xgboost::xgb.model.dt.tree(
  model = xgb_fit,
  feature_names = feat_names
)

trees_gain <- trees_dt %>%
  group_by(Tree) %>%
  summarize(total_gain = sum(Gain, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_gain))

best_tree <- trees_gain$Tree[1]

tree_graph <- xgboost::xgb.plot.tree(
  model = xgb_fit,
  feature_names = feat_names,
  n_first_tree = best_tree
)

tree_svg <- DiagrammeRsvg::export_svg(tree_graph)
tree_png <- paste0(out_prefix, "_tree.png")

rsvg::rsvg_png(
  svg = charToRaw(tree_svg),
  file = tree_png,
  width = plot_width_px,
  height = plot_height_px
)


# Combining panels (if interested)
img_imp  <- magick::image_read(imp_png)
img_tree <- magick::image_read(tree_png)

target_h <- max(
  magick::image_info(img_imp)$height,
  magick::image_info(img_tree)$height
)

img_imp  <- magick::image_resize(img_imp,  paste0("x", target_h))
img_tree <- magick::image_resize(img_tree, paste0("x", target_h))

combined <- magick::image_append(c(img_imp, img_tree))
combined_out <- paste0(out_prefix, "_combined.png")

magick::image_write(combined, combined_out)


# Finishing steps
cat(
  "Saved figures:\n",
  " -", imp_png, "\n",
  " -", tree_png, "\n",
  " -", combined_out, "\n"
)

### END ###
