################################################
# SHAP with Bootstrapping XGBoost Evaluation   #
#                 Version 1.0                  #
#       Author: Adam A. Capoferri, PhD         #
#      Contact: adam.capoferri@nih.gov         #
################################################

# Shapley Addiditive explanation (SHAP) values provide more consistent, theoretically grounded measure of feature contribution at both gloabl and sample levels for a final XGBoost classifier. Based on the inteded use, the preprocessed training predictors, after recipe steps (factor conversion, imputation, dummy encoding, and SMOTE) are extracted and converted to a numeric matrix. SHAP values computed for each training sample and feature using xgboost predict function.

# Configuration
rdata_path <- "filename_xgb_fit2.RData" #change name as needed
out_prefix <- "fig_xgb_shap"
top_n <- 15 #change as needed
n_boot <- 1000 #bootstrap replicates, change as needed
set.seed(1234)

# load libraries
suppressPackageStartupMessages({
  library(tidymodels)
  library(xgboost)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(stringr)
})

# confirming order of features
feature_groups <- tibble(
  FeatureIndex = colnames(X_train),
  BioFeature = case_when(
    FeatureIndex == "net_charge" ~ "net_charge",
    str_starts(FeatureIndex, "glyco_status_") ~ "glyco_status",
    FeatureIndex == "features_x4_TRUE." ~ "features_x4",
    str_starts(FeatureIndex, "pos25_") ~ "pos25",
    TRUE ~ "other"
  )
)


## load data and model
objs <- load(rdata_path)

find_by_class <- function(cls) {
  for (o in objs) {
    obj <- get(o, envir = .GlobalEnv)
    if (inherits(obj, cls)) return(obj)
  }
  NULL
}

wf <- find_by_class("workflow")
if (is.null(wf)) stop("Workflow not found")

xgb_fit <- workflows::extract_fit_engine(wf)

# recover training matrix 
mold <- workflows::extract_mold(wf)

X_train <- mold$predictors %>% as.matrix()

# compute SHAP values
shap_raw <- predict(
  xgb_fit,
  newdata = X_train,
  predcontrib = TRUE
)

if (!is.list(shap_raw)) {
  shap_raw <- list(Class1 = shap_raw)
}

shap_long <- imap_dfr(shap_raw, function(mat, cls) {
  mat <- as_tibble(mat)
  
  if ("BIAS" %in% colnames(mat)) {
    mat <- select(mat, -BIAS)
  }
  
  mat %>%
    mutate(row = row_number()) %>%
    pivot_longer(
      cols = -row,
      names_to = "FeatureIndex",
      values_to = "SHAP"
    ) %>%
    left_join(feature_groups, by = "FeatureIndex") %>%
    mutate(Class = cls)
})

boot_ci <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(tibble(mean = 0, lo = 0, hi = 0))
  }
  
  boots <- replicate(
    n_boot,
    mean(abs(sample(x, replace = TRUE)))
  )
  
  tibble(
    mean = mean(boots),
    lo = quantile(boots, 0.025, na.rm = TRUE),
    hi = quantile(boots, 0.975, na.rm = TRUE)
  )
}

shap_ci <- shap_long %>%
  group_by(Class, BioFeature) %>%
  summarise(boot_ci(SHAP), .groups = "drop")


top_feats <- shap_ci %>%
  filter(is.finite(mean), mean > 0) %>%
  group_by(Class) %>%
  arrange(desc(mean), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup()

# bootstrap confidence intervals
boot_ci <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(tibble(mean = 0, lo = 0, hi = 0))
  }
  
  boots <- replicate(
    n_boot,
    mean(abs(sample(x, replace = TRUE)))
  )
  
  tibble(
    mean = mean(boots),
    lo   = quantile(boots, 0.025, na.rm = TRUE),
    hi   = quantile(boots, 0.975, na.rm = TRUE)
  )
}


# top features per class
top_feats <- shap_ci %>%
  filter(is.finite(mean), mean > 0) %>%
  group_by(Class) %>%
  arrange(desc(mean), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup()

if (nrow(top_feats) == 0) {
  stop("No SHAP features available to plot")
}

# generating plot
p <- ggplot(top_feats, aes(x = reorder(BioFeature, mean), y = mean)) +
  geom_col(fill = "#2563EB") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25) +
  coord_flip() +
  facet_wrap(~ Class, scales = "free_y") +
  labs(
    title = "Per-class SHAP feature importance (aggregated biologically)",
    y = "Mean |SHAP| (bootstrap 95% CI)",
    x = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "bold"))

ggsave(
  paste0(out_prefix, "_shap.png"),
  p,
  width = 9,
  height = 6,
  dpi = 300
)


# export table as csv
# helpful if wanting to generate figure outside
write.csv(
  shap_ci,
  paste0(out_prefix, "_shap_table.csv"),
  row.names = FALSE
)

cat("✓ SHAP analysis complete\n")

### END ###
