## ============================================================
## Hydra AEP – Genome Feature Modeling
## Response variables:
##   1. is_duplicate   → Logistic regression (binomial GLM)
##   2. dist_near_gene → Gaussian GLM (log-transformed)
##   3. gene_density   → Negative Binomial GLM (count)
##
## Predictors: gene_length, scaffold_chr_length, is_opsin,
##             strand, (+ cross-predictors where appropriate)
##
## Outputs:
##   - Coefficient tables with odds ratios / rate ratios
##   - Correlation / feature importance summary table
##   - Plots: coefficient plots, correlation heatmap
## ============================================================

## ── 0. Packages ───────────────────────────────────────────────
pkgs_cran <- c("dplyr", "ggplot2", "MASS", "broom", "corrplot",
               "patchwork", "scales", "tibble", "stringr", "readr")
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(dplyr)
library(ggplot2)
library(MASS)        # glm.nb
library(broom)
library(corrplot)
library(patchwork)
library(scales)
library(tibble)
library(stringr)
library(readr)

## ── 1. Load data ──────────────────────────────────────────────
# Expects hydra_aep_gene_table.csv from the companion data-building script.
# If you have the RDS, swap to: gene_table <- readRDS("hydra_aep_gene_table.rds")
if (file.exists("hydra_aep_gene_table.rds")) {
  gene_table <- readRDS("hydra_aep_gene_table.rds")
} else if (file.exists("hydra_aep_gene_table.csv")) {
  gene_table <- read_csv("hydra_aep_gene_table.csv", show_col_types = FALSE)
} else {
  stop("Run hydra_aep_gene_table.R first to generate the gene table.")
}

## ── 2. Pre-processing ─────────────────────────────────────────
model_df <- gene_table %>%
  mutate(
    is_duplicate        = as.integer(is_duplicate),          # 0/1
    is_opsin            = as.integer(is_opsin),              # 0/1
    strand_plus         = as.integer(strand == "+"),         # 0/1
    log_gene_length     = log1p(gene_length),
    log_scaffold_length = log1p(scaffold_chr_length),
    log_dist_near_gene  = log1p(dist_near_gene),
    gene_density        = as.integer(gene_density)
  ) %>%
  filter(
    !is.na(dist_near_gene),
    !is.na(gene_density),
    !is.na(scaffold_chr_length)
  )

cat("Modelling dataset:", nrow(model_df), "genes\n")
cat("Duplicate genes:  ", sum(model_df$is_duplicate), "\n")
cat("Opsin genes:      ", sum(model_df$is_opsin), "\n\n")

## ── 3. Predictor set ──────────────────────────────────────────
# Core predictors used across all models
PREDICTORS <- c("log_gene_length", "log_scaffold_length",
                "is_opsin", "strand_plus")

formula_base <- function(response) {
  as.formula(paste(response, "~", paste(PREDICTORS, collapse = " + ")))
}

## ═══════════════════════════════════════════════════════════════
## MODEL 1: is_duplicate ~ predictors  (Logistic regression)
## ═══════════════════════════════════════════════════════════════
cat("── Model 1: Logistic regression for is_duplicate ──\n")

m1 <- glm(formula_base("is_duplicate"),
          data   = model_df,
          family = binomial(link = "logit"))

m1_tidy <- tidy(m1, conf.int = TRUE, exponentiate = TRUE) %>%
  rename(odds_ratio = estimate, OR_lo = conf.low, OR_hi = conf.high) %>%
  mutate(model = "Logistic (is_duplicate)", significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    p.value < 0.1   ~ ".",
    TRUE            ~ ""
  ))

cat("\nOdds Ratios (exponentiated coefficients):\n")
print(m1_tidy %>% select(term, odds_ratio, OR_lo, OR_hi, p.value, significance))

cat("\nModel 1 AIC:", AIC(m1), "\n")
cat("Null deviance:", m1$null.deviance,
    "| Residual deviance:", m1$deviance,
    "| McFadden R²:", round(1 - m1$deviance / m1$null.deviance, 4), "\n\n")

## ═══════════════════════════════════════════════════════════════
## MODEL 2: dist_near_gene ~ predictors  (Gaussian GLM, log link)
## ═══════════════════════════════════════════════════════════════
cat("── Model 2: Gaussian GLM for dist_near_gene (log-transformed) ──\n")

# We model log1p(dist_near_gene) as response with Gaussian/identity,
# equivalent to a log-linear model; easy to interpret coefficients.
m2 <- glm(log_dist_near_gene ~ log_gene_length + log_scaffold_length +
            is_opsin + strand_plus + is_duplicate,
          data   = model_df,
          family = gaussian(link = "identity"))

m2_tidy <- tidy(m2, conf.int = TRUE) %>%
  mutate(model = "Gaussian GLM (log dist_near_gene)", significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    p.value < 0.1   ~ ".",
    TRUE            ~ ""
  ))

cat("\nCoefficients (additive effect on log distance):\n")
print(m2_tidy %>% select(term, estimate, conf.low, conf.high, p.value, significance))
cat("\nModel 2 AIC:", AIC(m2), "\n\n")

## ═══════════════════════════════════════════════════════════════
## MODEL 3: gene_density ~ predictors  (Negative Binomial GLM)
## ═══════════════════════════════════════════════════════════════
cat("── Model 3: Negative Binomial GLM for gene_density ──\n")

# Poisson first to check for overdispersion
m3_pois <- glm(gene_density ~ log_gene_length + log_scaffold_length +
                 is_opsin + strand_plus + is_duplicate,
               data   = model_df,
               family = poisson(link = "log"))

dispersion_ratio <- m3_pois$deviance / m3_pois$df.residual
cat("Poisson dispersion ratio (deviance/df):", round(dispersion_ratio, 2),
    "–", ifelse(dispersion_ratio > 1.5, "overdispersed → use NegBin", "acceptable"), "\n")

# Fit Negative Binomial regardless (robust to both cases)
m3 <- glm.nb(gene_density ~ log_gene_length + log_scaffold_length +
               is_opsin + strand_plus + is_duplicate,
             data = model_df)

m3_tidy <- tidy(m3, conf.int = TRUE, exponentiate = TRUE) %>%
  rename(rate_ratio = estimate, RR_lo = conf.low, RR_hi = conf.high) %>%
  mutate(model = "Neg. Binomial (gene_density)", significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    p.value < 0.1   ~ ".",
    TRUE            ~ ""
  ))

cat("\nRate Ratios (exponentiated coefficients):\n")
print(m3_tidy %>% select(term, rate_ratio, RR_lo, RR_hi, p.value, significance))
cat("\nModel 3 AIC:", AIC(m3), "| Theta:", round(m3$theta, 4), "\n\n")

## ═══════════════════════════════════════════════════════════════
## FEATURE IMPORTANCE SUMMARY TABLE
## ═══════════════════════════════════════════════════════════════
cat("── Feature Importance Summary ──\n")

# Standardised coefficients: re-fit each model on z-scored predictors
# so effect sizes are directly comparable across predictors.
scale_predictors <- c("log_gene_length", "log_scaffold_length",
                      "is_opsin", "strand_plus", "is_duplicate")

model_df_scaled <- model_df %>%
  mutate(across(all_of(scale_predictors), scale))

# M1 scaled (no is_duplicate as predictor here)
m1s <- glm(is_duplicate ~ log_gene_length + log_scaffold_length +
             is_opsin + strand_plus,
           data = model_df_scaled, family = binomial())

m2s <- glm(log_dist_near_gene ~ log_gene_length + log_scaffold_length +
             is_opsin + strand_plus + is_duplicate,
           data = model_df_scaled, family = gaussian())

m3s <- glm.nb(gene_density ~ log_gene_length + log_scaffold_length +
                is_opsin + strand_plus + is_duplicate,
              data = model_df_scaled)

extract_importance <- function(model, response_label) {
  tidy(model, conf.int = FALSE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      response       = response_label,
      abs_std_effect = abs(estimate),
      direction      = ifelse(estimate > 0, "positive", "negative")
    ) %>%
    select(response, predictor = term, std_estimate = estimate,
           abs_std_effect, direction, p.value)
}

importance_table <- bind_rows(
  extract_importance(m1s, "is_duplicate"),
  extract_importance(m2s, "log(dist_near_gene)"),
  extract_importance(m3s, "gene_density")
) %>%
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    p.value < 0.1   ~ ".",
    TRUE            ~ ""
  )) %>%
  arrange(response, desc(abs_std_effect))

cat("\nStandardised effect sizes (predictors scaled to mean=0, SD=1):\n")
print(importance_table)

## ═══════════════════════════════════════════════════════════════
## CORRELATION MATRIX (numeric features)
## ═══════════════════════════════════════════════════════════════
corr_vars <- model_df %>%
  select(gene_length, scaffold_chr_length, is_opsin, strand_plus,
         is_duplicate, dist_near_gene, gene_density) %>%
  rename(
    gene_len    = gene_length,
    scaff_len   = scaffold_chr_length,
    opsin       = is_opsin,
    strand_p    = strand_plus,
    duplicate   = is_duplicate,
    dist_near   = dist_near_gene,
    gene_dens   = gene_density
  )

corr_matrix <- cor(corr_vars, use = "pairwise.complete.obs", method = "spearman")

## ═══════════════════════════════════════════════════════════════
## PLOTS
## ═══════════════════════════════════════════════════════════════

theme_set(theme_minimal(base_size = 12))

## — Plot 1: Coefficient plot – Model 1 (OR) ────────────────────
p1 <- m1_tidy %>%
  filter(term != "(Intercept)") %>%
  mutate(term = str_replace_all(term, c(
    "log_gene_length"     = "Gene length (log)",
    "log_scaffold_length" = "Scaffold length (log)",
    "is_opsin"            = "Is opsin",
    "strand_plus"         = "Strand (+)"
  ))) %>%
  ggplot(aes(x = odds_ratio, y = reorder(term, odds_ratio))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0.2, colour = "#2c7bb6") +
  geom_point(size = 3, colour = "#2c7bb6") +
  scale_x_log10(labels = label_number(accuracy = 0.01)) +
  labs(title = "Model 1: Predictors of gene duplication",
       subtitle = "Logistic regression – Odds Ratios (95% CI, log scale)",
       x = "Odds Ratio", y = NULL) +
  theme(plot.title = element_text(face = "bold"))

## — Plot 2: Standardised effect sizes across all models ─────────
p2 <- importance_table %>%
  mutate(
    predictor = str_replace_all(predictor, c(
      "log_gene_length"     = "Gene length (log)",
      "log_scaffold_length" = "Scaffold length (log)",
      "is_opsin"            = "Is opsin",
      "strand_plus"         = "Strand (+)",
      "is_duplicate"        = "Is duplicate"
    )),
    response = factor(response,
                      levels = c("is_duplicate",
                                 "log(dist_near_gene)",
                                 "gene_density"))
  ) %>%
  ggplot(aes(x = std_estimate, y = reorder(predictor, abs_std_effect),
             fill = direction)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = significance,
                x = std_estimate + sign(std_estimate) * 0.01),
            hjust = 0.5, size = 3.5) +
  facet_wrap(~response, scales = "free_x") +
  scale_fill_manual(values = c("positive" = "#d7191c", "negative" = "#2c7bb6")) +
  labs(title = "Feature importance: standardised coefficients",
       subtitle = "Predictors scaled to SD=1; effect sizes comparable within each model",
       x = "Standardised coefficient", y = NULL, fill = "Direction") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

## — Plot 3: Correlation heatmap ─────────────────────────────────
png("hydra_aep_correlation_heatmap.png", width = 700, height = 650, res = 110)
corrplot(corr_matrix,
         method   = "color",
         type     = "upper",
         order    = "hclust",
         addCoef.col = "black",
         number.cex  = 0.75,
         tl.col   = "black",
         tl.srt   = 45,
         col      = colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(200),
         title    = "Spearman correlations – Hydra AEP genome features",
         mar      = c(0, 0, 2, 0))
dev.off()
cat("Saved: hydra_aep_correlation_heatmap.png\n")

## — Combine coefficient + importance plots ──────────────────────
combined_plot <- p1 / p2 + plot_layout(heights = c(1, 1.4))

ggsave("hydra_aep_model_plots.png", combined_plot,
       width = 10, height = 10, dpi = 150)
cat("Saved: hydra_aep_model_plots.png\n")

## ═══════════════════════════════════════════════════════════════
## SAVE OUTPUTS
## ═══════════════════════════════════════════════════════════════
write_csv(m1_tidy,        "model1_logistic_is_duplicate.csv")
write_csv(m2_tidy,        "model2_gaussian_dist_near_gene.csv")
write_csv(m3_tidy,        "model3_negbin_gene_density.csv")
write_csv(importance_table, "feature_importance_summary.csv")

cat("\n── All outputs saved ──\n")
cat("  model1_logistic_is_duplicate.csv\n")
cat("  model2_gaussian_dist_near_gene.csv\n")
cat("  model3_negbin_gene_density.csv\n")
cat("  feature_importance_summary.csv\n")
cat("  hydra_aep_model_plots.png\n")
cat("  hydra_aep_correlation_heatmap.png\n")