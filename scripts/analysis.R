# ============================================================
# NTHU AIS User Satisfaction Analysis — IMBA Thesis
# Author:  Dale John Baltazar
# Method:  PLS-SEM via seminr
# Seed:    2024 | Bootstrap: 5,000 resamples
# All analyses are fully reproducible with this script.
#
# Usage:   Open the .Rproj file, then source this script.
#          All paths are relative to the project root.
# ============================================================

# ── 1. Libraries ────────────────────────────────────────────
suppressPackageStartupMessages({
  library(psych)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(seminr)
  library(gridExtra)
})

# ── 2. Load Data ────────────────────────────────────────────
raw <- read.csv(
  "data/Evaluating User Satisfaction with NTHUs AIS.csv",
  stringsAsFactors = FALSE,
  check.names      = FALSE
)
cat("Rows:", nrow(raw), "| Cols:", ncol(raw), "\n")

# ── 3. Recode Likert & Build Construct Data Frames ──────────
recode_likert <- function(x) {
  x2 <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
  x2 <- trimws(tolower(x2))
  dplyr::case_when(
    x2 == "strongly agree"    ~ 5L,
    x2 == "agree"             ~ 4L,
    x2 == "neutral"           ~ 3L,
    x2 == "disagree"          ~ 2L,
    x2 == "strongly disagree" ~ 1L,
    TRUE ~ NA_integer_
  )
}

IQ    <- data.frame(lapply(raw[, 11:14], recode_likert))
SQ    <- data.frame(lapply(raw[, 15:18], recode_likert))
SERVQ <- data.frame(lapply(raw[, 19:21], recode_likert))
MUX   <- data.frame(lapply(raw[, 22:25], recode_likert))
CA    <- data.frame(lapply(raw[, 26:29], recode_likert))
US    <- data.frame(lapply(raw[, 30:33], recode_likert))

colnames(IQ)    <- paste0("IQ",    1:4)
colnames(SQ)    <- paste0("SQ",    1:4)
colnames(SERVQ) <- paste0("SERVQ", 1:3)
colnames(MUX)   <- paste0("MUX",   1:4)
colnames(CA)    <- paste0("CA",    1:4)
colnames(US)    <- paste0("US",    1:4)

cat("NA counts — IQ:", sum(is.na(IQ)),
    "SQ:", sum(is.na(SQ)),
    "SERVQ:", sum(is.na(SERVQ)),
    "MUX:", sum(is.na(MUX)),
    "CA:", sum(is.na(CA)),
    "US:", sum(is.na(US)), "\n")

# ── 4. Respondent Profile ──────────────────────────────────
gender_clean <- trimws(iconv(raw[, 8],  to = "ASCII//TRANSLIT", sub = ""))
age_clean    <- trimws(iconv(raw[, 7],  to = "ASCII//TRANSLIT", sub = ""))
lev_raw      <- trimws(iconv(raw[, 9],  to = "ASCII//TRANSLIT", sub = ""))
freq_clean   <- trimws(iconv(raw[, 10], to = "ASCII//TRANSLIT", sub = ""))

lev_simple <- dplyr::case_when(
  grepl("PhD",      lev_raw, ignore.case = TRUE) ~ "PhD",
  grepl("Master",   lev_raw, ignore.case = TRUE) ~ "Master's",
  grepl("ndergrad", lev_raw, ignore.case = TRUE) ~ "Undergraduate",
  grepl("Exchange", lev_raw, ignore.case = TRUE) ~ "Exchange",
  TRUE ~ "Other"
)

freq_std <- dplyr::case_when(
  grepl("daily",   freq_clean, ignore.case = TRUE) ~ "Daily",
  grepl("weekly",  freq_clean, ignore.case = TRUE) ~ "Weekly",
  grepl("monthly", freq_clean, ignore.case = TRUE) ~ "Monthly",
  grepl("rarely",  freq_clean, ignore.case = TRUE) ~ "Rarely",
  TRUE ~ "Other"
)

cat("Gender:\n");         print(table(gender_clean))
cat("Age group:\n");      print(table(age_clean))
cat("Level of study:\n"); print(table(lev_simple))
cat("AIS frequency:\n");  print(table(freq_std))

# ── 5. Item-Level Descriptive Statistics ────────────────────
all_items <- cbind(IQ, SQ, SERVQ, MUX, CA, US)
desc      <- psych::describe(all_items)
print(round(desc[, c("n", "mean", "sd", "median", "min", "max", "skew", "kurtosis")], 3))

CS <- data.frame(
  IQ    = rowMeans(IQ,    na.rm = TRUE),
  SQ    = rowMeans(SQ,    na.rm = TRUE),
  SERVQ = rowMeans(SERVQ, na.rm = TRUE),
  MUX   = rowMeans(MUX,   na.rm = TRUE),
  CA    = rowMeans(CA,    na.rm = TRUE),
  US    = rowMeans(US,    na.rm = TRUE)
)
cat("\nConstruct-level descriptives:\n")
print(round(psych::describe(CS)[, c("n", "mean", "sd", "min", "max")], 3))

# ── 6. Reliability (Cronbach's α, ρ_A, CR, AVE) ───────────
constructs_lst <- list(IQ = IQ, SQ = SQ, SERVQ = SERVQ,
                       MUX = MUX, CA = CA, US = US)

alpha_vals <- sapply(constructs_lst, function(d)
  round(suppressWarnings(
    psych::alpha(d, check.keys = FALSE, warnings = FALSE)$total$raw_alpha), 3))
cat("\nCronbach's Alpha:\n"); print(alpha_vals)

cr_ave_df <- data.frame(
  Construct = character(), Alpha = numeric(),
  rho_A = numeric(), CR = numeric(),
  AVE = numeric(), sqrt_AVE = numeric(),
  stringsAsFactors = FALSE
)

for (nm in names(constructs_lst)) {
  fa_res <- suppressWarnings(tryCatch(
    psych::fa(constructs_lst[[nm]], nfactors = 1, rotate = "none", fm = "minres"),
    error = function(e)
      psych::fa(constructs_lst[[nm]], nfactors = 1, rotate = "none", fm = "pa")
  ))
  lam <- as.numeric(fa_res$loadings)
  lam[is.na(lam)] <- 0
  CR  <- round(sum(lam)^2 / (sum(lam)^2 + sum(1 - lam^2)), 3)
  AVE <- round(mean(lam^2), 3)

  # rho_A approximation (Dijkstra-Henseler)
  om <- tryCatch(
    suppressWarnings(psych::omega(constructs_lst[[nm]], nfactors = 1,
                                   plot = FALSE, warnings = FALSE)),
    error = function(e) NULL)
  rho_a <- if (!is.null(om)) round(om$omega.tot, 3) else NA_real_

  cr_ave_df <- rbind(cr_ave_df, data.frame(
    Construct = nm, Alpha = alpha_vals[nm], rho_A = rho_a,
    CR = CR, AVE = AVE, sqrt_AVE = round(sqrt(AVE), 3),
    stringsAsFactors = FALSE))
}
cat("\nReliability & Convergent Validity:\n"); print(cr_ave_df)

# ── 7. Discriminant Validity ───────────────────────────────
# Fornell-Larcker
cor_mat  <- cor(CS, use = "pairwise.complete.obs")
fl       <- cor_mat
diag(fl) <- cr_ave_df$sqrt_AVE
cat("\nFornell-Larcker Matrix (diagonal = sqrt(AVE)):\n"); print(round(fl, 3))

# HTMT
get_htmt <- function(a_df, b_df) {
  a <- scale(a_df); b <- scale(b_df)
  cross_cors <- cor(a, b, use = "pairwise.complete.obs")
  within_a   <- cor(a, use = "pairwise.complete.obs")
  within_b   <- cor(b, use = "pairwise.complete.obs")
  mc  <- mean(abs(cross_cors), na.rm = TRUE)
  wa  <- if (ncol(a) > 1) mean(abs(within_a[lower.tri(within_a)])) else 1
  wb  <- if (ncol(b) > 1) mean(abs(within_b[lower.tri(within_b)])) else 1
  round(mc / sqrt(wa * wb), 3)
}

clist  <- list(IQ = IQ, SQ = SQ, SERVQ = SERVQ, MUX = MUX, CA = CA, US = US)
cnames <- names(clist)
htmt_mat <- matrix(NA_real_, 6, 6, dimnames = list(cnames, cnames))
for (i in 1:6) for (j in 1:6)
  if (i != j) htmt_mat[i, j] <- get_htmt(clist[[i]], clist[[j]])
cat("\nHTMT Matrix:\n"); print(round(htmt_mat, 3))

# ── 8. PLS-SEM Measurement Model (seminr) ─────────────────
data_items <- cbind(IQ, SQ, SERVQ, MUX, CA, US)

mm <- constructs(
  composite("IQ",    multi_items("IQ",    1:4), weights = mode_A),
  composite("SQ",    multi_items("SQ",    1:4), weights = mode_A),
  composite("SERVQ", multi_items("SERVQ", 1:3), weights = mode_A),
  composite("MUX",   multi_items("MUX",   1:4), weights = mode_A),
  composite("CA",    multi_items("CA",    1:4), weights = mode_A),
  composite("US",    multi_items("US",    1:4), weights = mode_A)
)

sm_pls <- relationships(
  paths(from = c("IQ", "SQ", "SERVQ", "MUX", "CA"), to = "US")
)

fit         <- estimate_pls(data = data_items,
                            measurement_model = mm,
                            structural_model  = sm_pls)
summary_fit <- summary(fit)

cat("\n=== Outer Loadings ===\n")
print(round(summary_fit$loadings, 3))

cat("\n=== R² for US ===\n")
print(summary_fit$paths)

# ── 9. Bootstrap Inference ─────────────────────────────────
set.seed(2024)
boot         <- bootstrap_model(fit, nboot = 5000, cores = 1, seed = 2024)
summary_boot <- summary(boot)
bp_raw       <- as.data.frame(summary_boot$bootstrapped_paths)
cat("\n=== Bootstrapped Paths (raw) ===\n"); print(bp_raw)

# Robust column matching
col_beta  <- grep("Original",          colnames(bp_raw), value = TRUE)[1]
col_bmean <- grep("Bootstrap.*Mean|Mean", colnames(bp_raw), value = TRUE)[1]
col_sd    <- grep("Bootstrap.*SD|SE|SD",  colnames(bp_raw), value = TRUE)[1]
col_t     <- grep("T Stat",            colnames(bp_raw), value = TRUE)[1]
col_ci_lo <- grep("2\\.5",             colnames(bp_raw), value = TRUE)[1]
col_ci_hi <- grep("97\\.5",            colnames(bp_raw), value = TRUE)[1]
col_pval  <- grep("[Pp].*[Vv]al|p.val", colnames(bp_raw), value = TRUE)[1]

# If p-value column absent, compute from t-statistic (two-tailed)
if (is.na(col_pval)) {
  df_resid  <- nrow(data_items) - 6
  pvals_vec <- round(2 * pt(abs(bp_raw[[col_t]]), df = df_resid,
                             lower.tail = FALSE), 4)
} else {
  pvals_vec <- round(as.numeric(bp_raw[[col_pval]]), 4)
}

hyp <- data.frame(
  Hypothesis = c("H1: IQ -> US", "H2: SQ -> US", "H3: SERVQ -> US",
                 "H4: MUX -> US", "H5: CA -> US"),
  Beta      = round(as.numeric(bp_raw[[col_beta]]),  3),
  Boot_Mean = round(as.numeric(bp_raw[[col_bmean]]), 3),
  SE        = round(as.numeric(bp_raw[[col_sd]]),    3),
  T_Stat    = round(as.numeric(bp_raw[[col_t]]),     3),
  CI_Low    = round(as.numeric(bp_raw[[col_ci_lo]]), 3),
  CI_High   = round(as.numeric(bp_raw[[col_ci_hi]]), 3),
  p_value   = pvals_vec,
  Decision  = ifelse(pvals_vec < 0.05, "Supported", "Not Supported"),
  stringsAsFactors = FALSE
)
rownames(hyp) <- NULL
cat("\n=== Hypothesis Testing Results ===\n"); print(hyp)

# ── 10. Effect Sizes (f²) & VIF ───────────────────────────
CS_std <- as.data.frame(scale(CS))
m_ols  <- lm(US ~ IQ + SQ + SERVQ + MUX + CA, data = CS_std)
r2full <- summary(m_ols)$r.squared
preds  <- c("IQ", "SQ", "SERVQ", "MUX", "CA")

f2_vif <- data.frame(
  Predictor = character(), f2 = numeric(),
  Magnitude = character(), VIF = numeric(),
  stringsAsFactors = FALSE
)

pls_vif <- summary_fit$vif_antecedents$US
for (p in preds) {
  others <- setdiff(preds, p)
  r2red  <- summary(lm(as.formula(paste("US ~", paste(others, collapse = "+"))),
                        data = CS_std))$r.squared
  f2v    <- round((r2full - r2red) / (1 - r2full), 3)
  vif_v  <- round(as.numeric(pls_vif[p]), 3)
  mag    <- ifelse(f2v >= 0.35, "Large",
                   ifelse(f2v >= 0.15, "Medium",
                          ifelse(f2v >= 0.02, "Small", "Negligible")))
  f2_vif <- rbind(f2_vif, data.frame(
    Predictor = p, f2 = f2v, Magnitude = mag, VIF = vif_v,
    stringsAsFactors = FALSE))
}
cat("\nEffect Sizes and VIF:\n"); print(f2_vif)

r2_us  <- summary_fit$paths["R^2",    "US"]
r2_adj <- summary_fit$paths["AdjR^2", "US"]
cat(sprintf("\nR² = %.3f | Adjusted R² = %.3f\n", r2_us, r2_adj))


# ══════════════════════════════════════════════════════════════
# SECTION 11: FIGURES
# ══════════════════════════════════════════════════════════════
theme_thesis <- theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle   = element_text(hjust = 0.5, size = 10, color = "grey40"),
    axis.text       = element_text(size = 11),
    axis.title      = element_text(size = 11),
    legend.position = "none",
    plot.margin     = margin(10, 15, 10, 15)
  )

# ── Figure 2: Demographic Bar Charts (2×2 panel) ──────────
make_demo_bar <- function(var_vec, title, level_order) {
  df  <- data.frame(Category = factor(var_vec, levels = level_order))
  cnt <- df |>
    dplyr::count(Category) |>
    dplyr::mutate(pct   = round(n / sum(n) * 100, 1),
                  label = paste0(n, "\n(", pct, "%)"))
  ggplot(cnt, aes(x = Category, y = n, fill = Category)) +
    geom_col(width = 0.65, show.legend = FALSE) +
    geom_text(aes(label = label), vjust = -0.4, size = 3.3) +
    scale_fill_manual(values = rep(c("#4472C4", "#70AD47", "#ED7D31",
                                      "#FFC000", "#A5A5A5"),
                                    length.out = nrow(cnt))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = title, x = NULL, y = "Count") +
    theme_thesis
}

p_gender <- make_demo_bar(gender_clean, "Gender",
                          c("Female", "Male"))
p_age    <- make_demo_bar(age_clean, "Age Group",
                          c("Under 20", "20-24", "25-29", "30+"))
p_level  <- make_demo_bar(lev_simple, "Level of Study",
                          c("Undergraduate", "Master's", "PhD", "Exchange"))
p_freq   <- make_demo_bar(freq_std, "AIS Use Frequency",
                          c("Daily", "Weekly", "Monthly", "Rarely"))

fig2 <- gridExtra::arrangeGrob(
  p_gender, p_age, p_level, p_freq,
  ncol = 2,
  top  = grid::textGrob(
    "Figure 2. Demographic Characteristics of Respondents (N = 77)",
    gp = grid::gpar(fontface = "bold", fontsize = 13)
  )
)
grid::grid.draw(fig2)
ggsave("figures/Figure2_Demographic_Profile.png", fig2,
       width = 10, height = 8, dpi = 300, bg = "white")

# ── Figure 3: Construct Means with 95% CI ─────────────────
construct_labels <- c(
  IQ    = "Information\nQuality (IQ)",
  SQ    = "System\nQuality (SQ)",
  SERVQ = "Service\nQuality (SERVQ)",
  MUX   = "Mobile User\nExperience (MUX)",
  CA    = "Complementary\nAssets (CA)",
  US    = "User\nSatisfaction (US)"
)

ci_df <- CS |>
  tidyr::pivot_longer(everything(), names_to = "Construct", values_to = "Score") |>
  dplyr::group_by(Construct) |>
  dplyr::summarise(
    Mean    = mean(Score, na.rm = TRUE),
    SE      = sd(Score, na.rm = TRUE) / sqrt(n()),
    CI95_lo = Mean - 1.96 * SE,
    CI95_hi = Mean + 1.96 * SE,
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Construct = factor(Construct,
                       levels = c("IQ", "SQ", "SERVQ", "MUX", "CA", "US"),
                       labels = construct_labels)
  )

fig3 <- ggplot(ci_df, aes(x = Construct, y = Mean, fill = Construct)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = CI95_lo, ymax = CI95_hi),
                width = 0.25, linewidth = 0.7, color = "grey30") +
  geom_text(aes(y = CI95_hi, label = sprintf("%.2f", Mean)),
            vjust = -0.5, size = 3.8, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, 5.6),
                     breaks = seq(0, 5, 1),
                     expand = expansion(mult = c(0, 0))) +
  geom_hline(yintercept = 3, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  labs(
    title    = "Figure 3. Mean Construct Scores with 95% Confidence Intervals",
    subtitle = "N = 77 | Scale: 1 (Strongly Disagree) – 5 (Strongly Agree)",
    x = "Construct", y = "Mean Score"
  ) +
  theme_thesis +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 9, lineheight = 1.2),
    plot.margin     = margin(t = 10, r = 15, b = 20, l = 15)
  )
print(fig3)
ggsave("figures/Figure3_Construct_Means_CI.png", fig3,
       width = 12, height = 7, dpi = 300, bg = "white")

# ── Figure 4: Outer Loadings Bar Chart ────────────────────
loadings_df <- data.frame(
  Item      = c("IQ1", "IQ2", "IQ3", "IQ4",
                "SQ1", "SQ2", "SQ3", "SQ4",
                "SERVQ1", "SERVQ2", "SERVQ3",
                "MUX1", "MUX2", "MUX3", "MUX4",
                "CA1", "CA2", "CA3", "CA4",
                "US1", "US2", "US3", "US4"),
  Construct = c(rep("IQ", 4), rep("SQ", 4), rep("SERVQ", 3),
                rep("MUX", 4), rep("CA", 4), rep("US", 4)),
  Loading   = c(0.618, 0.726, 0.792, 0.673,
                0.852, 0.559, 0.707, 0.852,
                0.832, 0.843, 0.833,
                0.816, 0.893, 0.768, 0.866,
                0.609, 0.369, 0.866, 0.839,
                0.901, 0.897, 0.830, 0.880),
  stringsAsFactors = FALSE
)
loadings_df$Item <- factor(loadings_df$Item,
                           levels = rev(unique(loadings_df$Item)))
loadings_df$Construct <- factor(loadings_df$Construct,
                                levels = c("IQ", "SQ", "SERVQ", "MUX", "CA", "US"))

fig4 <- ggplot(loadings_df, aes(y = Item, x = Loading, fill = Construct)) +
  geom_col(width = 0.65) +
  geom_vline(xintercept = 0.70, color = "red",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", y = 1.5, x = 0.72,
           label = "Threshold = 0.70", color = "red",
           size = 3.2, hjust = 0) +
  geom_text(aes(label = sprintf("%.3f", Loading)),
            hjust = -0.1, size = 3.1) +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1.0, 0.2)) +
  labs(
    title    = "Figure 4. Outer Loadings for All Measurement Items (PLS-SEM)",
    subtitle = "N = 77 | Red dashed line = threshold \u03bb \u2265 0.70 (Hair et al., 2019)",
    y = NULL, x = "Standardized Loading (\u03bb)", fill = "Construct"
  ) +
  theme_thesis +
  theme(legend.position = "right")
print(fig4)
ggsave("figures/Figure4_Outer_Loadings.png", fig4,
       width = 10, height = 8, dpi = 300, bg = "white")

# ── Figure 5: Reliability & Validity Grouped Bar Chart ────
rel_df <- data.frame(
  Construct = rep(c("IQ", "SQ", "SERVQ", "MUX", "CA", "US"), each = 4),
  Metric    = rep(c("Cronbach's \u03b1", "\u03c1A", "CR", "AVE"), 6),
  Value     = c(
    0.673, 0.687, 0.692, 0.376,
    0.745, 0.749, 0.753, 0.443,
    0.774, 0.790, 0.790, 0.559,
    0.856, 0.862, 0.863, 0.611,
    0.694, 0.700, 0.701, 0.376,
    0.900, 0.901, 0.901, 0.695
  ),
  stringsAsFactors = FALSE
)
rel_df$Construct <- factor(rel_df$Construct,
                           levels = c("IQ", "SQ", "SERVQ", "MUX", "CA", "US"))
rel_df$Metric    <- factor(rel_df$Metric,
                           levels = c("Cronbach's \u03b1", "\u03c1A", "CR", "AVE"))

fig5 <- ggplot(rel_df, aes(x = Construct, y = Value, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, alpha = 0.9) +
  geom_hline(yintercept = 0.70, color = "#C00000", linetype = "dashed",
             linewidth = 0.9) +
  geom_hline(yintercept = 0.50, color = "#375623", linetype = "dashed",
             linewidth = 0.9) +
  annotate("text", x = 6.55, y = 0.72,
           label = "\u03b1/\u03c1A/CR threshold = 0.70",
           color = "#C00000", size = 3.6, hjust = 0, fontface = "bold") +
  annotate("text", x = 6.55, y = 0.52,
           label = "AVE threshold = 0.50",
           color = "#375623", size = 3.6, hjust = 0, fontface = "bold") +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1.0, 0.2)) +
  coord_cartesian(clip = "off") +
  labs(
    title    = "Figure 5. Construct Reliability and Convergent Validity",
    subtitle = "N = 77 | Dashed lines = thresholds (Hair et al., 2019)",
    x = NULL, y = "Value", fill = "Metric"
  ) +
  theme_thesis +
  theme(
    legend.position = "bottom",
    plot.margin     = margin(t = 10, r = 160, b = 10, l = 10)
  )
print(fig5)
ggsave("figures/Figure5_Reliability_Validity.png", fig5,
       width = 11, height = 6, dpi = 300, bg = "white")

# ── Figure 6: Fornell-Larcker Heatmap ─────────────────────
constructs <- c("IQ", "SQ", "SERVQ", "MUX", "CA", "US")

fl_correct <- matrix(NA, nrow = 6, ncol = 6,
                     dimnames = list(constructs, constructs))

fl_correct["IQ",    "IQ"]    <- 0.613
fl_correct["SQ",    "IQ"]    <- 0.609; fl_correct["SQ",    "SQ"]    <- 0.666
fl_correct["SERVQ", "IQ"]    <- 0.515; fl_correct["SERVQ", "SQ"]   <- 0.801; fl_correct["SERVQ", "SERVQ"] <- 0.748
fl_correct["MUX",   "IQ"]    <- 0.383; fl_correct["MUX",   "SQ"]   <- 0.703; fl_correct["MUX",   "SERVQ"] <- 0.630; fl_correct["MUX", "MUX"] <- 0.782
fl_correct["CA",    "IQ"]    <- 0.436; fl_correct["CA",    "SQ"]   <- 0.452; fl_correct["CA",    "SERVQ"] <- 0.482; fl_correct["CA",  "MUX"] <- 0.407; fl_correct["CA", "CA"] <- 0.613
fl_correct["US",    "IQ"]    <- 0.653; fl_correct["US",    "SQ"]   <- 0.798; fl_correct["US",    "SERVQ"] <- 0.818; fl_correct["US",  "MUX"] <- 0.627; fl_correct["US", "CA"] <- 0.456; fl_correct["US", "US"] <- 0.834

fl_long <- as.data.frame(fl_correct) |>
  tibble::rownames_to_column("Row") |>
  tidyr::pivot_longer(-Row, names_to = "Col", values_to = "Value") |>
  dplyr::mutate(
    Row  = factor(Row, levels = rev(constructs)),
    Col  = factor(Col, levels = constructs),
    Diag = (as.character(Row) == as.character(Col))
  )

fig6 <- ggplot(fl_long, aes(x = Col, y = Row, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(is.na(Value), "", sprintf("%.3f", Value)),
                fontface = ifelse(Diag, "bold", "plain")),
            size = 3.8) +
  scale_fill_gradient2(
    low = "#FFFFFF", mid = "#9DC3E6", high = "#2E75B6",
    midpoint = 0.5, na.value = "grey90", name = "Value"
  ) +
  labs(
    title    = "Figure 6. Fornell\u2013Larcker Discriminant Validity Matrix",
    subtitle = "Diagonal = \u221aAVE; Off-diagonal = inter-construct correlations",
    x = NULL, y = NULL
  ) +
  theme_thesis +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right")
print(fig6)
ggsave("figures/Figure6_Fornell_Larcker_Heatmap.png", fig6,
       width = 8, height = 6, dpi = 300, bg = "white")

# ── Figure 7: Structural Path Model ──────────────────────
make_hex <- function(cx, cy, rx = 0.88, ry = 0.68) {
  angles <- (0:5) * 60 * pi / 180
  data.frame(x = cx + rx * cos(angles), y = cy + ry * sin(angles))
}

cx_exog <- 5.5;  cx_us <- 9.8
rx <- 0.88;      ry <- 0.68
ix_exog <- 1.3;  ix_us <- 14.0
ibw <- 0.52;     ibh <- 0.23

cy_vals <- c(SQ = 9.5, SERVQ = 7.9, MUX = 6.3, CA = 4.7, IQ = 3.1)
us_y <- 6.3

items_list <- list(
  SQ    = data.frame(item = c("SQ1", "SQ2", "SQ3", "SQ4"),
                     lam = c(0.852, 0.559, 0.707, 0.852)),
  SERVQ = data.frame(item = c("SERVQ1", "SERVQ2", "SERVQ3"),
                     lam = c(0.832, 0.843, 0.833)),
  MUX   = data.frame(item = c("MUX1", "MUX2", "MUX3", "MUX4"),
                     lam = c(0.816, 0.893, 0.768, 0.866)),
  CA    = data.frame(item = c("CA1", "CA2", "CA3", "CA4"),
                     lam = c(0.609, 0.369, 0.866, 0.839)),
  IQ    = data.frame(item = c("IQ1", "IQ2", "IQ3", "IQ4"),
                     lam = c(0.618, 0.726, 0.792, 0.673)),
  US    = data.frame(item = c("US1", "US2", "US3", "US4"),
                     lam = c(0.901, 0.897, 0.830, 0.880))
)

item_df <- do.call(rbind, lapply(names(items_list), function(cn) {
  df  <- items_list[[cn]]
  n   <- nrow(df)
  by  <- if (cn == "US") us_y else cy_vals[cn]
  sp  <- 0.43
  offs <- if (n == 1) 0 else seq((n - 1) * sp / 2, -(n - 1) * sp / 2, length.out = n)
  data.frame(item = df$item, lam = df$lam,
             x = if (cn == "US") ix_us else ix_exog,
             y = by + offs, construct = cn, stringsAsFactors = FALSE)
}))

all_cn <- c(names(cy_vals), "US")
cx_all <- c(rep(cx_exog, 5), cx_us)
cy_all <- c(as.numeric(cy_vals), us_y)

hex_poly <- do.call(rbind, lapply(seq_along(all_cn), function(i) {
  h <- make_hex(cx_all[i], cy_all[i], rx, ry)
  h$construct <- all_cn[i]; h$gid <- i; h
}))
constr_pos <- data.frame(name = all_cn, cx = cx_all, cy = cy_all,
                         stringsAsFactors = FALSE)

arrow_exog <- do.call(rbind, lapply(names(cy_vals), function(cn) {
  its <- item_df[item_df$construct == cn, ]
  data.frame(xs = cx_exog - rx, ys = cy_vals[cn],
             xe = ix_exog + ibw, ye = its$y,
             lam_label = paste0("\u03bb = ", sprintf("%.3f", its$lam)),
             stringsAsFactors = FALSE)
}))

arrow_us <- {
  its <- item_df[item_df$construct == "US", ]
  data.frame(xs = cx_us + rx, ys = us_y, xe = ix_us - ibw, ye = its$y,
             lam_label = paste0("\u03bb = ", sprintf("%.3f", its$lam)),
             stringsAsFactors = FALSE)
}

p1 <- 0.65; p2 <- 0.55
arrow_exog$lx <- arrow_exog$xs + p1 * (arrow_exog$xe - arrow_exog$xs)
arrow_exog$ly <- arrow_exog$ys + p1 * (arrow_exog$ye - arrow_exog$ys) + 0.09
arrow_us$lx   <- arrow_us$xs + p2 * (arrow_us$xe - arrow_us$xs)
arrow_us$ly   <- arrow_us$ys + p2 * (arrow_us$ye - arrow_us$ys) + 0.09

struct_df <- data.frame(
  beta = c(0.211, 0.413, 0.103, 0.057, 0.250),
  sig  = c(FALSE, TRUE, FALSE, FALSE, TRUE),
  xs = cx_exog + rx, ys = as.numeric(cy_vals),
  xe = cx_us - rx,   ye = us_y, stringsAsFactors = FALSE
)
struct_df$blabel <- paste0("\u03b2 = ", sprintf("%.3f", struct_df$beta))
struct_df$lx <- 7.55
struct_df$ly <- (struct_df$ys + struct_df$ye) / 2 + 0.16

fig7 <- ggplot() +
  geom_segment(data = arrow_exog, aes(x = xs, y = ys, xend = xe, yend = ye),
               arrow = arrow(length = unit(0.12, "cm"), type = "open"),
               color = "grey55", linewidth = 0.4) +
  geom_segment(data = arrow_us, aes(x = xs, y = ys, xend = xe, yend = ye),
               arrow = arrow(length = unit(0.12, "cm"), type = "open"),
               color = "grey55", linewidth = 0.4) +
  geom_text(data = arrow_exog, aes(x = lx, y = ly, label = lam_label),
            size = 2.3, hjust = 1.05) +
  geom_text(data = arrow_us, aes(x = lx, y = ly, label = lam_label),
            size = 2.3, hjust = 0.0) +
  geom_segment(data = struct_df[!struct_df$sig, ],
               aes(x = xs, y = ys, xend = xe, yend = ye),
               arrow = arrow(length = unit(0.13, "cm"), type = "closed"),
               color = "black", linewidth = 0.55) +
  geom_segment(data = struct_df[struct_df$sig, ],
               aes(x = xs, y = ys, xend = xe, yend = ye),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               color = "black", linewidth = 1.6) +
  geom_text(data = struct_df, aes(x = lx, y = ly, label = blabel),
            size = 2.9, hjust = 0) +
  geom_polygon(data = hex_poly, aes(x = x, y = y, group = gid),
               fill = "white", color = "black", linewidth = 0.65) +
  geom_text(data = constr_pos[constr_pos$name != "US", ],
            aes(x = cx, y = cy, label = name), size = 4.0, fontface = "bold") +
  annotate("text", x = cx_us, y = us_y + 0.20, label = "US",
           size = 4.0, fontface = "bold") +
  annotate("text", x = cx_us, y = us_y - 0.20,
           label = paste0("R\u00b2 = ", round(r2_us, 3)), size = 2.8) +
  geom_rect(data = item_df,
            aes(xmin = x - ibw, xmax = x + ibw,
                ymin = y - ibh, ymax = y + ibh),
            fill = "white", color = "black", linewidth = 0.45) +
  geom_text(data = item_df, aes(x = x, y = y, label = item), size = 2.6) +
  coord_cartesian(xlim = c(-0.3, 15.5), ylim = c(2.3, 10.3)) +
  labs(
    title    = "Figure 7. Structural Path Model with Path Coefficients (\u03b2) & Outer Loadings (\u03bb)",
    subtitle = "N = 77 | Thick arrows = significant (p < .05) | Thin arrows = non-significant",
    x = NULL, y = NULL
  ) +
  theme_thesis +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), legend.position = "none")
print(fig7)
ggsave("figures/Figure7_Path_Model.png", fig7,
       width = 14, height = 7.5, dpi = 300, bg = "white")

# ── Figure 8: Path Coefficients Bar Chart ─────────────────
path_df <- hyp |>
  dplyr::mutate(
    Sig      = ifelse(p_value < 0.05, "Significant", "Non-Significant"),
    HypLabel = sub(".*: ", "", Hypothesis)
  )

fig8 <- ggplot(path_df, aes(y = reorder(HypLabel, Beta),
                            x = Beta, fill = Sig)) +
  geom_col(width = 0.60, alpha = 0.9) +
  geom_errorbar(aes(xmin = CI_Low, xmax = CI_High),
                height = 0.22, linewidth = 0.7, color = "grey30") +
  geom_text(aes(label = paste0("\u03b2 = ", sprintf("%.3f", Beta),
                               "\np = ", sprintf("%.4f", p_value))),
            hjust = -0.1, size = 3.3) +
  scale_fill_manual(values = c("Significant"     = "#4472C4",
                               "Non-Significant" = "#A5A5A5")) +
  scale_x_continuous(limits = c(-0.15, 0.70),
                     breaks = seq(-0.1, 0.6, 0.1)) +
  labs(
    title    = "Figure 8. Standardized Path Coefficients to User Satisfaction",
    subtitle = paste0("R\u00b2 = ", round(r2_us, 3),
                      " | Adjusted R\u00b2 = ", round(r2_adj, 3),
                      " | Bootstrap N = 5,000 | Seed = 2024"),
    y = NULL, x = "Standardized Path Coefficient (\u03b2)", fill = NULL
  ) +
  theme_thesis +
  theme(legend.position = "bottom")
print(fig8)
ggsave("figures/Figure8_Path_Coefficients.png", fig8,
       width = 10, height = 6, dpi = 300, bg = "white")

# ── 12. Save Analysis Objects ──────────────────────────────
save(raw, IQ, SQ, SERVQ, MUX, CA, US,
     all_items, CS, desc, cr_ave_df,
     alpha_vals, cor_mat, fl_correct, htmt_mat,
     m_ols, f2_vif, hyp,
     fit, boot, summary_fit, summary_boot,
     rel_df, loadings_df,
     file = "data/NTHU_AIS_analysis_final.RData")
cat("\nAll results saved to data/NTHU_AIS_analysis_final.RData\n")
cat("Script complete.\n")
