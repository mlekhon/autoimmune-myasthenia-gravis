# CASANOVA-style analysis in R (MANCOVA for multiple outcomes with covariates)
# Research question:
# Does Seronegativity in Myasthenia Gravis patients correlate with delayed diagnosis and worse functional outcomes?

# ============================
# How to use this script
# ============================
# 1) Put this script and the CSV "casanova_data.csv" in the same folder, or change `csv_path` below.
# 2) Open this file in RStudio and click "Run" (or run line by line).
# 3) Outputs will print to the Console, and a few CSVs will be saved alongside your data.

# ----------------------------
# Package setup
# ----------------------------
need <- c("car", "heplots", "effectsize", "broom")
to_install <- need[!(need %in% installed.packages()[,"Package"])]
if(length(to_install)) install.packages(to_install, dependencies = TRUE)

library(car)         # Anova utilities
library(heplots)     # boxM test
library(effectsize)  # effect sizes
library(broom)       # tidy summaries

# ----------------------------
# Load data
# ----------------------------
# If your CSV lives elsewhere, replace with the absolute path.
csv_path <- "casanova_data.csv"
df <- read.csv(csv_path, stringsAsFactors = FALSE)

# ----------------------------
# Variable preparation
# ----------------------------
# Convert categorical variables to factors and label them
df$seronegative <- factor(df$seronegative, levels = c(0,1), labels = c("Seropositive","Seronegative"))
df$sex          <- factor(df$sex, levels = c(0,1), labels = c("Male","Female"))
df$thymom       <- factor(df$thymom, levels = c(0,1), labels = c("No","Yes"))
df$ee_autoimmunerkrankungen_rbzu <- factor(df$ee_autoimmunerkrankungen_rbzu)

# Outcomes (delayed diagnosis + functional outcomes)
Y <- cbind(
  time_erst_to_diag,
  mgfa_aktuell_1bis5,
  scoreqmg_neu,
  scoreadl_neu,
  scoreqol_neu
)

# ----------------------------
# Quick sanity checks
# ----------------------------
message("N (rows) in dataset: ", nrow(df))
message("Seronegativity counts:")
print(table(df$seronegative, useNA = "ifany"))
message("Missing values check (should be all FALSE because we cleaned earlier):")
print(sapply(df[,c('time_erst_to_diag','mgfa_aktuell_1bis5','scoreqmg_neu','scoreadl_neu','scoreqol_neu',
                   'seronegative','age_erst','sex','thymom','dauer','ee_autoimmunerkrankungen_rbzu')], function(x) any(is.na(x))))

# ----------------------------
# Core CASANOVA-style model (MANCOVA)
# ----------------------------
manova_model <- manova(Y ~ seronegative + age_erst + sex + thymom + dauer + ee_autoimmunerkrankungen_rbzu, data = df)

# Multivariate tests (report both Wilks and Pillai for robustness)
cat("\n========== Multivariate Tests (Wilks' Lambda) ==========\n")
print(summary(manova_model, test = "Wilks"))

cat("\n========== Multivariate Tests (Pillai's Trace) ==========\n")
print(summary(manova_model, test = "Pillai"))

# ----------------------------
# Box's M test for homogeneity of covariance matrices (assumption check)
# ----------------------------
cat("\n========== Box's M Test ==========\n")
# Grouping factor is the main group: seronegative
bm <- boxM(Y, df$seronegative)
print(bm)

# ----------------------------
# Univariate follow-up ANCOVAs
# ----------------------------
cat("\n========== Univariate ANCOVA (per outcome) ==========\n")
univ <- summary.aov(manova_model)   # list of ANOVA tables for each outcome
print(univ)

# Save tidy univariate tables to CSV for convenience
univ_tidy_list <- lapply(names(univ), function(nm) {
  tt <- univ[[nm]]
  # Convert to data.frame and add outcome name
  df_t <- as.data.frame(tt)
  df_t$Term <- rownames(df_t)
  rownames(df_t) <- NULL
  df_t$Outcome <- nm
  df_t[,c("Outcome","Term",colnames(df_t)[colnames(df_t)!="Outcome" & colnames(df_t)!="Term"])]
})

univ_tidy <- do.call(rbind, univ_tidy_list)
write.csv(univ_tidy, file = "casanova_univariate_tables.csv", row.names = FALSE)

# ----------------------------
# Partial eta^2 effect sizes for each outcome & term
# ----------------------------
# Compute partial eta^2 from the ANOVA tables: SS_effect / (SS_effect + SS_residual)
# Different outcomes can have different residual SS; extract per outcome.
get_partial_eta2 <- function(aov_table) {
  # aov_table is a summary.aov single-element table
  df_tab <- as.data.frame(aov_table)
  df_tab$Term <- rownames(df_tab)
  rownames(df_tab) <- NULL
  # locate residual row (usually "Residuals")
  res_row <- df_tab[df_tab$Term %in% c("Residuals", "Residual"), , drop = FALSE]
  if(nrow(res_row) != 1) return(NULL)
  ss_res <- res_row[["Sum Sq"]]
  # exclude residual row for effect size calc
  eff <- df_tab[!(df_tab$Term %in% c("Residuals","Residual")), , drop = FALSE]
  if(nrow(eff) == 0) return(NULL)
  eff$partial_eta2 <- eff[["Sum Sq"]] / (eff[["Sum Sq"]] + ss_res)
  eff
}

effect_sizes <- lapply(univ, get_partial_eta2)
names(effect_sizes) <- names(univ)
# Bind and save
es_tidy_list <- lapply(names(effect_sizes), function(nm) {
  es <- effect_sizes[[nm]]
  if(is.null(es)) return(NULL)
  es$Outcome <- nm
  es[, c("Outcome","Term","Df","Sum Sq","Mean Sq","F value","Pr(>F)","partial_eta2")]
})

es_tidy <- do.call(rbind, es_tidy_list)
if(!is.null(es_tidy)) {
  write.csv(es_tidy, file = "casanova_partial_eta2.csv", row.names = FALSE)
  cat("\n========== Partial eta^2 (saved to casanova_partial_eta2.csv) ==========\n")
  print(es_tidy[order(es_tidy$Outcome, -es_tidy$partial_eta2), ])
} else {
  cat("\n(No effect sizes computed; check ANOVA tables structure.)\n")
}

# ----------------------------
# Multiple-comparison correction for univariate p-values (BH/FDR)
# ----------------------------
if(!is.null(univ_tidy)) {
  # one row per term per outcome; filter out Residuals
  univ_terms <- subset(univ_tidy, !(Term %in% c("Residuals", "Residual")))
  # Rename p-value column if necessary
  pv_col <- if("Pr..F." %in% names(univ_terms)) "Pr..F." else if("Pr(>F)" %in% names(univ_terms)) "Pr(>F)" else NA
  if(!is.na(pv_col)) {
    univ_terms$padj_BH <- p.adjust(univ_terms[[pv_col]], method = "BH")
    write.csv(univ_terms, file = "casanova_univariate_pvals_BH.csv", row.names = FALSE)
    cat("\n========== Univariate p-values with BH correction (saved to casanova_univariate_pvals_BH.csv) ==========\n")
    print(univ_terms[, c("Outcome","Term", pv_col, "padj_BH")])
  }
}

cat("\n\nAnalysis complete. Key files written:\n",
    "- casanova_univariate_tables.csv\n",
    "- casanova_partial_eta2.csv\n",
    "- casanova_univariate_pvals_BH.csv\n", sep="")
