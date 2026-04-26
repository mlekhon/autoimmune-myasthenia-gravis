# Delayed diagnosis via CASANOVA (GFDsurv)

# install.packages("GFDsurv")   # uncomment if not yet installed
# install.packages("survival")  # for Surv objects if you want plots later

library(GFDsurv)

# ---- Load data ----
data_dir <- "/Users/mdmuhtashimlekhon/Desktop/R project"
df <- read.csv(file.path(data_dir, "df_original.csv"), check.names = FALSE)

# ---- Keep required columns & enforce types ----
dat <- subset(df, select = c(time_erst_to_diag, status_diagnosis, seronegative, sex))

dat$time_erst_to_diag <- as.numeric(dat$time_erst_to_diag)
dat$status_diagnosis  <- as.integer(dat$status_diagnosis)           # 1 = event occurred (diagnosis), 0 = censored
dat$seronegative      <- factor(dat$seronegative, levels = c(0,1), labels = c("Seropositive","Seronegative"))
dat$sex               <- factor(dat$sex,          levels = c(0,1), labels = c("Female","Male"))

# Optional safety: drop non-finite/zero times; add tiny epsilon to zeros if any slipped through
dat <- dat[is.finite(dat$time_erst_to_diag), ]
dat$time_erst_to_diag[dat$time_erst_to_diag <= 0] <- 1e-6

# ---- Run CASANOVA ----
set.seed(123)  # for permutation reproducibility
fit <- casanova(
  formula = "time_erst_to_diag ~ seronegative*sex",
  event   = "status_diagnosis",
  data    = dat,
  nperm   = 4999,                 # increase permutations for stable perm p-values
  cross   = TRUE,                 # include crossing weight (handles NPH / crossing hazards)
  rg      = list(c(0,0), c(1,0), c(0,1))  # log-rank + early/late emphasis weights
)

# ---- Results ----
print(fit)
summary(fit)

# Cell sizes (sanity check)
with(dat, table(seronegative, sex))

# Packages
# install.packages(c("survival","survminer","dplyr"))
library(survival)
library(survminer)
library(dplyr)

# Load your clean dataset
data_dir <- "/Users/mdmuhtashimlekhon/Desktop/R project"
df <- read.csv(file.path(data_dir, "df_original.csv"), check.names = FALSE)

# Keep/format needed columns
dat <- df[, c("time_erst_to_diag","status_diagnosis","seronegative","sex")]
dat$time_erst_to_diag <- as.numeric(dat$time_erst_to_diag)
dat$status_diagnosis  <- as.integer(dat$status_diagnosis)    # should be 1 for all kept rows
dat$seronegative      <- factor(dat$seronegative, levels=c(0,1), labels=c("Seropositive","Seronegative"))
dat$sex               <- factor(dat$sex,          levels=c(0,1), labels=c("Female","Male"))

# --- 1) KM curves by Serostatus ---
fit_sero <- survfit(Surv(time_erst_to_diag, status_diagnosis) ~ seronegative, data = dat)
p_sero <- ggsurvplot(
  fit_sero, data = dat, conf.int = TRUE, surv.median.line = "hv",
  xlab = "Years from symptom onset to diagnosis",
  ylab = "Not yet diagnosed (S(t))",
  legend.title = "Serostatus", legend.labs = levels(dat$seronegative),
  ggtheme = theme_minimal()
)

# --- 2) KM curves by Sex ---
fit_sex <- survfit(Surv(time_erst_to_diag, status_diagnosis) ~ sex, data = dat)
p_sex <- ggsurvplot(
  fit_sex, data = dat, conf.int = TRUE, surv.median.line = "hv",
  xlab = "Years from symptom onset to diagnosis",
  ylab = "Not yet diagnosed (S(t))",
  legend.title = "Sex", legend.labs = levels(dat$sex),
  ggtheme = theme_minimal()
)

# --- 3) Optional: 4-curve KM (Sex × Serostatus) ---
dat$group <- interaction(dat$sex, dat$seronegative, sep=" × ", drop=TRUE)
fit_4 <- survfit(Surv(time_erst_to_diag, status_diagnosis) ~ group, data = dat)
p_4 <- ggsurvplot(
  fit_4, data = dat, conf.int = FALSE,
  xlab = "Years from symptom onset to diagnosis",
  ylab = "Not yet diagnosed (S(t))",
  legend.title = "Group", legend.labs = levels(dat$group),
  ggtheme = theme_minimal()
)

# --- 4) Nelson–Aalen cumulative hazard (aligns with CASANOVA) ---
p_ch <- ggsurvplot(
  fit_sero, data = dat, fun = "cumhaz",
  xlab = "Years from symptom onset to diagnosis",
  ylab = "Cumulative hazard H(t)",
  legend.title = "Serostatus", legend.labs = levels(dat$seronegative),
  ggtheme = theme_minimal()
)

# Print plots
print(p_sero); print(p_sex); print(p_4); print(p_ch)

# --- 5) Descriptive table (report beside CASANOVA p-values) ---
tab <- dat %>%
  group_by(seronegative, sex) %>%
  summarise(
    n      = n(),
    mean   = round(mean(time_erst_to_diag), 2),
    median = round(median(time_erst_to_diag), 2),
    q1     = round(quantile(time_erst_to_diag, 0.25), 2),
    q3     = round(quantile(time_erst_to_diag, 0.75), 2),
    .groups = "drop"
  )
print(tab)

# --- Nelson–Aalen cumulative hazard for 4 groups ----------------------------
# Requires: survival, survminer, dplyr

library(survival)
library(survminer)
library(dplyr)

# --- Load / prepare (adjust path if needed) ---
data_dir <- "/Users/mdmuhtashimlekhon/Desktop/R project"
df <- read.csv(file.path(data_dir, "df_original.csv"), check.names = FALSE)

dat <- df[, c("time_erst_to_diag","status_diagnosis","seronegative","sex")]
dat$time_erst_to_diag <- as.numeric(dat$time_erst_to_diag)
dat$status_diagnosis  <- as.integer(dat$status_diagnosis)  # 1 = diagnosed, 0 = censored
dat$seronegative      <- factor(dat$seronegative, levels = c(0,1),
                                labels = c("Seropositive","Seronegative"))
dat$sex               <- factor(dat$sex, levels = c(0,1),
                                labels = c("Female","Male"))

# Optional safety
dat <- dat[is.finite(dat$time_erst_to_diag), ]
dat$time_erst_to_diag[dat$time_erst_to_diag <= 0] <- 1e-6

# --- Build 4-group factor (Sex × Serostatus) with readable labels/order ---
dat$group <- interaction(dat$sex, dat$seronegative, sep = " × ", drop = TRUE)

# Relevel to a clear order (change if you prefer a different order)
dat$group <- factor(dat$group,
                    levels = c("Male × Seropositive",
                               "Female × Seropositive",
                               "Male × Seronegative",
                               "Female × Seronegative")
)

# --- Fit survival object per 4 groups ---
fit_4 <- survfit(Surv(time_erst_to_diag, status_diagnosis) ~ group, data = dat)

# --- Plot Nelson–Aalen cumulative hazard ---
p_ch_4 <- ggsurvplot(
  fit_4, data = dat, fun = "cumhaz",
  xlab = "Years from symptom onset to diagnosis",
  ylab = "Cumulative hazard H(t)",
  legend.title = "Group",
  legend.labs = levels(dat$group),
  ggtheme = theme_minimal(),
  conf.int = FALSE
)

print(p_ch_4)

# --- (Optional) Save the plot ---
ggsave("NA_cumhaz_4groups.png", p_ch_4$plot, width = 7, height = 5, dpi = 300)

# --- Export underlying Nelson–Aalen curves to CSV for each group -----------
# survminer::surv_summary gives a tidy data.frame of the survfit object
haz_df <- survminer::surv_summary(fit_4, data = dat)
# For fun="cumhaz", surv_summary stores cumulative hazard in 'cumhaz'
# (If your version stores it differently, it may be in 'n.event'/'n.risk' and 'cumhaz' computed internally.)
# Keep relevant columns and rename nicely if present
keep_cols <- intersect(c(".strata","time","cumhaz","n.risk","n.event"), names(haz_df))
haz_out <- haz_df[, keep_cols]
# Map strata to our group labels (remove "group=" prefix if present)
haz_out$Group <- gsub("^group=", "", haz_df$.strata)

write.csv(haz_out, "NA_cumhaz_4groups.csv", row.names = FALSE)
cat("Saved: NA_cumhaz_4groups.csv\n")

# --- Hazard summary table by serostatus and sex -----------------------
