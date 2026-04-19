#!/usr/bin/env Rscript
# =============================================================================
# STEP 4b: GWAS — Survival Analysis (Cox Proportional Hazards)
# AOO GWAS Pipeline
# =============================================================================
# PURPOSE:
#   Model AOO as a time-to-event outcome using Cox PH regression.
#   This complements the linear model because:
#     - Incorporates censored controls (subjects who have not developed AD)
#     - Handles the right-skewed distribution of AOO
#     - Provides hazard ratios — clinically interpretable effect sizes
#
# MODEL:
#   For each SNP: Cox ~ SNP + sex + age_at_assessment + APOE_e4_count +
#                       PC1..PC10 + site
#   Time  = AOO (cases) or age at assessment (controls, censored)
#   Event = 1 (AD/MCI) or 0 (NCI)
#
# NOTE: Full genome-wide Cox is computationally expensive.
#       We use PLINK2 --cox for efficiency, then validate top hits with
#       the R survival package for full model diagnostics.
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(survival)
    library(survminer)
    library(ggplot2)
    library(gridExtra)
})

# ─── CONFIG ──────────────────────────────────────────────────────────────────
PIPELINE_DIR <- Sys.getenv("PIPELINE_DIR", ".")
GWAS_DIR     <- file.path(PIPELINE_DIR, "04_gwas", "cox")
FIG_DIR      <- file.path(PIPELINE_DIR, "09_figures")
dir.create(GWAS_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIG_DIR,  showWarnings=FALSE)

log_file <- file.path(PIPELINE_DIR, "logs", "step04b_cox.log")
con <- file(log_file, "a")
sink(con, append=TRUE, type="output")
sink(con, append=TRUE, type="message")

cat("============================================================\n")
cat(sprintf(" STEP 4b: COX GWAS  [%s]\n", Sys.time()))
cat("============================================================\n")

checkpoint <- function(name, condition, note="") {
    if (isTRUE(condition)) {
        cat(sprintf("  ✅ CHECKPOINT: %s\n", name))
    } else {
        cat(sprintf("  ❌ CHECKPOINT FAILED: %s %s\n", name, note))
        stop(paste("Pipeline halted at checkpoint:", name))
    }
}

# ─── LOAD PHENOTYPE & PCA ────────────────────────────────────────────────────
cat("\n--- Loading phenotype data ---\n")
pheno <- fread(file.path(PIPELINE_DIR, "01_phenotype", "pheno_cox.txt"))
pca   <- fread(file.path(PIPELINE_DIR, "03_pca", "pca.eigenvec"))
setnames(pca, "#FID", "FID")

df <- merge(pheno, pca[, c("FID","IID",paste0("PC",1:10)), with=FALSE],
            by=c("FID","IID"), all.x=TRUE)

cat(sprintf("  Loaded: %d samples\n", nrow(df)))
cat(sprintf("  Events (AD/MCI): %d  Controls (NCI): %d\n",
            sum(df$cox_event==1), sum(df$cox_event==0)))

checkpoint("Phenotype loaded", nrow(df) > 500)
checkpoint("PCs merged",       "PC1" %in% names(df))

# ─── KAPLAN-MEIER CURVES ─────────────────────────────────────────────────────
cat("\n--- Generating Kaplan-Meier survival curves ---\n")

# KM by APOE4 carrier status
df_km <- df %>% filter(!is.na(APOE_e4_count))
df_km$APOE4_status <- ifelse(df_km$APOE_e4_count > 0, "APOE4 Carrier", "Non-carrier")

surv_obj <- Surv(time=df_km$cox_time, event=df_km$cox_event)
km_fit   <- survfit(surv_obj ~ APOE4_status, data=df_km)

p_km_apoe <- ggsurvplot(
    km_fit, data=df_km,
    palette=c("#d73027","#4393c3"),
    conf.int=TRUE,
    pval=TRUE, pval.method=TRUE,
    risk.table=TRUE,
    risk.table.height=0.25,
    xlab="Age (years)",
    ylab="AD/MCI-free Probability",
    title="Age of Onset by APOE4 Status",
    legend.title="APOE4 Status",
    ggtheme=theme_bw(base_size=11),
    surv.median.line="hv"
)

# KM by sex
km_fit_sex <- survfit(surv_obj ~ sex_bin, data=df)
p_km_sex <- ggsurvplot(
    km_fit_sex, data=df,
    palette=c("#d7191c","#2c7bb6"),
    conf.int=TRUE, pval=TRUE,
    risk.table=TRUE, risk.table.height=0.25,
    xlab="Age (years)", ylab="AD/MCI-free Probability",
    title="Age of Onset by Sex",
    legend.labs=c("Female","Male"),
    legend.title="Sex",
    ggtheme=theme_bw(base_size=11)
)

# KM by ancestry region
df_region <- df %>% filter(region != "Missing")
km_fit_region <- survfit(Surv(cox_time, cox_event) ~ region, data=df_region)
p_km_region <- ggsurvplot(
    km_fit_region, data=df_region,
    palette=c("#2166ac","#d6604d","#4dac26"),
    conf.int=FALSE, pval=TRUE,
    risk.table=TRUE, risk.table.height=0.25,
    xlab="Age (years)", ylab="AD/MCI-free Probability",
    title="Age of Onset by Region of Origin",
    ggtheme=theme_bw(base_size=11)
)

pdf(file.path(FIG_DIR, "fig_km_curves.pdf"), width=10, height=7)
print(p_km_apoe)
print(p_km_sex)
print(p_km_region)
dev.off()
cat("  Written: fig_km_curves.pdf\n")

# ─── COX NULL MODEL ───────────────────────────────────────────────────────────
cat("\n--- Fitting null Cox model (covariates only) ---\n")

cox_null <- coxph(
    Surv(cox_time, cox_event) ~
        sex_bin + age + APOE_e4_count +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
    data=df
)

cox_summary <- summary(cox_null)
cat("  Null Cox model summary:\n")
print(cox_summary$coefficients)
cat(sprintf("  Concordance: %.3f (se=%.3f)\n",
            cox_summary$concordance[1], cox_summary$concordance[2]))
cat(sprintf("  Likelihood ratio test p: %.2e\n", cox_summary$logtest[3]))

checkpoint("Cox null model fitted",    !is.null(cox_null))
checkpoint("Cox concordance > 0.5",    cox_summary$concordance[1] > 0.5)

# Test proportional hazards assumption
cat("\n--- Testing proportional hazards assumption ---\n")
ph_test <- cox.zph(cox_null)
print(ph_test)
if (any(ph_test$table[,"p"] < 0.05 & rownames(ph_test$table) != "GLOBAL")) {
    cat("  ⚠️  WARNING: Some covariates violate PH assumption\n")
    cat("      Consider time-varying covariates or stratified Cox model\n")
} else {
    cat("  ✅ PH assumption satisfied for all covariates\n")
}

pdf(file.path(FIG_DIR, "fig_cox_forest.pdf"), width=8, height=6)
ggforest(cox_null, data=df, main="Cox PH Model — Covariates")
dev.off()
cat("  Written: fig_cox_forest.pdf\n")

# ─── PLINK2 COX GWAS COMMAND GENERATOR ────────────────────────────────────────
cat("\n--- Writing PLINK2 Cox GWAS command ---\n")

# Phenotype file in PLINK2 format for Cox
cox_pheno_plink <- df %>%
    select(FID, IID, cox_time, cox_event) %>%
    rename(TIME=cox_time, EVENT=cox_event)
fwrite(cox_pheno_plink, file.path(GWAS_DIR, "cox_pheno_plink.txt"), sep=" ")

# Covariate file
cox_covar_plink <- df %>%
    select(FID, IID, sex_bin, age, APOE_e4_count,
           all_of(paste0("PC", 1:10)),
           starts_with("site_")) %>%
    rename(sex=sex_bin)
fwrite(cox_covar_plink, file.path(GWAS_DIR, "cox_covar_plink.txt"),
       sep=" ", na="NA")

# Generate PLINK2 Cox GWAS bash command
cox_cmd <- sprintf(
'#!/bin/bash
# PLINK2 Cox PH GWAS
source ${PIPELINE_DIR}/pipeline_config.sh

GWAS_DIR="${PIPELINE_DIR}/04_gwas/cox"

%s \\
    --bfile %s \\
    --pheno %s \\
    --pheno-name EVENT \\
    --covar %s \\
    --1 \\
    --covar-variance-standardize \\
    --cox TIME \\
    --out ${GWAS_DIR}/plink2_cox_aoo \\
    --threads ${THREADS} \\
    --memory ${MEMORY}
',
    Sys.getenv("PLINK2"),
    file.path(Sys.getenv("PIPELINE_DIR"), "02_genotype_qc", "step2b_sampleQC"),
    file.path(GWAS_DIR, "cox_pheno_plink.txt"),
    file.path(GWAS_DIR, "cox_covar_plink.txt")
)

writeLines(cox_cmd, file.path(GWAS_DIR, "run_plink2_cox.sh"))
Sys.chmod(file.path(GWAS_DIR, "run_plink2_cox.sh"), "0755")
cat("  Written: run_plink2_cox.sh\n")
cat("  Run this script after confirming PLINK2 Cox support in your version\n")

cat("\n")
cat("============================================================\n")
cat(sprintf(" STEP 4b COMPLETE ✅  [%s]\n", Sys.time()))
cat(" KM curves:  09_figures/fig_km_curves.pdf\n")
cat(" Cox forest: 09_figures/fig_cox_forest.pdf\n")
cat(" Next: Run step05_postGWAS.py after both GWAS complete\n")
cat("============================================================\n")

sink()
sink(type="message")
