#!/usr/bin/env python3
# =============================================================================
# STEP 1: Phenotype Preparation & QC
# AOO GWAS Pipeline
# =============================================================================
# PURPOSE:
#   - Load raw phenotype CSV
#   - Clean and recode all variables
#   - Create analysis-ready phenotype files for REGENIE and R Cox model
#   - Generate Table 1 (demographic summary) for paper
#   - Run comprehensive QC checkpoints at each stage
#
# OUTPUT FILES:
#   01_phenotype/pheno_aoo_regenie.txt    — REGENIE quantitative phenotype
#   01_phenotype/pheno_cox.txt            — Survival analysis phenotype
#   01_phenotype/covariates_noPCs.txt     — Covariates (PCs added after Step 3)
#   01_phenotype/sample_qc_log.txt        — QC summary log
#   01_phenotype/table1_demographics.csv  — Table 1 for paper
# =============================================================================

import pandas as pd
import numpy as np
import os, sys, logging
from datetime import datetime

# ─── CONFIG ──────────────────────────────────────────────────────────────────
PHENO_CSV = os.environ.get("PHENO_CSV", "phenotype.csv")
OUT_DIR   = os.path.join(os.environ.get("PIPELINE_DIR", "."), "01_phenotype")
os.makedirs(OUT_DIR, exist_ok=True)

# ─── LOGGING ─────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(f"{OUT_DIR}/sample_qc_log.txt"),
        logging.StreamHandler(sys.stdout)
    ]
)
log = logging.getLogger()

def checkpoint(name, condition, df_before=None, df_after=None):
    """Validate a pipeline checkpoint. Exits on failure."""
    if condition:
        msg = f"✅  CHECKPOINT PASSED: {name}"
        if df_before is not None and df_after is not None:
            msg += f"  [{len(df_before)} → {len(df_after)} samples]"
        log.info(msg)
    else:
        log.error(f"❌  CHECKPOINT FAILED: {name}")
        sys.exit(1)

# ─── STEP 1.1: LOAD DATA ─────────────────────────────────────────────────────
log.info("=" * 60)
log.info(" STEP 1: PHENOTYPE PREPARATION")
log.info("=" * 60)
log.info(f"Loading: {PHENO_CSV}")

df_raw = pd.read_csv(PHENO_CSV)
checkpoint("File loaded", len(df_raw) > 0)

# Update this list to match the column names in your phenotype file
required_cols = ['SAMPLE_ID', 'site', 'sex', 'age_at_subject',
                 'AOO', 'CDX', 'APOE', 'country_of_birth']
checkpoint("Expected columns present",
           all(c in df_raw.columns for c in required_cols))
log.info(f"  Raw shape: {df_raw.shape}")
df = df_raw.copy()

# ─── STEP 1.2: CLEAN NUMERIC FIELDS ─────────────────────────────────────────
log.info("\n--- Cleaning numeric fields ---")

# Age at assessment
df['age'] = pd.to_numeric(df['age_at_subject'].replace('.', np.nan), errors='coerce')

# Age of onset (AOO) — replace '.' sentinel with NaN
df['AOO_raw'] = df['AOO'].replace('.', np.nan)
df['AOO']     = pd.to_numeric(df['AOO_raw'], errors='coerce')

# Flag biologically implausible AOO values for late-onset AD
df.loc[df['AOO'] < 40, 'AOO'] = np.nan
df.loc[df['AOO'] > 95, 'AOO'] = np.nan

log.info(f"  Age: {df['age'].notna().sum()} non-missing, "
         f"range {df['age'].min():.1f}–{df['age'].max():.1f}")
log.info(f"  AOO: {df['AOO'].notna().sum()} non-missing (implausible values removed), "
         f"range {df['AOO'].min():.1f}–{df['AOO'].max():.1f}")

# ─── STEP 1.3: RECODE CATEGORICAL VARIABLES ──────────────────────────────────
log.info("\n--- Recoding categorical variables ---")

# Sex: F=0, M=1 (standard genetic convention)
sex_map = {'F': 0, 'M': 1}
df['sex_bin'] = df['sex'].map(sex_map)
checkpoint("Sex coding valid", df['sex_bin'].notna().all())

# Diagnosis (CDX)
df['CDX_clean'] = df['CDX'].str.strip()
dx_map = {
    'AD':  2,   # Alzheimer's Disease (case)
    'MCI': 1,   # Mild Cognitive Impairment
    'NCI': 0,   # No Cognitive Impairment (control)
    '0':   0    # Unaffected, treat as control
}
df['CDX_num'] = df['CDX_clean'].map(dx_map)

# Binary case/control: AD+MCI = case, NCI = control
df['case_control'] = np.where(df['CDX_clean'].isin(['AD','MCI']), 2,
                     np.where(df['CDX_clean'] == 'NCI', 1, np.nan))

log.info(f"  AD:  {(df['CDX_clean']=='AD').sum()}")
log.info(f"  MCI: {(df['CDX_clean']=='MCI').sum()}")
log.info(f"  NCI: {(df['CDX_clean']=='NCI').sum()}")
log.info(f"  Other/unknown: {df['CDX_num'].isna().sum()}")

# APOE: parse XX format → ε4 and ε2 allele counts
def parse_apoe(val):
    """Convert APOE genotype code (e.g. 33, 34, 44) to ε4 allele count."""
    if pd.isna(val):
        return np.nan, np.nan, np.nan
    code = str(int(val))
    if len(code) != 2:
        return np.nan, np.nan, np.nan
    a1, a2 = int(code[0]), int(code[1])
    e4 = (a1 == 4) + (a2 == 4)
    e2 = (a1 == 2) + (a2 == 2)
    e4_carrier = 1 if e4 > 0 else 0
    return e4, e2, e4_carrier

apoe_parsed = df['APOE'].apply(parse_apoe)
df['APOE_e4_count']   = [x[0] for x in apoe_parsed]
df['APOE_e2_count']   = [x[1] for x in apoe_parsed]
df['APOE_e4_carrier'] = [x[2] for x in apoe_parsed]

log.info(f"  APOE4 carriers: {df['APOE_e4_carrier'].sum():.0f} / "
         f"{df['APOE_e4_carrier'].notna().sum():.0f}")

# Site dummy encoding for covariates
site_dummies = pd.get_dummies(df['site'], prefix='site', drop_first=True)
df = pd.concat([df, site_dummies], axis=1)
site_cols = site_dummies.columns.tolist()
log.info(f"  Site dummies: {site_cols}")

# Country-of-birth region grouping
# Update these sets to match your cohort's country codes
caribbean = {'PRI', 'DOM', 'CUB'}
df['region'] = df['country_of_birth'].apply(
    lambda x: 'Caribbean' if x in caribbean
    else ('USA' if x == 'USA'
    else ('Missing' if pd.isna(x) or x == 'Missing' else 'LatAm')))

# ─── STEP 1.4: DEFINE ANALYSIS SUBSETS ───────────────────────────────────────
log.info("\n--- Defining analysis subsets ---")

# Primary AOO analysis: cases (AD/MCI) with valid AOO
df_aoo_cases = df[
    df['CDX_clean'].isin(['AD', 'MCI']) &
    df['AOO'].notna() &
    df['age'].notna() &
    df['sex_bin'].notna()
].copy()
log.info(f"  AOO linear model (cases with valid AOO): n={len(df_aoo_cases)}")

# Cox model: all cases + controls with age and diagnosis
df_cox = df[
    df['CDX_clean'].isin(['AD', 'MCI', 'NCI']) &
    df['age'].notna() &
    df['sex_bin'].notna()
].copy()

# Cox time: AOO if case, age at assessment if control (censored)
df_cox['cox_time']  = np.where(df_cox['CDX_clean'].isin(['AD','MCI']),
                                df_cox['AOO'].fillna(df_cox['age']),
                                df_cox['age'])
df_cox['cox_event'] = np.where(df_cox['CDX_clean'].isin(['AD','MCI']), 1, 0)
log.info(f"  Cox model dataset: n={len(df_cox)} ({df_cox['cox_event'].sum()} events)")

checkpoint("Sufficient cases for AOO GWAS", len(df_aoo_cases) >= 200)
checkpoint("Sufficient samples for Cox GWAS", len(df_cox) >= 500)

# ─── STEP 1.5: PLINK/REGENIE ID FORMAT ───────────────────────────────────────
log.info("\n--- Preparing IDs for PLINK/REGENIE format ---")
# REGENIE expects FID IID — use SAMPLE_ID as both
# IMPORTANT: These must exactly match the IDs in your .fam file
df_aoo_cases['FID'] = df_aoo_cases['SAMPLE_ID']
df_aoo_cases['IID'] = df_aoo_cases['SAMPLE_ID']
df_cox['FID'] = df_cox['SAMPLE_ID']
df_cox['IID'] = df_cox['SAMPLE_ID']
df['FID'] = df['SAMPLE_ID']
df['IID'] = df['SAMPLE_ID']

# ─── STEP 1.6: WRITE REGENIE PHENOTYPE FILE ──────────────────────────────────
log.info("\n--- Writing REGENIE phenotype file ---")
# Format: FID IID PHENO1 ... (space-delimited, NA for missing)
pheno_regenie = df[['FID','IID']].copy()
pheno_regenie['AOO_quant'] = np.where(
    df['CDX_clean'].isin(['AD','MCI']), df['AOO'], np.nan)

# Z-score normalize AOO
aoo_mean = pheno_regenie['AOO_quant'].mean()
aoo_std  = pheno_regenie['AOO_quant'].std()
pheno_regenie['AOO_zscore'] = (pheno_regenie['AOO_quant'] - aoo_mean) / aoo_std

pheno_regenie.to_csv(f"{OUT_DIR}/pheno_aoo_regenie.txt",
                     sep=' ', index=False, na_rep='NA')
log.info(f"  Written: pheno_aoo_regenie.txt  "
         f"(n non-NA = {pheno_regenie['AOO_quant'].notna().sum()})")

# ─── STEP 1.7: WRITE COVARIATES FILE ─────────────────────────────────────────
log.info("--- Writing covariates file ---")
# Covariates: sex, age, APOE_e4_count, site dummies
# PCs will be added after Step 3 (PCA)
cov_cols = ['FID','IID','sex_bin','age','APOE_e4_count'] + site_cols
covariates = df[cov_cols].copy()
covariates.columns = ['FID','IID','sex','age','APOE_e4_count'] + site_cols
covariates.to_csv(f"{OUT_DIR}/covariates_noPCs.txt",
                  sep=' ', index=False, na_rep='NA')
log.info(f"  Written: covariates_noPCs.txt  (PCs to be added after Step 3)")

# ─── STEP 1.8: WRITE COX MODEL PHENOTYPE FILE ────────────────────────────────
log.info("--- Writing Cox model phenotype file ---")
cox_out = df_cox[['FID','IID','cox_time','cox_event','sex_bin','age',
                   'APOE_e4_count','CDX_clean','region'] + site_cols].copy()
cox_out.to_csv(f"{OUT_DIR}/pheno_cox.txt", sep='\t', index=False, na_rep='NA')
log.info(f"  Written: pheno_cox.txt  "
         f"(n={len(cox_out)}, events={int(cox_out['cox_event'].sum())})")

# ─── STEP 1.9: TABLE 1 — DEMOGRAPHICS ────────────────────────────────────────
log.info("\n--- Generating Table 1 ---")

def fmt_mean_sd(series):
    s = pd.to_numeric(series, errors='coerce').dropna()
    return f"{s.mean():.1f} ± {s.std():.1f}"

def fmt_n_pct(count, total):
    return f"{count} ({100*count/total:.1f}%)"

groups = {
    'NCI': df[df['CDX_clean']=='NCI'],
    'MCI': df[df['CDX_clean']=='MCI'],
    'AD':  df[df['CDX_clean']=='AD']
}

table1_rows = []
for grp, gdf in groups.items():
    n = len(gdf)
    row = {
        'Group': grp,
        'N': n,
        'Age (mean±SD)': fmt_mean_sd(gdf['age']),
        'Female (n,%)': fmt_n_pct((gdf['sex']=='F').sum(), n),
        'AOO (mean±SD)': fmt_mean_sd(gdf['AOO']),
        'APOE4 carrier (n,%)': fmt_n_pct(gdf['APOE_e4_carrier'].sum(), n),
        'APOE missing (n)': int(gdf['APOE'].isna().sum()),
    }
    table1_rows.append(row)

table1 = pd.DataFrame(table1_rows)
table1.to_csv(f"{OUT_DIR}/table1_demographics.csv", index=False)
log.info(f"  Written: table1_demographics.csv")
print("\nTable 1 preview:")
print(table1.to_string(index=False))

# ─── FINAL CHECKPOINTS ───────────────────────────────────────────────────────
log.info("\n--- Final checkpoints ---")
checkpoint("REGENIE phenotype file written",
           os.path.exists(f"{OUT_DIR}/pheno_aoo_regenie.txt"))
checkpoint("Covariates file written",
           os.path.exists(f"{OUT_DIR}/covariates_noPCs.txt"))
checkpoint("Cox phenotype file written",
           os.path.exists(f"{OUT_DIR}/pheno_cox.txt"))
checkpoint("Table 1 written",
           os.path.exists(f"{OUT_DIR}/table1_demographics.csv"))
checkpoint("No duplicate SAMPLE_IDs", df['SAMPLE_ID'].nunique() == len(df))
checkpoint("AOO non-missing in expected range",
           df['AOO'].dropna().between(40, 96).all())

log.info("\n" + "=" * 60)
log.info(" STEP 1 COMPLETE ✅")
log.info(f" Outputs in: {OUT_DIR}")
log.info(" Next: Run step02_genotype_qc.sh")
log.info("=" * 60)
