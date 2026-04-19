#!/bin/bash
# =============================================================================
# STEP 4a: GWAS — Linear Regression on AOO (REGENIE)
# AOO GWAS Pipeline
# =============================================================================
# PURPOSE:
#   Run genome-wide association of age-of-onset (quantitative) using REGENIE.
#   REGENIE is preferred over PLINK2 for:
#     - Handling cryptically related samples via whole-genome regression
#     - Better calibration in admixed/structured populations
#     - Efficient handling of large imputed datasets
#
# MODEL: AOO_zscore ~ SNP + sex + age + APOE_e4_count + PC1..PC10 + site
#        (cases: AD and MCI with valid AOO only)
#
# TWO-STEP REGENIE APPROACH:
#   Step 1 (Ridge regression): fits null model on common variants
#   Step 2 (Score test):       tests each imputed variant against null residuals
# =============================================================================

set -euo pipefail
source ${PIPELINE_DIR}/pipeline_config.sh

GWAS_DIR="${PIPELINE_DIR}/04_gwas/regenie"
LOG="${PIPELINE_DIR}/logs/step04a_regenie.log"
mkdir -p ${GWAS_DIR}
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo " STEP 4a: REGENIE GWAS  [$(date)]"
echo "============================================================"

checkpoint() {
    [ -f "$2" ] && echo "  ✅ CHECKPOINT: $1" || \
        { echo "  ❌ CHECKPOINT FAILED: $1 — $2 not found"; exit 1; }
}

# ─── INPUT FILES ─────────────────────────────────────────────────────────────
PHENO_FILE="${PIPELINE_DIR}/01_phenotype/pheno_aoo_regenie.txt"
COVAR_FILE="${PIPELINE_DIR}/01_phenotype/covariates_withPCs.txt"
GENO_QC="${PIPELINE_DIR}/02_genotype_qc/step2b_sampleQC"

checkpoint "Phenotype file"         "${PHENO_FILE}"
checkpoint "Covariates with PCs"    "${COVAR_FILE}"
checkpoint "QC'd genotype .bed"     "${GENO_QC}.bed"

# ─── REGENIE STEP 1: Null model (whole-genome ridge regression) ───────────────
echo ""
echo "--- REGENIE Step 1: Fitting null model ---"
echo "    Accounts for population structure and relatedness."
echo "    Expected runtime: 30–60 min"

${REGENIE} \
    --step 1 \
    --bed ${GENO_QC} \
    --phenoFile ${PHENO_FILE} \
    --covarFile ${COVAR_FILE} \
    --phenoColList AOO_zscore \
    --covarColList sex,age,APOE_e4_count,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix ${GWAS_DIR}/regenie_tmp \
    --out ${GWAS_DIR}/regenie_step1 \
    --threads ${THREADS} \
    --gz

checkpoint "REGENIE Step 1 predictions" "${GWAS_DIR}/regenie_step1_pred.list"
echo "  ✅ REGENIE Step 1 complete"

# ─── REGENIE STEP 2: Association testing ─────────────────────────────────────
echo ""
echo "--- REGENIE Step 2: Association testing ---"
echo "    Tests each imputed variant against null model residuals."
echo "    Expected runtime: 1–4 hrs"

${REGENIE} \
    --step 2 \
    --bed ${GENO_QC} \
    --phenoFile ${PHENO_FILE} \
    --covarFile ${COVAR_FILE} \
    --phenoColList AOO_zscore \
    --covarColList sex,age,APOE_e4_count,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --pred ${GWAS_DIR}/regenie_step1_pred.list \
    --bsize 400 \
    --qt \
    --gz \
    --out ${GWAS_DIR}/regenie_aoo \
    --threads ${THREADS}

checkpoint "REGENIE Step 2 results" "${GWAS_DIR}/regenie_aoo_AOO_zscore.regenie.gz"
echo "  ✅ REGENIE Step 2 complete"

# ─── POST-REGENIE: RESULTS CHECK ──────────────────────────────────────────────
echo ""
echo "--- Checking GWAS results ---"

python3 - <<'PYCHECK'
import pandas as pd, numpy as np, os, gzip

gwas_dir = os.environ['GWAS_DIR']
res_file = f"{gwas_dir}/regenie_aoo_AOO_zscore.regenie.gz"

df = pd.read_csv(res_file, sep=' ', compression='gzip')
print(f"  Total variants tested: {len(df):,}")
print(f"  Columns: {list(df.columns)}")

# REGENIE output columns: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P
df['P'] = 10 ** (-df['LOG10P'])

# Genomic inflation factor (lambda GC)
chisq = df['CHISQ'].dropna()
lambda_gc = np.median(chisq) / 0.4549
print(f"\n  λ GC (genomic inflation factor): {lambda_gc:.4f}")
if lambda_gc > 1.15:
    print(f"  ⚠️  WARNING: λ GC = {lambda_gc:.3f} — may indicate residual stratification")
    print(f"     Consider adding more PCs or checking covariate specification")
elif lambda_gc < 0.95:
    print(f"  ⚠️  WARNING: λ GC = {lambda_gc:.3f} — unexpectedly low")
else:
    print(f"  ✅ λ GC acceptable")

# Significant hits
gws      = df[df['LOG10P'] >= 7.301]   # p < 5e-8
suggestive = df[(df['LOG10P'] >= 5.301) & (df['LOG10P'] < 7.301)]  # p < 5e-6
print(f"\n  Genome-wide significant (p<5e-8): {len(gws)} variants")
print(f"  Suggestive (p<5e-6):              {len(suggestive)} variants")

if len(gws) > 0:
    print(f"\n  Top 10 hits:")
    top = df.nlargest(10, 'LOG10P')[['CHROM','GENPOS','ID','A1FREQ','BETA','SE','P']]
    print(top.to_string(index=False))

df.to_csv(f"{gwas_dir}/regenie_aoo_clean.txt.gz", sep='\t', index=False, compression='gzip')
print(f"\n  Saved: regenie_aoo_clean.txt.gz")
PYCHECK

echo ""
echo "============================================================"
echo " STEP 4a COMPLETE ✅  [$(date)]"
echo " Next: Run step04b_cox_gwas.R"
echo "============================================================"

# ─── SLURM SUBMISSION TEMPLATE ───────────────────────────────────────────────
cat > ${GWAS_DIR}/submit_regenie.slurm << 'SLURM'
#!/bin/bash
#SBATCH --job-name=aoo_gwas_regenie
#SBATCH --output=logs/regenie_%j.log
#SBATCH --error=logs/regenie_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --partition=compute

source pipeline_config.sh
bash 00_scripts/step04a_gwas_regenie.sh
SLURM
echo "  SLURM script: ${GWAS_DIR}/submit_regenie.slurm"
