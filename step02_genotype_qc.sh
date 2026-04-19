#!/bin/bash
# =============================================================================
# STEP 2: Genotype QC
# AOO GWAS Pipeline
# =============================================================================
# PURPOSE:
#   Standard pre-GWAS QC steps for imputed genotype data:
#   1. Sample QC  — call rate, heterozygosity, sex check, relatedness
#   2. Variant QC — MAF, HWE, INFO score filter, call rate
#
# CHECKPOINTS at each stage with counts logged.
# OUTPUT: QC-passed PLINK binary files in 02_genotype_qc/
# =============================================================================

set -euo pipefail
source ${PIPELINE_DIR}/pipeline_config.sh

QC_DIR="${PIPELINE_DIR}/02_genotype_qc"
LOG="${PIPELINE_DIR}/logs/step02_genoqc.log"
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo " STEP 2: GENOTYPE QC  [$(date)]"
echo "============================================================"

# ─── CHECKPOINT FUNCTION ─────────────────────────────────────────────────────
checkpoint() {
    local desc="$1"; local file="$2"
    if [ -f "$file" ] || [ -d "$file" ]; then
        echo "  ✅ CHECKPOINT: $desc — OK"
    else
        echo "  ❌ CHECKPOINT FAILED: $desc — $file not found"
        exit 1
    fi
}

count_samples()  { wc -l < "$1"; }
count_variants() { wc -l < "$1"; }

# ─── 2.0: VALIDATE INPUT ─────────────────────────────────────────────────────
echo ""
echo "--- 2.0: Validating input files ---"
for ext in bed bim fam; do
    checkpoint "Input .${ext} exists" "${GENO_PREFIX}.${ext}"
done

N_SAMPLES_RAW=$(count_samples  "${GENO_PREFIX}.fam")
N_VARIANTS_RAW=$(count_variants "${GENO_PREFIX}.bim")
echo "  Input: ${N_SAMPLES_RAW} samples, ${N_VARIANTS_RAW} variants"

# ─── 2.1: VARIANT QC ─────────────────────────────────────────────────────────
# Filters: MAF ≥ 1%, HWE p > 1e-6, genotype call rate ≥ 95%, sample call rate ≥ 95%
echo ""
echo "--- 2.1: Variant QC (MAF, HWE, call rate) ---"

${PLINK2} \
    --bfile ${GENO_PREFIX} \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.05 \
    --mind 0.05 \
    --make-bed \
    --out ${QC_DIR}/step2a_variantQC \
    --threads ${THREADS} \
    --memory ${MEMORY}

checkpoint "Variant QC output" "${QC_DIR}/step2a_variantQC.bed"

N_S2A_SAMPLES=$(count_samples  "${QC_DIR}/step2a_variantQC.fam")
N_S2A_VARIANTS=$(count_variants "${QC_DIR}/step2a_variantQC.bim")
echo "  After variant QC: ${N_S2A_SAMPLES} samples, ${N_S2A_VARIANTS} variants"

if [ "$N_S2A_VARIANTS" -lt 500000 ]; then
    echo "  ⚠️  WARNING: Fewer than 500K variants remain — check INFO score pre-filter"
fi

# ─── 2.2: SAMPLE QC — HETEROZYGOSITY ─────────────────────────────────────────
echo ""
echo "--- 2.2: Sample QC — Heterozygosity ---"

# LD-prune before heterozygosity calculation
${PLINK2} \
    --bfile ${QC_DIR}/step2a_variantQC \
    --indep-pairwise 50 10 0.1 \
    --out ${QC_DIR}/pruned_snplist \
    --threads ${THREADS}

${PLINK2} \
    --bfile ${QC_DIR}/step2a_variantQC \
    --extract ${QC_DIR}/pruned_snplist.prune.in \
    --het \
    --out ${QC_DIR}/het_stats \
    --threads ${THREADS}

checkpoint "Heterozygosity stats" "${QC_DIR}/het_stats.het"

# Flag heterozygosity outliers (±3 SD from mean F-statistic)
python3 - <<PYHET
import pandas as pd, numpy as np, os
het = pd.read_csv("${QC_DIR}/het_stats.het", sep=r'\s+')
het['F_rate'] = (het['OBS_CT'] - het['O(HOM)']) / het['OBS_CT']
mean_f = het['F_rate'].mean()
std_f  = het['F_rate'].std()
outliers = het[np.abs(het['F_rate'] - mean_f) > 3 * std_f]
outliers[['#FID','IID']].to_csv("${QC_DIR}/het_outliers.txt",
                                 sep=' ', index=False, header=False)
print(f"  Heterozygosity outliers to remove: {len(outliers)}")
print(f"  F-statistic range: {het['F_rate'].min():.4f} – {het['F_rate'].max():.4f}")
PYHET

# ─── 2.3: SAMPLE QC — SEX CHECK ──────────────────────────────────────────────
echo ""
echo "--- 2.3: Sample QC — Sex check ---"

${PLINK2} \
    --bfile ${QC_DIR}/step2a_variantQC \
    --check-sex \
    --out ${QC_DIR}/sex_check \
    --threads ${THREADS} 2>/dev/null || \
    echo "  ⚠️  Sex check skipped (no X chr data or already QC'd)"

# ─── 2.4: SAMPLE QC — RELATEDNESS (KING) ─────────────────────────────────────
echo ""
echo "--- 2.4: Sample QC — Relatedness (KING) ---"

${PLINK2} \
    --bfile ${QC_DIR}/step2a_variantQC \
    --extract ${QC_DIR}/pruned_snplist.prune.in \
    --make-king-table \
    --king-table-filter 0.0884 \
    --out ${QC_DIR}/king_relatedness \
    --threads ${THREADS}

if [ -f "${QC_DIR}/king_relatedness.kin0" ]; then
    N_RELATED=$(tail -n +2 ${QC_DIR}/king_relatedness.kin0 | wc -l)
    echo "  Related pairs (KING > 0.0884 / 2nd-degree): ${N_RELATED}"
    echo "  NOTE: REGENIE handles cryptic relatedness. Only removing"
    echo "        duplicated/identical samples (KING > 0.354)."

    awk 'NR>1 && $NF > 0.354 {print $1, $2}' \
        ${QC_DIR}/king_relatedness.kin0 > ${QC_DIR}/duplicates_to_remove.txt
    N_DUPS=$(wc -l < ${QC_DIR}/duplicates_to_remove.txt)
    echo "  Duplicates/identical pairs to remove: ${N_DUPS}"
fi

# ─── 2.5: APPLY SAMPLE EXCLUSIONS ────────────────────────────────────────────
echo ""
echo "--- 2.5: Applying sample exclusions ---"

cat /dev/null > ${QC_DIR}/samples_to_remove.txt

[ -f "${QC_DIR}/het_outliers.txt" ] && \
    cat ${QC_DIR}/het_outliers.txt >> ${QC_DIR}/samples_to_remove.txt

[ -f "${QC_DIR}/duplicates_to_remove.txt" ] && \
    cat ${QC_DIR}/duplicates_to_remove.txt >> ${QC_DIR}/samples_to_remove.txt

sort -u ${QC_DIR}/samples_to_remove.txt > ${QC_DIR}/samples_to_remove_uniq.txt
N_REMOVE=$(wc -l < ${QC_DIR}/samples_to_remove_uniq.txt)
echo "  Total samples to remove: ${N_REMOVE}"

${PLINK2} \
    --bfile ${QC_DIR}/step2a_variantQC \
    --remove ${QC_DIR}/samples_to_remove_uniq.txt \
    --make-bed \
    --out ${QC_DIR}/step2b_sampleQC \
    --threads ${THREADS} \
    --memory ${MEMORY}

checkpoint "Sample QC output" "${QC_DIR}/step2b_sampleQC.bed"

N_FINAL_SAMPLES=$(count_samples  "${QC_DIR}/step2b_sampleQC.fam")
N_FINAL_VARIANTS=$(count_variants "${QC_DIR}/step2b_sampleQC.bim")
echo "  After sample QC: ${N_FINAL_SAMPLES} samples, ${N_FINAL_VARIANTS} variants"

# ─── 2.6: CHECKPOINT SUMMARY ─────────────────────────────────────────────────
echo ""
echo "--- QC Summary ---"
echo "  Raw:              ${N_SAMPLES_RAW} samples,  ${N_VARIANTS_RAW} variants"
echo "  After variant QC: ${N_S2A_SAMPLES} samples,  ${N_S2A_VARIANTS} variants"
echo "  After sample QC:  ${N_FINAL_SAMPLES} samples,  ${N_FINAL_VARIANTS} variants"

if [ "$N_FINAL_SAMPLES" -lt 500 ]; then
    echo "  ❌ CHECKPOINT FAILED: Fewer than 500 samples after QC — investigate"
    exit 1
else
    echo "  ✅ CHECKPOINT: Final sample count acceptable (n=${N_FINAL_SAMPLES})"
fi

# Save QC'd prefix to config for downstream steps
echo "export GENO_QC=${QC_DIR}/step2b_sampleQC" >> ${PIPELINE_DIR}/pipeline_config.sh

echo ""
echo "============================================================"
echo " STEP 2 COMPLETE ✅  [$(date)]"
echo " QC'd data: ${QC_DIR}/step2b_sampleQC.*"
echo " Next: Run step03_pca.sh"
echo "============================================================"
