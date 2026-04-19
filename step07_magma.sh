#!/bin/bash
# =============================================================================
# STEP 7: MAGMA — Gene-Based & Pathway Analysis
# AOO GWAS Pipeline
# =============================================================================
# PURPOSE:
#   Aggregate SNP-level p-values to gene-level tests, then test for
#   enrichment in biological pathways (MSigDB, GO terms).
#
# REQUIRES (download before running):
#   - MAGMA LD reference panel: https://ctg.cncr.nl/software/magma
#   - Gene location file (NCBI37.3.gene.loc)
#   - Gene sets (MSigDB .gmt file): https://www.gsea-msigdb.org/gsea/msigdb/
# =============================================================================

set -euo pipefail
source ${PIPELINE_DIR}/pipeline_config.sh

MAGMA_DIR="${PIPELINE_DIR}/07_magma"
mkdir -p ${MAGMA_DIR}
LOG="${PIPELINE_DIR}/logs/step07_magma.log"
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo " STEP 7: MAGMA GENE & PATHWAY ANALYSIS  [$(date)]"
echo "============================================================"

# ─── REFERENCE FILES ─────────────────────────────────────────────────────────
MAGMA_REF_DIR="${MAGMA_DIR}/reference"
mkdir -p ${MAGMA_REF_DIR}

# Use a population-matched LD reference if available
# For admixed populations, European LD panel is an approximation —
# state this limitation in any manuscript methods section
MAGMA_REF="${MAGMA_REF_DIR}/g1000_eur"
GENE_LOC="${MAGMA_REF_DIR}/NCBI37.3.gene.loc"
GENE_SETS="${MAGMA_REF_DIR}/msigdb_v7.2.entrez.gmt"

if [ ! -f "${MAGMA_REF}.bed" ]; then
    echo ""
    echo "  ⚠️  MAGMA reference files not found."
    echo "  Download and extract to ${MAGMA_REF_DIR}:"
    echo "    wget https://ctg.cncr.nl/software/MAGMA/aux_files/g1000_eur.zip"
    echo "    wget https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.gene.loc.zip"
    echo "  Then re-run this step."
    exit 0
fi

# ─── 7.1: PREPARE SNP P-VALUE FILE ───────────────────────────────────────────
echo ""
echo "--- 7.1: Preparing SNP location file ---"

python3 - <<'PYPREP'
import pandas as pd, os

gwas_dir  = os.environ['PIPELINE_DIR'] + "/04_gwas"
magma_dir = os.environ['MAGMA_DIR']

df = pd.read_csv(f"{gwas_dir}/regenie/regenie_aoo_clean.txt.gz",
                 sep='\t', compression='gzip',
                 usecols=['CHR','POS','SNP','P','N'])
df = df.rename(columns={'POS':'BP'})
df = df.dropna(subset=['P'])
df['P'] = df['P'].clip(lower=1e-300)

# MAGMA SNP pval file format: SNP CHR BP P N
df[['SNP','CHR','BP','P','N']].to_csv(f"{magma_dir}/snp_pvals.txt",
                                       sep=' ', index=False)
print(f"  Written: snp_pvals.txt  ({len(df):,} SNPs)")
PYPREP

# ─── 7.2: ANNOTATION STEP ────────────────────────────────────────────────────
echo ""
echo "--- 7.2: Annotating SNPs to genes (±10 kb window) ---"

${MAGMA} \
    --annotate window=10,10 \
    --snp-loc ${MAGMA_DIR}/snp_pvals.txt \
    --gene-loc ${GENE_LOC} \
    --out ${MAGMA_DIR}/magma_annot

[ -f "${MAGMA_DIR}/magma_annot.genes.annot" ] && \
    echo "  ✅ Annotation complete" || \
    { echo "  ❌ Annotation failed"; exit 1; }

# ─── 7.3: GENE-BASED ASSOCIATION TEST ────────────────────────────────────────
echo ""
echo "--- 7.3: Gene-based association test ---"

${MAGMA} \
    --bfile ${MAGMA_REF} \
    --gene-annot ${MAGMA_DIR}/magma_annot.genes.annot \
    --pval ${MAGMA_DIR}/snp_pvals.txt use=SNP,P ncol=N \
    --gene-model snp-wise=mean \
    --out ${MAGMA_DIR}/magma_gene

[ -f "${MAGMA_DIR}/magma_gene.genes.out" ] && \
    echo "  ✅ Gene-based test complete" || \
    { echo "  ❌ Gene-based test failed"; exit 1; }

echo ""
echo "  Top 20 gene-based associations:"
sort -k9 -g ${MAGMA_DIR}/magma_gene.genes.out | head -21

# Count significant genes (Bonferroni: 0.05 / ~20,000 genes)
N_SIG_GENES=$(awk 'NR>1 && $9 < 2.5e-6 {count++} END {print count+0}' \
              ${MAGMA_DIR}/magma_gene.genes.out)
echo "  Gene-level significant (p<2.5e-6): ${N_SIG_GENES}"

# ─── 7.4: PATHWAY ANALYSIS ───────────────────────────────────────────────────
echo ""
echo "--- 7.4: Pathway/gene-set enrichment ---"

if [ -f "${GENE_SETS}" ]; then
    ${MAGMA} \
        --gene-results ${MAGMA_DIR}/magma_gene.genes.raw \
        --set-annot ${GENE_SETS} \
        --out ${MAGMA_DIR}/magma_pathways

    [ -f "${MAGMA_DIR}/magma_pathways.gsa.out" ] && \
        echo "  ✅ Pathway analysis complete" || \
        echo "  ⚠️  Pathway analysis failed"

    echo ""
    echo "  Top 20 enriched pathways:"
    sort -k8 -g ${MAGMA_DIR}/magma_pathways.gsa.out | head -21
else
    echo "  ⚠️  Gene set file not found — download MSigDB .gmt from:"
    echo "      https://www.gsea-msigdb.org/gsea/msigdb/"
fi

echo ""
echo "============================================================"
echo " STEP 7 COMPLETE ✅  [$(date)]"
echo "============================================================"
