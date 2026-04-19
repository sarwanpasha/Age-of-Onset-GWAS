#!/bin/bash
# =============================================================================
# STEP 3: Population Stratification — PCA
# AOO GWAS Pipeline
# =============================================================================
# PURPOSE:
#   Compute principal components to control for population stratification.
#   For admixed cohorts, stratification correction is critical.
#
# APPROACH:
#   1. LD-prune QC'd genotypes
#   2. Exclude high-LD / inversion regions (chr6 MHC, chr8 inversion)
#   3. Compute top 20 PCs with PLINK2
#   4. Generate PCA plots colored by ancestry group and diagnosis
#   5. Append top 10 PCs to covariates file
#
# OUTPUT:
#   03_pca/pca.eigenvec               — PC vectors (IID + PC1..PC20)
#   03_pca/pca.eigenval               — Eigenvalues
#   01_phenotype/covariates_withPCs.txt — Full covariates for REGENIE
#   09_figures/fig_pca.pdf            — PCA scatter + scree plot
# =============================================================================

set -euo pipefail
source ${PIPELINE_DIR}/pipeline_config.sh

PCA_DIR="${PIPELINE_DIR}/03_pca"
FIG_DIR="${PIPELINE_DIR}/09_figures"
LOG="${PIPELINE_DIR}/logs/step03_pca.log"
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo " STEP 3: PCA  [$(date)]"
echo "============================================================"

mkdir -p ${PCA_DIR} ${FIG_DIR}

checkpoint() {
    [ -f "$2" ] && echo "  ✅ CHECKPOINT: $1" || \
        { echo "  ❌ CHECKPOINT FAILED: $1 — $2"; exit 1; }
}

# ─── 3.1: EXCLUDE HIGH-LD REGIONS ────────────────────────────────────────────
echo ""
echo "--- 3.1: Excluding high-LD regions ---"

# Standard exclusion regions (hg38 coordinates)
cat > ${PCA_DIR}/high_ld_regions.txt << 'HIGHLD'
6 25500000 33500000 MHC
8 8000000 12000000 inversion_8p23
HIGHLD

${PLINK2} \
    --bfile ${GENO_QC} \
    --exclude range ${PCA_DIR}/high_ld_regions.txt \
    --indep-pairwise 200 100 0.1 \
    --out ${PCA_DIR}/pca_prune \
    --threads ${THREADS}

checkpoint "LD prune list" "${PCA_DIR}/pca_prune.prune.in"
N_PCA_SNPS=$(wc -l < ${PCA_DIR}/pca_prune.prune.in)
echo "  SNPs for PCA: ${N_PCA_SNPS}"

if [ "$N_PCA_SNPS" -lt 50000 ]; then
    echo "  ⚠️  WARNING: Fewer than 50K SNPs for PCA — results may be unreliable"
fi

# ─── 3.2: RUN PCA ─────────────────────────────────────────────────────────────
echo ""
echo "--- 3.2: Running PCA (top 20 PCs) ---"

${PLINK2} \
    --bfile ${GENO_QC} \
    --extract ${PCA_DIR}/pca_prune.prune.in \
    --pca 20 approx \
    --out ${PCA_DIR}/pca \
    --threads ${THREADS} \
    --memory ${MEMORY}

checkpoint "PCA eigenvec" "${PCA_DIR}/pca.eigenvec"
checkpoint "PCA eigenval" "${PCA_DIR}/pca.eigenval"

# ─── 3.3: MERGE PCs INTO COVARIATES FILE ─────────────────────────────────────
echo ""
echo "--- 3.3: Merging PCs with covariates ---"

python3 - <<'PYMERGE'
import pandas as pd, os

pipeline_dir = os.environ['PIPELINE_DIR']
pca_dir      = os.environ['PCA_DIR']

pca = pd.read_csv(f"{pca_dir}/pca.eigenvec", sep=r'\s+')
pca.rename(columns={'#FID': 'FID'}, inplace=True)
pc_cols = [c for c in pca.columns if c.startswith('PC')]

cov = pd.read_csv(f"{pipeline_dir}/01_phenotype/covariates_noPCs.txt", sep=' ')
cov_pca = cov.merge(pca[['FID','IID'] + pc_cols], on=['FID','IID'], how='left')

use_pcs = pc_cols[:10]   # top 10 PCs as covariates
missing_pcs = cov_pca[use_pcs].isna().sum().sum()
print(f"  PCs merged. Missing PC entries for {missing_pcs} rows "
      f"(samples not in genotype data).")

cov_pca.to_csv(f"{pipeline_dir}/01_phenotype/covariates_withPCs.txt",
               sep=' ', index=False, na_rep='NA')
print(f"  Written: covariates_withPCs.txt  (n={len(cov_pca)})")

# Report variance explained
eig = pd.read_csv(f"{pca_dir}/pca.eigenval", header=None, names=['eigenval'])
eig['var_pct'] = 100 * eig['eigenval'] / eig['eigenval'].sum()
eig['cumvar']  = eig['var_pct'].cumsum()
print(f"\n  Variance explained by top 10 PCs: {eig['var_pct'][:10].sum():.1f}%")
for i, row in eig.head(10).iterrows():
    print(f"    PC{i+1}: {row['var_pct']:.2f}%  (cumulative: {row['cumvar']:.1f}%)")
PYMERGE

checkpoint "Covariates with PCs" "${PIPELINE_DIR}/01_phenotype/covariates_withPCs.txt"
echo "export COVAR_FILE=${PIPELINE_DIR}/01_phenotype/covariates_withPCs.txt" \
    >> ${PIPELINE_DIR}/pipeline_config.sh

# ─── 3.4: GENERATE PCA PLOTS ─────────────────────────────────────────────────
echo ""
echo "--- 3.4: Generating PCA plots ---"

${RSCRIPT} - <<'RPLOT'
library(ggplot2); library(dplyr); library(data.table)

pipeline_dir <- Sys.getenv("PIPELINE_DIR")
pca_dir      <- Sys.getenv("PCA_DIR")
fig_dir      <- Sys.getenv("FIG_DIR")

pca <- fread(paste0(pca_dir, "/pca.eigenvec"))
setnames(pca, "#FID", "FID")

pheno <- read.csv(paste0(pipeline_dir, "/01_phenotype/pheno_cox.txt"), sep="\t")
pca_ann <- merge(pca, pheno[, c("IID","CDX_clean","region","sex_bin")],
                 by="IID", all.x=TRUE)

eig <- read.table(paste0(pca_dir, "/pca.eigenval"), header=FALSE)
var_exp <- round(100 * eig$V1 / sum(eig$V1), 2)

region_colors <- c("Caribbean"="#2166ac","LatAm"="#d6604d","USA"="#4dac26","Missing"="grey60")
dx_colors     <- c("NCI"="#1a9850","MCI"="#fee08b","AD"="#d73027")

# Plot 1: PC1 vs PC2 colored by region
p1 <- ggplot(pca_ann, aes(PC1, PC2, color=region)) +
    geom_point(alpha=0.5, size=0.8) +
    scale_color_manual(values=region_colors, name="Region of origin") +
    labs(x=paste0("PC1 (", var_exp[1], "% variance)"),
         y=paste0("PC2 (", var_exp[2], "% variance)"),
         title="Population Structure — PC1 vs PC2",
         subtitle="Colored by region of birth") +
    theme_bw(base_size=11)

# Plot 2: PC1 vs PC2 colored by diagnosis
p2 <- ggplot(pca_ann %>% filter(!is.na(CDX_clean)),
             aes(PC1, PC2, color=CDX_clean)) +
    geom_point(alpha=0.5, size=0.8) +
    scale_color_manual(values=dx_colors, name="Diagnosis") +
    labs(x=paste0("PC1 (", var_exp[1], "% variance)"),
         y=paste0("PC2 (", var_exp[2], "% variance)"),
         title="Population Structure — PC1 vs PC2",
         subtitle="Colored by diagnosis") +
    theme_bw(base_size=11)

# Plot 3: Scree plot
scree_df <- data.frame(PC=1:10, var=var_exp[1:10],
                        cumvar=cumsum(var_exp[1:10]))
p3 <- ggplot(scree_df, aes(x=factor(PC), y=var)) +
    geom_col(fill="#4393c3", alpha=0.8) +
    geom_line(aes(group=1, y=cumvar), color="red", linewidth=0.8) +
    geom_point(aes(y=cumvar), color="red", size=2) +
    labs(x="Principal Component", y="Variance Explained (%)",
         title="Scree Plot — Top 10 PCs") +
    theme_bw(base_size=11)

pdf(paste0(fig_dir, "/fig_pca.pdf"), width=12, height=5)
gridExtra::grid.arrange(p1, p2, p3, ncol=3)
dev.off()
cat("  Written: fig_pca.pdf\n")
RPLOT

checkpoint "PCA figure" "${FIG_DIR}/fig_pca.pdf"

echo ""
echo "============================================================"
echo " STEP 3 COMPLETE ✅  [$(date)]"
echo " Next: Run step04a_gwas_regenie.sh and step04b_cox_gwas.R"
echo "============================================================"
