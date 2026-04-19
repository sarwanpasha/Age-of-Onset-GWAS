#!/bin/bash
# =============================================================================
# MASTER RUN SCRIPT — AOO GWAS Pipeline
# =============================================================================
# USAGE:
#   Option A — Run interactively (not recommended for Steps 4+):
#     bash run_pipeline.sh
#
#   Option B — Submit to SLURM (recommended):
#     sbatch run_pipeline.slurm
#
#   Option C — Run individual steps:
#     bash 00_scripts/step01_phenotype_prep.py  (etc.)
#
# ESTIMATED RUNTIME ON HPC (16-core node):
#   Step 0–1:  ~5 min
#   Step 2:    20–60 min (genotype QC)
#   Step 3:    15–30 min (PCA)
#   Step 4a:   2–6 hrs   (REGENIE — parallelise by chr for speed)
#   Step 4b:   1–3 hrs   (Cox GWAS)
#   Step 5:    ~10 min   (post-GWAS)
#   Step 6:    ~10 min   (figures)
#   Step 7:    ~30 min   (MAGMA)
#   Total:     ~8–12 hrs on 16-core HPC node
# =============================================================================

set -euo pipefail

# ─── EDIT THESE BEFORE RUNNING ───────────────────────────────────────────────
export PHENO_CSV="/path/to/phenotype.csv"
export GENO_PREFIX="/path/to/genotype/all_chr"
export PIPELINE_DIR="/path/to/output/directory"
export THREADS=16
export MEMORY=32000   # MB
# ─────────────────────────────────────────────────────────────────────────────

SCRIPTS="${PIPELINE_DIR}/00_scripts"
mkdir -p ${PIPELINE_DIR}/logs

echo "========================================================"
echo " AOO GWAS PIPELINE — STARTING  [$(date)]"
echo "========================================================"

run_step() {
    local step_num="$1"
    local step_name="$2"
    local cmd="$3"
    echo ""
    echo "────────────────────────────────────────────────────"
    echo "  RUNNING STEP ${step_num}: ${step_name}"
    echo "  Started: $(date)"
    echo "────────────────────────────────────────────────────"
    eval "$cmd"
    echo "  ✅ STEP ${step_num} COMPLETE  [$(date)]"
}

run_step 0  "Setup & Software Check"      "bash ${SCRIPTS}/step00_setup.sh"
run_step 1  "Phenotype Preparation"       "python3 ${SCRIPTS}/step01_phenotype_prep.py"
run_step 2  "Genotype QC"                 "bash ${SCRIPTS}/step02_genotype_qc.sh"
run_step 3  "PCA"                         "bash ${SCRIPTS}/step03_pca.sh"
run_step 4a "GWAS Linear (REGENIE)"       "bash ${SCRIPTS}/step04a_gwas_regenie.sh"
run_step 4b "GWAS Cox (PLINK2 + R)"       "Rscript ${SCRIPTS}/step04b_cox_gwas.R && bash ${PIPELINE_DIR}/04_gwas/cox/run_plink2_cox.sh"
run_step 5  "Post-GWAS Annotation"        "python3 ${SCRIPTS}/step05_postGWAS.py"
run_step 6  "Figures & Tables"            "Rscript ${SCRIPTS}/step06_visualization.R"
run_step 7  "MAGMA Gene/Pathway"          "bash ${SCRIPTS}/step07_magma.sh"

echo ""
echo "========================================================"
echo " PIPELINE COMPLETE ✅  [$(date)]"
echo " Results in: ${PIPELINE_DIR}/09_figures/"
echo "========================================================"

# ─── SLURM SUBMISSION SCRIPT ──────────────────────────────────────────────────
cat > ${PIPELINE_DIR}/run_pipeline.slurm << 'SLURM'
#!/bin/bash
#SBATCH --job-name=aoo_gwas
#SBATCH --output=logs/pipeline_%j.log
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your@email.com   # <-- change this

module load plink2/2.0
module load regenie/3.2
module load R/4.3.0
module load python/3.10

bash run_pipeline.sh
SLURM
echo "  SLURM script saved: run_pipeline.slurm"
echo "  Submit with: sbatch run_pipeline.slurm"
