# Age of Onset GWAS Pipeline — Admixed Population Cohort

A modular, HPC-ready pipeline for genome-wide association study (GWAS) of **Alzheimer's Disease Age of Onset (AOO)** in an admixed cohort. The pipeline combines linear and survival (Cox PH) regression approaches, with full quality control, visualization, and gene/pathway analysis steps.

---

## Overview

This pipeline performs a multi-modal GWAS of AOO — treating it as both a quantitative trait (linear regression via REGENIE) and a time-to-event outcome (Cox proportional hazards). It is designed for imputed genotype data from admixed populations and includes population stratification correction via PCA.

### Pipeline Architecture

```
Step 0  → Environment setup & software validation
Step 1  → Phenotype preparation & QC
Step 2  → Genotype QC (PLINK2)
Step 3  → Population stratification / PCA
Step 4a → GWAS: Linear regression on AOO (REGENIE)
Step 4b → GWAS: Survival analysis, AOO as time-to-event (Cox PH, R)
Step 5  → Post-GWAS filtering, LD clumping & annotation
Step 6  → Paper-ready figures & tables
Step 7  → Gene-based & pathway analysis (MAGMA)
```

---

## Features

- **Dual GWAS models**: linear regression (REGENIE) + Cox proportional hazards (PLINK2 + R `survival`)
- **Admixture-aware**: 10 PCs computed from LD-pruned variants, excluding high-LD regions (MHC, chr8 inversion)
- **REGENIE whole-genome regression**: handles cryptic relatedness and population structure in admixed samples
- **Full sample & variant QC**: call rate, heterozygosity, sex check, KING relatedness, MAF, HWE, INFO score
- **Kaplan-Meier survival curves**: stratified by APOE4 status, sex, and region of origin
- **Automated checkpoints**: every step validates outputs before proceeding; pipeline halts with informative messages on failure
- **HPC-ready**: SLURM submission scripts included for each long-running step
- **Paper-ready outputs**: Manhattan, QQ, Miami, locus zoom, PCA, KM plots; Table 1 and supplementary tables as `.docx`

---

## Requirements

### Software

| Tool | Version | Purpose |
|------|---------|---------|
| PLINK2 | ≥ 2.0 | Genotype QC, PCA, association testing |
| REGENIE | ≥ 3.2 | Whole-genome regression GWAS |
| R | ≥ 4.3 | Cox model, visualization, tables |
| Python | ≥ 3.10 | Phenotype prep, post-GWAS QC |
| MAGMA | ≥ 1.10 | Gene/pathway analysis |
| PLINK 1.9 | any | LD clumping |

### R Packages

```r
install.packages(c(
  "data.table", "ggplot2", "ggrepel", "dplyr", "tidyr",
  "survival", "survminer", "qqman", "CMplot", "coloc",
  "susieR", "SKAT", "RColorBrewer", "gridExtra",
  "tableone", "gtsummary", "flextable", "officer"
))
```

### Python Packages

```bash
pip install pandas numpy scipy matplotlib seaborn lifelines statsmodels
```

---

## Input Data

| File | Description |
|------|-------------|
| `phenotype.csv` | Sample-level phenotype and covariate file |
| `genotype.bed/.bim/.fam` | PLINK-format imputed genotype data (all chromosomes merged) |

### Required Phenotype Columns

| Column | Description |
|--------|-------------|
| `SAMPLE_ID` | Unique sample identifier matching genotype IDs |
| `sex` | Sex (`F` / `M`) |
| `age_at_subject` | Age at assessment |
| `AOO` | Age of onset (cases only; `.` for missing) |
| `CDX` | Diagnosis: `AD`, `MCI`, or `NCI` |
| `APOE` | APOE genotype code (e.g., `33`, `34`, `44`) |
| `site` | Recruitment site identifier |
| `country_of_birth` | Country/region of birth (used for ancestry grouping) |

---

## Quick Start

### 1. Clone and configure

```bash
git clone https://github.com/your-org/aoo-gwas-pipeline.git
cd aoo-gwas-pipeline
```

Edit the three required paths in `run_pipeline.sh`:

```bash
PHENO_CSV="/path/to/phenotype.csv"
GENO_PREFIX="/path/to/genotype/all_chr"
PIPELINE_DIR="/path/to/output/directory"
```

### 2. Run interactively (small test datasets)

```bash
bash run_pipeline.sh
```

### 3. Submit to SLURM (recommended for full datasets)

```bash
# Edit email address in the SLURM header first, then:
sbatch run_pipeline.slurm
```

### 4. Run individual steps

```bash
bash 00_scripts/step00_setup.sh
python3 00_scripts/step01_phenotype_prep.py
bash 00_scripts/step02_genotype_qc.sh
bash 00_scripts/step03_pca.sh
bash 00_scripts/step04a_gwas_regenie.sh
Rscript 00_scripts/step04b_cox_gwas.R
python3 00_scripts/step05_postGWAS.py
Rscript 00_scripts/step06_visualization.R
bash 00_scripts/step07_magma.sh
```

---

## Key Parameters

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| MAF filter | ≥ 1% | Standard for common variant GWAS |
| INFO score | ≥ 0.8 | Imputation quality threshold |
| HWE p-value | > 1×10⁻⁶ | Removes poorly genotyped variants |
| Heterozygosity | ±3 SD | Sample-level outlier removal |
| Relatedness | KING > 0.0884 | Flagged; > 0.354 removed |
| PCs included | 10 | Controls admixture stratification |
| GWS threshold | p < 5×10⁻⁸ | Standard genome-wide significance |
| Suggestive threshold | p < 5×10⁻⁶ | Reported for replication follow-up |
| LD clumping | r² < 0.1, 500 kb | Independent signal identification |

---

## Phenotype Coding

| Variable | Coding | Notes |
|----------|--------|-------|
| `AOO_zscore` | Z-score normalized AOO | REGENIE phenotype (cases only) |
| `cox_time` | AOO (cases) / age at assessment (controls) | Cox time variable |
| `cox_event` | 1 = AD/MCI, 0 = NCI | Cox event indicator |
| `sex_bin` | 0 = Female, 1 = Male | Standard genetic convention |
| `APOE_e4_count` | 0 / 1 / 2 | Number of ε4 alleles |

---

## Output Structure

```
Pipeline_Output/
├── 00_scripts/           ← All pipeline scripts
├── 01_phenotype/         ← Cleaned phenotype & covariate files
│   ├── pheno_aoo_regenie.txt
│   ├── pheno_cox.txt
│   ├── covariates_noPCs.txt
│   ├── covariates_withPCs.txt
│   └── table1_demographics.csv
├── 02_genotype_qc/       ← QC-passed PLINK binary files
├── 03_pca/               ← PCA eigenvectors & eigenvalues
├── 04_gwas/
│   ├── regenie/          ← REGENIE linear AOO results
│   └── cox/              ← Cox GWAS results
├── 05_annotation/        ← Annotated hits, clumped signals
├── 06_finemapping/       ← Fine-mapping results (SuSiE/COJO)
├── 07_magma/             ← Gene-level & pathway results
├── 08_ldsc/              ← Genetic correlation outputs
├── 09_figures/           ← All paper-ready figures & tables
└── logs/                 ← Per-step log files
```

### Key Output Files for Publication

```
09_figures/
├── fig1_manhattan.pdf/.tiff     ← Figure 1: Manhattan plot
├── fig2_qqplot.pdf/.tiff        ← Figure 2: QQ plot with λ GC
├── fig3_miami.pdf/.tiff         ← Figure 3: Miami plot (linear vs Cox)
├── fig4_aoo_distributions.pdf   ← Figure 4: AOO distributions
├── fig_pca.pdf                  ← Supplementary: PCA scatter + scree
├── fig_km_curves.pdf            ← Supplementary: Kaplan-Meier curves
├── fig_cox_forest.pdf           ← Supplementary: Cox forest plot
├── table1_demographics.docx     ← Table 1
└── tableS1_GWS_hits.docx        ← Supplementary Table S1
```

---

## MAGMA Reference Files

Download from the [MAGMA website](https://ctg.cncr.nl/software/magma):

```bash
# LD reference panel (use population-matched if available)
wget https://ctg.cncr.nl/software/MAGMA/aux_files/g1000_eur.zip
unzip g1000_eur.zip -d 07_magma/reference/

# Gene location file (GRCh37/hg19)
wget https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.gene.loc.zip
unzip NCBI37.3.gene.loc.zip -d 07_magma/reference/

# Gene sets (MSigDB)
# Download from: https://www.gsea-msigdb.org/gsea/msigdb/
```

> **Note for admixed populations:** The pipeline uses a European LD reference as an approximation. For best results, use a matched LD panel (e.g., 1000 Genomes AMR). This limitation should be stated in any manuscript.

---

## Troubleshooting

**Sample IDs don't match between phenotype and genotype files:**
```bash
# Check IDs in the .fam file
cut -f2 all_chr.fam | head -10
# Compare to SAMPLE_ID column in phenotype CSV
```

**REGENIE Step 1 runs out of memory:**
- Increase `--memory` in `step04a_gwas_regenie.sh`
- The `--lowmem` flag is already enabled by default

**Genomic inflation λ GC > 1.10:**
- Increase number of PCs (try 15–20)
- Check for residual batch effects by site
- Verify that site covariates are correctly specified

**Few or no genome-wide significant hits:**
- Expected for moderate sample sizes — report suggestive hits (p < 5×10⁻⁶)
- Consider meta-analysis with external cohorts

**PLINK2 Cox model not supported:**
- Use the R `survival` package fallback in `step04b_cox_gwas.R`
- The script pre-generates the Cox command for your PLINK2 version

---

## Runtime Estimates (16-core HPC node)

| Step | Runtime |
|------|---------|
| Step 0–1 | ~5 min |
| Step 2 (Genotype QC) | 20–60 min |
| Step 3 (PCA) | 15–30 min |
| Step 4a (REGENIE) | 2–6 hrs |
| Step 4b (Cox GWAS) | 1–3 hrs |
| Step 5–6 | ~20 min |
| Step 7 (MAGMA) | ~30 min |
| **Total** | **~8–12 hrs** |

---

## License

MIT License. See `LICENSE` for details.
