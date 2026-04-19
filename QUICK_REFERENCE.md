# AOO GWAS Pipeline — Quick Reference

---

## BEFORE YOU START — 3 Things to Edit

In `run_pipeline.sh`, set:
```bash
PHENO_CSV="/path/to/phenotype.csv"
GENO_PREFIX="/path/to/genotype/all_chr"
PIPELINE_DIR="/path/to/output/directory"
```

---

## PIPELINE STEPS AT A GLANCE

| Step | Script | Language | Runtime | Key Output |
|------|--------|----------|---------|------------|
| 0 | `step00_setup.sh` | Bash | 2 min | Config file, dir structure |
| 1 | `step01_phenotype_prep.py` | Python | 1 min | Cleaned phenotype files, Table 1 |
| 2 | `step02_genotype_qc.sh` | Bash/PLINK2 | 30–60 min | QC'd .bed/.bim/.fam |
| 3 | `step03_pca.sh` | Bash/R | 20 min | PCA eigenvectors, PCA plot |
| 4a | `step04a_gwas_regenie.sh` | Bash/REGENIE | 2–6 hrs | GWAS summary stats |
| 4b | `step04b_cox_gwas.R` | R | 1–3 hrs | Cox GWAS, KM curves |
| 5 | `step05_postGWAS.py` | Python | 10 min | Annotated hits, λ GC |
| 6 | `step06_visualization.R` | R | 10 min | All paper figures |
| 7 | `step07_magma.sh` | Bash/MAGMA | 30 min | Gene/pathway results |

---

## KEY PARAMETERS

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| MAF filter | ≥ 1% | Standard for common variant GWAS |
| INFO score | ≥ 0.8 | Imputation quality threshold |
| HWE p-value | > 1e-6 | Remove poorly genotyped variants |
| PCs included | 10 | Controls for population stratification in admixed cohorts |
| GWS threshold | p < 5×10⁻⁸ | Standard genome-wide significance |
| Suggestive | p < 5×10⁻⁶ | Report for replication |
| LD clumping | r² < 0.1, 500 kb | Independent signal identification |

---

## PHENOTYPE CODING

| Variable | Coding | Notes |
|----------|--------|-------|
| AOO_zscore | Z-score normalized | REGENIE phenotype (cases only) |
| cox_time | AOO (cases) / age (controls) | Cox time-to-event |
| cox_event | 1 = AD/MCI, 0 = NCI | Cox event indicator |
| sex_bin | 0 = F, 1 = M | Standard genetic convention |
| APOE_e4_count | 0/1/2 | Number of ε4 alleles |

---

## OUTPUT FILES FOR PAPER

```
09_figures/
├── fig1_manhattan.pdf/.tiff      # Figure 1
├── fig2_qqplot.pdf/.tiff         # Figure 2
├── fig3_miami.pdf/.tiff          # Figure 3 (if Cox converges)
├── fig4_aoo_distributions.pdf    # Figure 4
├── fig_pca.pdf                   # Supplementary
├── fig_km_curves.pdf             # Supplementary
├── fig_cox_forest.pdf            # Supplementary
├── table1_demographics.docx      # Table 1
└── tableS1_GWS_hits.docx         # Table S1
```

---

## TROUBLESHOOTING

**IDs don't match between phenotype and genotype:**
→ Check that SAMPLE_ID in phenotype CSV matches IDs in .fam file exactly
→ Run: `cut -f2 all_chr.fam | head -5` and compare to SAMPLE_ID column

**REGENIE Step 1 runs out of memory:**
→ Increase `--memory` or use `--lowmem` flag (already included)

**λ GC > 1.10:**
→ Add more PCs (try 15–20), check for residual batch effects by site

**Few GWS hits:**
→ Report suggestive hits (p < 5e-6) — expected for moderate AOO case counts
→ Consider meta-analysis with external cohorts

**PLINK2 Cox not supported:**
→ Use R `survival` package for top hits validation
→ The Cox command in step04b is pre-written for your PLINK2 version
