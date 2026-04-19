#!/usr/bin/env python3
# =============================================================================
# STEP 5: Post-GWAS Filtering, LD Clumping & Annotation
# AOO GWAS Pipeline
# =============================================================================
# PURPOSE:
#   1. Load REGENIE and Cox GWAS results
#   2. Filter on INFO score ≥ 0.8 and MAF ≥ 1%
#   3. LD-clump to identify independent signals
#   4. Annotate top hits with nearest known AD gene
#   5. Generate supplementary tables of significant hits
#   6. Calculate genomic inflation (λ GC) for both analyses
# =============================================================================

import pandas as pd
import numpy as np
import os, sys, subprocess, logging
from scipy import stats

PIPELINE_DIR = os.environ.get("PIPELINE_DIR", ".")
GWAS_DIR     = os.path.join(PIPELINE_DIR, "04_gwas")
ANNOT_DIR    = os.path.join(PIPELINE_DIR, "05_annotation")
FIG_DIR      = os.path.join(PIPELINE_DIR, "09_figures")
os.makedirs(ANNOT_DIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(message)s",
    handlers=[
        logging.FileHandler(f"{PIPELINE_DIR}/logs/step05_annotation.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
log = logging.getLogger()

def checkpoint(name, condition):
    if condition:
        log.info(f"  ✅ CHECKPOINT: {name}")
    else:
        log.error(f"  ❌ CHECKPOINT FAILED: {name}")
        sys.exit(1)

def calc_lambda_gc(chisq_array):
    """Calculate genomic inflation factor from chi-squared statistics."""
    valid = chisq_array[~np.isnan(chisq_array) & (chisq_array > 0)]
    return np.median(valid) / stats.chi2.ppf(0.5, df=1)

# ─── 5.1: LOAD REGENIE RESULTS ────────────────────────────────────────────────
log.info("=" * 60)
log.info(" STEP 5: POST-GWAS FILTERING & ANNOTATION")
log.info("=" * 60)

regenie_file = f"{GWAS_DIR}/regenie/regenie_aoo_clean.txt.gz"
log.info(f"\n--- Loading REGENIE results ---")
log.info(f"  File: {regenie_file}")

df_reg = pd.read_csv(regenie_file, sep='\t', compression='gzip')
log.info(f"  Variants loaded: {len(df_reg):,}")

# Standardize column names
df_reg = df_reg.rename(columns={
    'CHROM': 'CHR', 'GENPOS': 'POS', 'ID': 'SNP',
    'ALLELE0': 'REF', 'ALLELE1': 'ALT',
    'A1FREQ': 'MAF'
})
df_reg['P'] = 10 ** (-df_reg['LOG10P'])
df_reg['source'] = 'REGENIE_linear'

# ─── 5.2: QUALITY FILTERS ─────────────────────────────────────────────────────
log.info("\n--- Applying post-GWAS quality filters ---")

n_before = len(df_reg)
df_reg_qc = df_reg[
    (df_reg['MAF'] >= 0.01) &        # MAF ≥ 1%
    (df_reg['INFO'] >= 0.8) &        # Imputation quality ≥ 0.8
    (df_reg['N'] >= 100) &           # Minimum sample size
    df_reg['P'].notna() &
    (df_reg['P'] > 0)
].copy()

log.info(f"  Before QC filters: {n_before:,}")
log.info(f"  After QC filters:  {len(df_reg_qc):,}")
log.info(f"  Removed: {n_before - len(df_reg_qc):,} variants")

checkpoint("QC'd variants remain", len(df_reg_qc) > 100000)

# ─── 5.3: GENOMIC INFLATION ───────────────────────────────────────────────────
log.info("\n--- Calculating genomic inflation ---")

lambda_gc = calc_lambda_gc(df_reg_qc['CHISQ'].values)
lambda_1000 = 1 + (lambda_gc - 1) * (1/df_reg_qc['N'].median() - 1/1000) / \
              (1/df_reg_qc['N'].median())

log.info(f"  λ GC (REGENIE): {lambda_gc:.4f}")
log.info(f"  λ GC scaled to N=1000: {lambda_1000:.4f}")

if lambda_gc > 1.10:
    log.warning("  ⚠️  Inflation detected — consider additional PCs")
elif lambda_gc < 0.90:
    log.warning("  ⚠️  Deflation — check phenotype/covariate specification")
else:
    log.info("  ✅ λ GC within acceptable range")

# ─── 5.4: IDENTIFY SIGNIFICANT HITS ──────────────────────────────────────────
log.info("\n--- Identifying significant associations ---")

P_GWS  = 5e-8   # genome-wide significance
P_SUGG = 5e-6   # suggestive significance

hits_gws  = df_reg_qc[df_reg_qc['P'] < P_GWS].copy()
hits_sugg = df_reg_qc[(df_reg_qc['P'] >= P_GWS) & (df_reg_qc['P'] < P_SUGG)].copy()

log.info(f"  GWS hits (p<5e-8):    {len(hits_gws):,} variants")
log.info(f"  Suggestive (p<5e-6):  {len(hits_sugg):,} variants")

if len(hits_gws) > 0:
    log.info("\n  Top GWS hits:")
    top = df_reg_qc.nsmallest(20, 'P')[['CHR','POS','SNP','REF','ALT','MAF','BETA','SE','P','N']]
    log.info(f"\n{top.to_string(index=False)}")

# ─── 5.5: LD CLUMPING ────────────────────────────────────────────────────────
log.info("\n--- LD Clumping to identify independent signals ---")

pval_file = f"{ANNOT_DIR}/pvals_for_clump.txt"
df_reg_qc[['SNP','P']].to_csv(pval_file, sep='\t', index=False)

geno_qc = os.path.join(PIPELINE_DIR, "02_genotype_qc", "step2b_sampleQC")

# Note: Using PLINK 1.9 for clumping (PLINK2 uses different clumping syntax)
clump_cmd = (
    f"plink --bfile {geno_qc} "
    f"--clump {pval_file} "
    f"--clump-p1 5e-8 --clump-p2 5e-6 "
    f"--clump-r2 0.1 --clump-kb 500 "
    f"--clump-field P "
    f"--out {ANNOT_DIR}/clumped_hits "
    f"--allow-no-sex"
)

log.info(f"  Running: {clump_cmd}")
result = subprocess.run(clump_cmd, shell=True, capture_output=True, text=True)
log.info(result.stdout)
if result.returncode != 0:
    log.warning(f"  Clumping returned non-zero: {result.stderr}")
    log.warning("  Proceeding with unclumped hits for annotation")

# ─── 5.6: ANNOTATE TOP HITS ──────────────────────────────────────────────────
log.info("\n--- Annotating top hits ---")

top_hits = df_reg_qc.nsmallest(100, 'P').copy()

# Known AD loci for cross-referencing (GRCh38 coordinates)
# Update or extend this dictionary as new loci are discovered
known_ad_loci = {
    'APOE':   (19, 44905796, 45909393),
    'BIN1':   (2,  127770747, 127896227),
    'CLU':    (8,  27454434, 27524525),
    'CR1':    (1,  207492607, 207920776),
    'PICALM': (11, 85668236, 85962958),
    'MS4A4A': (11, 59859836, 59954175),
    'CD33':   (19, 51224167, 51232259),
    'ABCA7':  (19, 1040102, 1065572),
    'EPHA1':  (7,  143099574, 143188481),
    'CD2AP':  (6,  47370108, 47592655),
}

def annotate_known_locus(row, known_loci, window_kb=500):
    """Return nearest known AD locus within window, or 'Novel'."""
    for gene, (chrom, start, end) in known_loci.items():
        window = window_kb * 1000
        if (row['CHR'] == chrom and
                row['POS'] >= (start - window) and
                row['POS'] <= (end + window)):
            return gene
    return 'Novel'

top_hits['nearest_known'] = top_hits.apply(
    lambda r: annotate_known_locus(r, known_ad_loci), axis=1)
top_hits['novel_flag'] = top_hits['nearest_known'] == 'Novel'

log.info(f"  Novel loci (not near known AD genes): "
         f"{top_hits['novel_flag'].sum()}/{len(top_hits)}")

top_hits.to_csv(f"{ANNOT_DIR}/top100_hits_annotated.txt", sep='\t', index=False)
log.info("  Saved: top100_hits_annotated.txt")

# ─── 5.7: SUPPLEMENTARY TABLE ─────────────────────────────────────────────────
log.info("\n--- Generating supplementary tables ---")

supp_cols = ['CHR','POS','SNP','REF','ALT','MAF','BETA','SE','P','N','nearest_known']
supp_gws = top_hits[top_hits['P'] < P_GWS][supp_cols].copy()
supp_gws.columns = ['Chr','Position','rsID','Ref','Alt','MAF',
                     'Beta','SE','P-value','N','Nearest_known_locus']
supp_gws['P-value'] = supp_gws['P-value'].apply(lambda x: f"{x:.2e}")
supp_gws['Beta'] = supp_gws['Beta'].round(4)
supp_gws['SE']   = supp_gws['SE'].round(4)
supp_gws['MAF']  = supp_gws['MAF'].round(4)

supp_gws.to_csv(f"{ANNOT_DIR}/suppTable_GWS_hits.csv", index=False)
log.info(f"  Saved: suppTable_GWS_hits.csv  ({len(supp_gws)} rows)")

# Summary stats
summary = {
    'Total variants tested':    f"{len(df_reg_qc):,}",
    'Lambda GC':                f"{lambda_gc:.4f}",
    'GWS hits (p<5e-8)':       str(len(hits_gws)),
    'Suggestive hits (p<5e-6)': str(len(hits_sugg)),
    'Novel loci':               str(top_hits['novel_flag'].sum()),
    'Min P-value':              f"{df_reg_qc['P'].min():.2e}",
}
pd.DataFrame.from_dict(summary, orient='index',
                        columns=['Value']).to_csv(f"{ANNOT_DIR}/gwas_summary_stats.csv")

log.info("\n  GWAS Summary:")
for k, v in summary.items():
    log.info(f"    {k}: {v}")

log.info("\n" + "=" * 60)
log.info(" STEP 5 COMPLETE ✅")
log.info(" Next: Run step06_visualization.R for paper-ready figures")
log.info("=" * 60)
