#!/usr/bin/env Rscript
# =============================================================================
# STEP 6: Paper-Ready Figures & Tables
# AOO GWAS Pipeline
# =============================================================================
# GENERATES:
#   Figure 1:  Manhattan plot (REGENIE linear AOO)
#   Figure 2:  QQ plot with λ GC annotation
#   Figure 3:  Miami plot (REGENIE vs Cox, mirrored) — if Cox results available
#   Figure 4:  AOO distribution by diagnosis, sex, APOE4 status, and region
#   Supp:      PCA scatter + scree (from step03)
#   Supp:      Kaplan-Meier curves (from step04b)
#   Supp:      Cox forest plot (from step04b)
#   Table 1:   Demographics (formatted .docx)
#   Table S1:  All GWS hits (.docx)
#   Table S2:  Top suggestive hits (.docx)
#
# All figures saved as 300 DPI TIFF + PDF (journal submission ready)
# =============================================================================

suppressPackageStartupMessages({
    library(data.table);   library(ggplot2);      library(dplyr)
    library(ggrepel);      library(RColorBrewer);  library(gridExtra)
    library(CMplot);       library(grid);          library(scales)
    library(tidyr);        library(stringr);       library(flextable)
    library(officer)
})

# ─── CONFIG ──────────────────────────────────────────────────────────────────
PIPELINE_DIR <- Sys.getenv("PIPELINE_DIR", ".")
FIG_DIR      <- file.path(PIPELINE_DIR, "09_figures")
ANNOT_DIR    <- file.path(PIPELINE_DIR, "05_annotation")
dir.create(FIG_DIR, showWarnings=FALSE)

log_file <- file.path(PIPELINE_DIR, "logs", "step06_figures.log")
con <- file(log_file, "a"); sink(con, append=TRUE, type="output")

cat("============================================================\n")
cat(sprintf(" STEP 6: PAPER FIGURES  [%s]\n", Sys.time()))
cat("============================================================\n")

save_fig <- function(plot_obj=NULL, filename, width, height, dpi=300) {
    path_pdf  <- file.path(FIG_DIR, paste0(filename, ".pdf"))
    path_tiff <- file.path(FIG_DIR, paste0(filename, ".tiff"))
    if (!is.null(plot_obj)) {
        ggsave(path_pdf,  plot=plot_obj, width=width, height=height)
        ggsave(path_tiff, plot=plot_obj, width=width, height=height,
               dpi=dpi, compression="lzw")
    } else {
        dev.copy(pdf,  path_pdf,  width=width, height=height);          dev.off()
        dev.copy(tiff, path_tiff, width=width*dpi, height=height*dpi,
                 units="px", res=dpi, compression="lzw");               dev.off()
    }
    cat(sprintf("  Saved: %s  (.pdf + .tiff)\n", filename))
}

# ─── LOAD GWAS RESULTS ────────────────────────────────────────────────────────
cat("\n--- Loading GWAS summary statistics ---\n")

gwas <- fread(file.path(PIPELINE_DIR, "04_gwas", "regenie",
                         "regenie_aoo_clean.txt.gz"))
gwas <- gwas %>%
    rename_with(~gsub("^#","",.), everything()) %>%
    rename(CHR=CHROM, POS=GENPOS, SNP=ID, P_linear=P) %>%
    mutate(P_linear = 10^(-LOG10P))

# Load Cox results if available
cox_file <- file.path(PIPELINE_DIR, "04_gwas", "cox",
                       "plink2_cox_aoo.EVENT.glm.cox")
has_cox <- file.exists(cox_file)
if (has_cox) {
    cox <- fread(cox_file)
    cox <- cox %>% rename(CHR=`#CHROM`, P_cox=P) %>%
        select(CHR, POS, ID, P_cox)
    gwas <- left_join(gwas, cox, by=c("CHR","POS","SNP"="ID"))
    cat("  Cox results merged\n")
} else {
    cat("  Cox results not yet available — Miami plot will be skipped\n")
    gwas$P_cox <- NA
}

cat(sprintf("  Variants: %d\n", nrow(gwas)))

# ─── FIGURE 1: MANHATTAN PLOT ─────────────────────────────────────────────────
cat("\n--- Figure 1: Manhattan plot ---\n")

gwas_m <- gwas %>%
    filter(!is.na(P_linear), P_linear > 0, CHR %in% 1:22) %>%
    arrange(CHR, POS) %>%
    group_by(CHR) %>%
    mutate(chr_len = max(POS)) %>%
    ungroup()

chr_offset <- gwas_m %>%
    group_by(CHR) %>%
    summarise(chr_len = max(POS), .groups='drop') %>%
    mutate(offset = lag(cumsum(as.numeric(chr_len)), default=0))

gwas_m <- gwas_m %>%
    left_join(chr_offset %>% select(CHR, offset), by="CHR") %>%
    mutate(pos_cum = POS + offset)

axis_df <- gwas_m %>%
    group_by(CHR) %>%
    summarise(center = mean(pos_cum), .groups='drop')

chr_colors  <- rep(c("#2166ac","#b2182b"), 11)
gws_thresh  <- -log10(5e-8)
sugg_thresh <- -log10(5e-6)

top_hits <- gwas_m %>%
    filter(-log10(P_linear) >= gws_thresh) %>%
    group_by(CHR) %>%
    slice_min(P_linear, n=1) %>%
    ungroup()

p_manhattan <- ggplot(gwas_m, aes(x=pos_cum, y=-log10(P_linear),
                                   color=factor(CHR))) +
    geom_point(size=0.4, alpha=0.7, show.legend=FALSE) +
    geom_hline(yintercept=gws_thresh,  linetype="dashed",
               color="red",  linewidth=0.6) +
    geom_hline(yintercept=sugg_thresh, linetype="dashed",
               color="blue", linewidth=0.4, alpha=0.6) +
    geom_point(data=gwas_m[gwas_m$P_linear < 5e-8, ],
               aes(x=pos_cum, y=-log10(P_linear)),
               color="red3", size=1.5, shape=18) +
    geom_label_repel(data=top_hits,
                     aes(x=pos_cum, y=-log10(P_linear), label=SNP),
                     size=2.5, color="black", fill="white",
                     box.padding=0.5, max.overlaps=20,
                     segment.color="grey50", show.legend=FALSE) +
    scale_color_manual(values=chr_colors) +
    scale_x_continuous(label=axis_df$CHR, breaks=axis_df$center) +
    scale_y_continuous(expand=c(0,0),
                       limits=c(0, max(-log10(gwas_m$P_linear), na.rm=TRUE)+1)) +
    labs(x="Chromosome", y=expression(-log[10](italic(p))),
         title="GWAS of Age of Onset — Admixed Population Cohort",
         subtitle=sprintf("N = %d cases | Linear regression (REGENIE) | covariates: sex, age, APOE-ε4, PC1-10, site",
                          max(gwas$N, na.rm=TRUE))) +
    annotate("text", x=max(gwas_m$pos_cum)*0.98, y=gws_thresh+0.3,
             label="p = 5×10⁻⁸", color="red", size=2.8, hjust=1) +
    theme_bw(base_size=11) +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.minor  =element_blank(),
          axis.text.x=element_text(size=8),
          plot.title=element_text(face="bold", size=12))

save_fig(p_manhattan, "fig1_manhattan", width=14, height=5)

# ─── FIGURE 2: QQ PLOT ────────────────────────────────────────────────────────
cat("\n--- Figure 2: QQ plot ---\n")

lambda_gc <- median(gwas$CHISQ, na.rm=TRUE) / qchisq(0.5, df=1)

qq_data <- gwas %>%
    filter(!is.na(P_linear), P_linear > 0) %>%
    arrange(P_linear) %>%
    mutate(
        obs     = -log10(P_linear),
        exp     = -log10(ppoints(n())),
        CI_low  = -log10(qbeta(0.975, seq_len(n()), n():1)),
        CI_high = -log10(qbeta(0.025, seq_len(n()), n():1))
    )

p_qq <- ggplot(qq_data, aes(x=exp, y=obs)) +
    geom_ribbon(aes(ymin=CI_low, ymax=CI_high), fill="lightblue", alpha=0.4) +
    geom_abline(intercept=0, slope=1, color="red", linewidth=0.7) +
    geom_point(size=0.4, alpha=0.6, color="#2166ac") +
    geom_point(data=qq_data %>% filter(obs > gws_thresh),
               color="red3", size=1.5, shape=18) +
    annotate("text", x=1, y=max(qq_data$obs)*0.95,
             label=sprintf("λ GC = %.4f", lambda_gc),
             hjust=0, size=4, fontface="bold") +
    labs(x=expression("Expected " * -log[10](italic(p))),
         y=expression("Observed " * -log[10](italic(p))),
         title="QQ Plot — AOO GWAS (REGENIE linear)") +
    theme_bw(base_size=12)

save_fig(p_qq, "fig2_qqplot", width=6, height=6)

# ─── FIGURE 3: MIAMI PLOT (if Cox results available) ─────────────────────────
if (has_cox) {
    cat("\n--- Figure 3: Miami plot (Linear vs Cox) ---\n")

    miami_data <- gwas %>%
        filter(!is.na(P_linear), !is.na(P_cox)) %>%
        left_join(chr_offset %>% select(CHR, offset), by="CHR") %>%
        mutate(pos_cum = POS + offset)

    p_miami <- ggplot(miami_data) +
        geom_point(aes(x=pos_cum, y=-log10(P_linear), color=factor(CHR)),
                   size=0.4, alpha=0.7, show.legend=FALSE) +
        geom_point(aes(x=pos_cum, y=log10(P_cox), color=factor(CHR)),
                   size=0.4, alpha=0.7, show.legend=FALSE) +
        geom_hline(yintercept=c(gws_thresh, -gws_thresh),
                   linetype="dashed", color="red", linewidth=0.6) +
        scale_color_manual(values=chr_colors) +
        scale_x_continuous(label=axis_df$CHR, breaks=axis_df$center) +
        labs(x="Chromosome", y=expression(-log[10](italic(p))),
             title="Miami Plot — Linear AOO (top) vs Cox PH (bottom)") +
        annotate("text", x=max(miami_data$pos_cum)*0.02, y=3,
                 label="Linear regression", hjust=0, size=3) +
        annotate("text", x=max(miami_data$pos_cum)*0.02, y=-3,
                 label="Cox PH model", hjust=0, size=3) +
        theme_bw(base_size=11) +
        theme(panel.grid.major.x=element_blank())

    save_fig(p_miami, "fig3_miami", width=14, height=6)
}

# ─── FIGURE 4: AOO DISTRIBUTION PLOTS ────────────────────────────────────────
cat("\n--- Figure 4: AOO distribution plots ---\n")

pheno <- fread(file.path(PIPELINE_DIR, "01_phenotype", "pheno_cox.txt"))
pheno$CDX_clean   <- factor(pheno$CDX_clean, levels=c("NCI","MCI","AD"))
pheno$APOE4_status <- ifelse(pheno$APOE_e4_count > 0, "APOE4 Carrier", "Non-carrier")

dx_colors <- c("NCI"="#1a9850","MCI"="#fee08b","AD"="#d73027")

# 4a: AOO density by diagnosis
p4a <- pheno %>%
    filter(!is.na(cox_time), CDX_clean %in% c("MCI","AD")) %>%
    ggplot(aes(x=cox_time, fill=CDX_clean)) +
    geom_density(alpha=0.6, linewidth=0.5) +
    scale_fill_manual(values=dx_colors[c("MCI","AD")], name="Diagnosis") +
    labs(x="Age of Onset (years)", y="Density",
         title="Age of Onset Distribution") +
    theme_bw(base_size=11)

# 4b: AOO boxplot by APOE4 × diagnosis
p4b <- pheno %>%
    filter(!is.na(cox_time), CDX_clean %in% c("MCI","AD"),
           !is.na(APOE4_status)) %>%
    ggplot(aes(x=CDX_clean, y=cox_time, fill=APOE4_status)) +
    geom_boxplot(outlier.size=0.5, alpha=0.8) +
    scale_fill_manual(values=c("APOE4 Carrier"="#ef8a62","Non-carrier"="#67a9cf"),
                      name="APOE4 Status") +
    labs(x="Diagnosis", y="Age of Onset (years)",
         title="AOO by Diagnosis and APOE4 Status") +
    theme_bw(base_size=11)

# 4c: AOO violin by sex
p4c <- pheno %>%
    filter(!is.na(cox_time), CDX_clean %in% c("MCI","AD")) %>%
    mutate(sex_label=ifelse(sex_bin==0,"Female","Male")) %>%
    ggplot(aes(x=sex_label, y=cox_time, fill=CDX_clean)) +
    geom_violin(alpha=0.7) +
    geom_boxplot(width=0.1, fill="white", outlier.size=0.3) +
    scale_fill_manual(values=dx_colors[c("MCI","AD")], name="Diagnosis") +
    labs(x="Sex", y="Age of Onset (years)",
         title="AOO by Sex and Diagnosis") +
    theme_bw(base_size=11)

# 4d: AOO by ancestry region
p4d <- pheno %>%
    filter(!is.na(cox_time), CDX_clean %in% c("MCI","AD"),
           region != "Missing") %>%
    ggplot(aes(x=region, y=cox_time, fill=region)) +
    geom_boxplot(alpha=0.8, show.legend=FALSE) +
    scale_fill_manual(values=c("Caribbean"="#2166ac","LatAm"="#d6604d","USA"="#4dac26")) +
    labs(x="Region of Origin", y="Age of Onset (years)",
         title="AOO by Region of Origin") +
    theme_bw(base_size=11)

p4_combined <- grid.arrange(p4a, p4b, p4c, p4d, ncol=2)
save_fig(filename="fig4_aoo_distributions", width=12, height=10)

# ─── TABLE 1: DEMOGRAPHICS ────────────────────────────────────────────────────
cat("\n--- Generating formatted Table 1 ---\n")

table1_raw <- read.csv(file.path(PIPELINE_DIR, "01_phenotype",
                                  "table1_demographics.csv"))

ft <- flextable(table1_raw) %>%
    bold(part="header") %>%
    bg(bg="#f0f4f8", part="header") %>%
    border_outer(part="all") %>%
    border_inner_h(part="body") %>%
    autofit() %>%
    set_caption("Table 1. Baseline demographic and clinical characteristics")

doc <- read_docx() %>% body_add_flextable(ft)
print(doc, target=file.path(FIG_DIR, "table1_demographics.docx"))
cat("  Saved: table1_demographics.docx\n")

# ─── SUPPLEMENTARY TABLE S1 ──────────────────────────────────────────────────
cat("\n--- Supplementary Table S1: GWS hits ---\n")

if (file.exists(file.path(ANNOT_DIR, "suppTable_GWS_hits.csv"))) {
    s1 <- read.csv(file.path(ANNOT_DIR, "suppTable_GWS_hits.csv"))
    if (nrow(s1) > 0) {
        ft_s1 <- flextable(s1) %>%
            bold(part="header") %>%
            bg(bg="#f0f4f8", part="header") %>%
            border_outer(part="all") %>%
            autofit() %>%
            set_caption("Table S1. Genome-wide significant associations with age of onset")
        doc_s1 <- read_docx() %>% body_add_flextable(ft_s1)
        print(doc_s1, target=file.path(FIG_DIR, "tableS1_GWS_hits.docx"))
        cat("  Saved: tableS1_GWS_hits.docx\n")
    } else {
        cat("  No GWS hits to table (threshold p<5e-8 not reached)\n")
        cat("  Consider reporting suggestive hits in Table S1\n")
    }
}

cat("\n")
cat("============================================================\n")
cat(sprintf(" STEP 6 COMPLETE ✅  [%s]\n", Sys.time()))
cat(sprintf(" All figures in: %s\n", FIG_DIR))
cat(" Next: Run step07_magma.sh for gene/pathway analysis\n")
cat("============================================================\n")

sink(); sink(type="message")
