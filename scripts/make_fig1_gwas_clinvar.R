#!/usr/bin/env Rscript

defaults <- list(
  gwas_catalog = file.path("data", "raw", "gwas_catalog_v1.0.2-associations_e115_r2026-04-27_full.tsv"),
  clinvar = file.path("data", "raw", "variant_summary.txt.gz"),
  joined_clinvar_gnomad = file.path("data", "processed", "clinvar_gnomad_v4_joined_allele.tsv"),
  outdir = file.path("results", "figures"),
  prefix = "fig1_gwas_catalog_clinvar_gnomad",
  gwas_release = "GWAS Catalog e115 r2026-04-27",
  clinvar_label = "ClinVar May 2026",
  gnomad_label = "gnomAD v4",
  clinvar_af_column = "AF",
  clinvar_top_n = "10"
)

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/make_fig1_gwas_clinvar.R [options]\n\n",
    "Options:\n",
    "  --gwas-catalog PATH            GWAS Catalog associations TSV\n",
    "  --clinvar PATH                 ClinVar variant_summary.txt.gz\n",
    "  --joined-clinvar-gnomad PATH   Output from prep_gnomad_clinvar_join.py\n",
    "  --outdir DIR                   Output directory\n",
    "  --prefix NAME                  Output file prefix\n",
    "  --clinvar-af-column NAME       AF column to plot: AF or AF_GRPMAX\n",
    "  --clinvar-top-n N              Number of common ClinVar labels to retain\n",
    "  --help                         Show this message\n",
    sep = ""
  )
}

parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args || "-h" %in% args) {
    usage()
    quit(save = "no", status = 0)
  }
  opts <- defaults
  if (length(args) == 0) {
    return(opts)
  }
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    key <- gsub("-", "_", sub("^--", "", key))
    if (!key %in% names(opts)) {
      stop("Unknown option: --", gsub("_", "-", key), call. = FALSE)
    }
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      stop("Missing value for --", gsub("_", "-", key), call. = FALSE)
    }
    opts[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  opts
}

opts <- parse_args(defaults)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(viridis)
})

assert_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

assert_columns <- function(dat, cols, label) {
  missing <- setdiff(cols, names(dat))
  if (length(missing) > 0) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

numeric_column <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

assert_file(opts$gwas_catalog, "GWAS Catalog file")
assert_file(opts$clinvar, "ClinVar file")
assert_file(opts$joined_clinvar_gnomad, "ClinVar x gnomAD joined file")
dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

# Panel A: GWAS Catalog significant associations by minor allele frequency.
gwas_assoc <- fread(opts$gwas_catalog, quote = "")
assert_columns(
  gwas_assoc,
  c("CHR_ID", "CHR_POS", "RISK ALLELE FREQUENCY", "PVALUE_MLOG"),
  "GWAS Catalog"
)

gwas_sig <- gwas_assoc[, .(
  chrom = as.character(CHR_ID),
  pos = as.integer(numeric_column(CHR_POS)),
  risk_af = numeric_column(`RISK ALLELE FREQUENCY`),
  pvalue_mlog = numeric_column(PVALUE_MLOG)
)]

gwas_sig <- gwas_sig[
  !is.na(risk_af) &
    risk_af > 0 &
    risk_af < 1 &
    !is.na(pvalue_mlog) &
    pvalue_mlog > 8 &
    chrom != "" &
    !grepl(";", chrom)
]

gwas_sig[, maf := pmin(risk_af, 1 - risk_af)]
gwas_sig[, pval_range := fcase(
  pvalue_mlog > 60, "< 1e-60",
  pvalue_mlog > 30, "1e-60 to 1e-30",
  pvalue_mlog > 15, "1e-30 to 1e-15",
  default = "1e-15 to 1e-8"
)]
gwas_sig[, pval_range := factor(
  pval_range,
  levels = c("1e-15 to 1e-8", "1e-30 to 1e-15", "1e-60 to 1e-30", "< 1e-60"),
  ordered = TRUE
)]

axis_breaks <- c(0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5)
theme_fig <- theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

gwas_panel <- ggplot(gwas_sig, aes(x = maf, fill = pval_range)) +
  geom_histogram(position = "stack", color = "gray60", alpha = 0.85, linewidth = 0.2, bins = 30) +
  theme_fig +
  xlab("Minor Allele Frequency") +
  ylab("Count of associations with P < 1e-8") +
  ggtitle("GWAS Catalog allele frequency by association strength") +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 1e-6, base = 10),
    expand = c(0, 0),
    breaks = axis_breaks
  ) +
  scale_fill_viridis_d(option = "D", direction = -1, end = 0.9) +
  labs(fill = "P-value range")

# Panel B: ClinVar classifications by gnomAD allele frequency.
cv_gn <- fread(opts$joined_clinvar_gnomad)
assert_columns(
  cv_gn,
  c("CHROM", "POS", "REF", "ALT", opts$clinvar_af_column, "ClinicalSignificance"),
  "ClinVar x gnomAD joined file"
)

plot_af_col <- opts$clinvar_af_column
cv_gn[, plot_af_raw := numeric_column(get(plot_af_col))]
if (plot_af_col == "AF_GRPMAX" && "AF" %in% names(cv_gn)) {
  cv_gn[, total_af_raw := numeric_column(AF)]
  cv_gn[is.na(plot_af_raw) & !is.na(total_af_raw) & total_af_raw == 0, plot_af_raw := 0]
}
cv_gn[, plot_af := pmin(plot_af_raw, 1 - plot_af_raw)]
cv_gn <- cv_gn[!is.na(plot_af)]

clinvar_all <- fread(opts$clinvar)[Assembly == "GRCh38"]
assert_columns(
  clinvar_all,
  c("Chromosome", "Start", "ReferenceAllele", "AlternateAllele", "ClinicalSignificance"),
  "ClinVar GRCh38"
)

matched_key <- paste(cv_gn$CHROM, cv_gn$POS, cv_gn$REF, cv_gn$ALT)
clinvar_all[, join_key := paste0(
  "chr", Chromosome, " ", Start, " ", ReferenceAllele, " ", AlternateAllele
)]
clinvar_unmatched <- clinvar_all[!join_key %in% matched_key]
clinvar_unmatched[, plot_af := 0]

clinvar_full <- rbind(
  cv_gn[, .(ClinicalSignificance, plot_af)],
  clinvar_unmatched[, .(ClinicalSignificance, plot_af)],
  fill = TRUE
)

clinvar_full <- clinvar_full[
  !is.na(ClinicalSignificance) &
    ClinicalSignificance != "NA" &
    ClinicalSignificance != "-"
]

clinvar_top_n <- as.integer(opts$clinvar_top_n)
top_labels <- names(tail(sort(table(clinvar_full$ClinicalSignificance)), clinvar_top_n))
top_labels <- setdiff(
  top_labels,
  c("drug response", "not provided", "drug_response", "not_provided")
)
clinvar_sub <- clinvar_full[ClinicalSignificance %in% top_labels]

clinvar_sub[, category := fcase(
  ClinicalSignificance %in% c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"), "P/LP",
  ClinicalSignificance %in% c("Conflicting classifications of pathogenicity", "Uncertain significance"), "Conflicting/VUS",
  ClinicalSignificance %in% c("Benign", "Likely benign", "Benign/Likely benign"), "Benign/LB",
  default = "Other"
)]
clinvar_sub[, category := factor(category, levels = c("Benign/LB", "Conflicting/VUS", "P/LP", "Other"))]

af_label <- if (plot_af_col == "AF_GRPMAX") {
  "gnomAD v4 AF_GRPMAX"
} else {
  "gnomAD v4 Allele Frequency"
}

clinvar_panel <- ggplot(clinvar_sub, aes(x = plot_af, fill = category)) +
  geom_histogram(color = "gray60", alpha = 0.85, bins = 30, position = "stack", linewidth = 0.2) +
  theme_fig +
  xlab(af_label) +
  ylab("Count of ClinVar variants") +
  annotate(
    "text",
    x = 6.4e-6,
    y = Inf,
    label = "minimum\nobservable\nAF",
    size = 2.5,
    hjust = 0.5,
    vjust = 1.4,
    color = "gray40"
  ) +
  ggtitle(paste(opts$gnomad_label, "allele frequency for ClinVar variants")) +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 1e-6, base = 10),
    expand = c(0, 0),
    breaks = axis_breaks
  ) +
  scale_fill_viridis_d(option = "A", direction = -1) +
  labs(fill = "ClinVar classification")

clinvar_nonzero <- clinvar_sub[plot_af > 0]
if (nrow(clinvar_nonzero) > 0) {
  clinvar_inset <- ggplot(clinvar_nonzero, aes(x = plot_af, fill = category)) +
    geom_histogram(color = "gray60", alpha = 0.85, bins = 30, position = "stack", linewidth = 0.2) +
    theme_bw(base_size = 7) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8),
      plot.background = element_rect(color = "black", fill = NA, linewidth = 0.3)
    ) +
    xlab(af_label) +
    ylab("Count") +
    ggtitle("AF > 0") +
    scale_x_continuous(
      trans = pseudo_log_trans(sigma = 1e-6, base = 10),
      expand = c(0, 0),
      breaks = c(0, 1e-5, 1e-3, 1e-1, 0.5)
    ) +
    scale_fill_viridis_d(option = "A", direction = -1)

  clinvar_panel <- clinvar_panel +
    inset_element(clinvar_inset, left = 0.23, bottom = 0.15, right = 0.99, top = 0.99)
}

figure <- gwas_panel / clinvar_panel +
  plot_annotation(title = paste(opts$gwas_release, "and", opts$clinvar_label))

pdf_path <- file.path(opts$outdir, paste0(opts$prefix, ".pdf"))
png_path <- file.path(opts$outdir, paste0(opts$prefix, ".png"))
summary_path <- file.path(opts$outdir, paste0(opts$prefix, "_summary.tsv"))

ggsave(pdf_path, figure, width = 6.1, height = 5.9, units = "in")
ggsave(png_path, figure, width = 6.1, height = 5.9, units = "in", dpi = 300)

summary <- data.table(
  metric = c(
    "gwas_significant_associations_with_risk_af",
    "clinvar_joined_gnomad_rows",
    "clinvar_unmatched_grch38_rows_added_as_af0",
    "clinvar_rows_plotted",
    "clinvar_af_column"
  ),
  value = c(
    nrow(gwas_sig),
    nrow(cv_gn),
    nrow(clinvar_unmatched),
    nrow(clinvar_sub),
    plot_af_col
  )
)
fwrite(summary, summary_path, sep = "\t")

message("Saved PDF: ", pdf_path)
message("Saved PNG: ", png_path)
message("Saved summary: ", summary_path)
