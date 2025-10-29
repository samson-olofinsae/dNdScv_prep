#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  # dNdScv is usually on CRAN; if not installed, guide user
  if (!requireNamespace("dndscv", quietly = TRUE)) {
    stop("The 'dndscv' package is not installed. In R: install.packages('dndscv')", call.=FALSE)
  }
  library(dndscv)
})

opt_list <- list(
  make_option(c("-i", "--input"),   type="character", help="Path to dNdScv input CSV (sampleID,chr,pos,ref,mut)", metavar="FILE"),
  make_option(c("-a", "--assembly"),type="character", default="auto",
              help="Assembly for refdb: hg19 | hg38 | auto (tries to infer) | custom (then provide --refdb_rds)", metavar="STR"),
  make_option(c("-r", "--refdb_rds"), type="character", default=NA,
              help="Path to custom refdb RDS (used if --assembly=custom)", metavar="FILE"),
  make_option(c("-o", "--outdir"),  type="character", default="results/dndscv",
              help="Output directory (default: results/dndscv)", metavar="DIR")
)

parser <- OptionParser(option_list = opt_list)
opt <- parse_args(parser)

if (is.null(opt$input)) stop("Missing --input")
if (!file.exists(opt$input)) stop(paste("Input not found:", opt$input))

# Create outdir
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading: ", opt$input)
m <- read.csv(opt$input, stringsAsFactors = FALSE)

# Light validation
required_cols <- c("sampleID","chr","pos","ref","mut")
if (!all(required_cols %in% colnames(m))) {
  stop(sprintf("Input CSV must have columns: %s", paste(required_cols, collapse=",")))
}

# Infer assembly if requested
infer_assembly <- function(df) {
  # crude heuristic: if any chr starts with "chr", assume hg38/hg19 style; else plain numeric
  has_chr_prefix <- any(grepl("^chr", df$chr, ignore.case = FALSE))
  # You might refine with lengths, but for now prefer hg38 on chr-prefix datasets
  if (has_chr_prefix) "hg38" else "hg19"
}

assembly <- opt$assembly
if (assembly == "auto") {
  assembly <- infer_assembly(m)
  message("Auto-inferred assembly: ", assembly)
}

# Reference DB
refdb <- NULL
if (assembly %in% c("hg19","hg38")) {
  refdb <- if (assembly == "hg19") refcds_hg19 else refcds_hg38
} else if (assembly == "custom") {
  if (is.na(opt$refdb_rds)) stop("With --assembly=custom you must provide --refdb_rds")
  if (!file.exists(opt$refdb_rds)) stop(paste("Custom refdb not found:", opt$refdb_rds))
  refdb <- readRDS(opt$refdb_rds)
} else {
  stop("Unknown --assembly. Use hg19 | hg38 | auto | custom")
}

# Run dNdScv
message("Running dNdScv ...")
dndsout <- dndscv(m, refdb = refdb)

# Save core tables
write.csv(dndsout$globaldnds, file=file.path(opt$outdir, "globaldnds.csv"), row.names = FALSE)
write.csv(dndsout$sel_cv,     file=file.path(opt$outdir, "sel_cv.csv"),     row.names = FALSE)
write.csv(dndsout$genesmuts,  file=file.path(opt$outdir, "genesmuts.csv"),  row.names = FALSE)

# Convenience: a table of significant genes at q<0.1 (common choice; adjust as needed)
if (!is.null(dndsout$sel_cv) && nrow(dndsout$sel_cv) > 0 && "qglobal_cv" %in% names(dndsout$sel_cv)) {
  signif_genes <- subset(dndsout$sel_cv, qglobal_cv < 0.1, select=c("gene_name","qglobal_cv"))
  rownames(signif_genes) <- NULL
  write.csv(signif_genes, file=file.path(opt$outdir, "significant_genes_q0.1.csv"), row.names = FALSE)
}

# Plots
pdf(file.path(opt$outdir, "dndscv_plots.pdf"), width=7, height=6)
  plot(dndsout)           # built-in diagnostic plots
dev.off()

message("Done. Outputs written to: ", normalizePath(opt$outdir))
