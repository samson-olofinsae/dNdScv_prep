#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (!requireNamespace("dndscv", quietly = TRUE)) {
    install.packages("dndscv", repos="https://cloud.r-project.org")
  }
  library(dndscv)
})

# the vignette dataset (well-formed for human refdb)
data("dataset_simbreast", package = "dndscv")
out <- "R/demo/dndscv_demo_input.csv"
dir.create("R/demo", showWarnings = FALSE, recursive = TRUE)
write.csv(dataset_simbreast, out, row.names = FALSE)
message("Wrote: ", out)
