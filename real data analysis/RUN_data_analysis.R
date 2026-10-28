# Reproduce the current Table S9 from the original TCGA observations.
# Base R only. No internet, package installation, reconstruction, permutations,
# workers, or automatic result/figure saving. Run sections separately if desired.
# RStudio: open this file, choose Session > Set Working Directory > To Source
# File Location, then run the sections. source() and Rscript also work from
# another directory because the setup locates this saved script automatically.

# 1. SETUP ---------------------------------------------------------------
bundle_dir <- local({
  source_paths <- unlist(lapply(sys.frames(), function(f) f$ofile),
                         use.names = FALSE)
  if (length(source_paths)) {
    dirname(normalizePath(tail(source_paths, 1L), mustWork = TRUE))
  } else {
    script_arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(script_arg)) {
      script_path <- sub("^--file=", "", script_arg[1L])
      # Rscript may encode spaces as ~+~ in commandArgs().
      if (!file.exists(script_path))
        script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
      dirname(normalizePath(script_path, mustWork = TRUE))
    } else getwd()
  }
})
if (!file.exists(file.path(bundle_dir, "EA_unknown_covariance_realdata.R"))) {
  stop("Set the working directory to the real data analysis folder.")
}
data_dir <- file.path(bundle_dir, "data")
source(file.path(bundle_dir, "EA_manuscript_simulations.R"))
source(file.path(bundle_dir, "EA_unknown_covariance_realdata.R"))

# Check that the bundled input bytes match those used for the reported values.
input_checksums <- read.csv(file.path(bundle_dir, "INPUT_MD5.csv"),
                           stringsAsFactors = FALSE)
actual_md5 <- unname(tools::md5sum(file.path(data_dir, input_checksums$file)))
if (anyNA(actual_md5) || !identical(actual_md5, input_checksums$md5)) {
  stop("An input file is missing or has changed. Inspect INPUT_MD5.csv.")
}

# 2. COLORECTAL: CIN VERSUS GS --------------------------------------------
message("Calculating colorectal global p-values...")
colorectal_data <- load_tcga_unknown(data_dir, "coadread")
colorectal_groups <- tcga_cin_gs_groups(colorectal_data)
colorectal_result <- analytic_unknown_groups(colorectal_groups$Y1,
                                            colorectal_groups$Y2,
                                            backend = "base")
stopifnot(colorectal_result$m1 == 175L, colorectal_result$m2 == 45L,
          colorectal_result$d == 189L, colorectal_result$n == 218L)
print(colorectal_result$table[, c("method", "p_value_display", "log10_p_value")],
      row.names = FALSE)

# 3. STAD: N0 VERSUS N+ --------------------------------------------------
message("Calculating STAD global p-values...")
stad_data <- load_tcga_unknown(data_dir, "stad")
stad_groups <- tcga_nodal_groups(stad_data)
stad_result <- analytic_unknown_groups(stad_groups$Y1, stad_groups$Y2,
                                      backend = "base")
stopifnot(stad_result$m1 == 108L, stad_result$m2 == 230L,
          stad_result$d == 193L, stad_result$n == 336L)
print(stad_result$table[, c("method", "p_value_display", "log10_p_value")],
      row.names = FALSE)

# 4. COMBINED TABLE: THE SEVEN DISPLAYED METHODS ---------------------------
# Z-TopK and IHC remain NA in the full result objects, with reasons given.
make_display_table <- function(result, dataset) {
  z <- result$table[!is.na(result$table$log10_p_value), , drop = FALSE]
  data.frame(dataset = dataset, z, row.names = NULL, check.names = FALSE)
}
comparison_table <- rbind(make_display_table(colorectal_result, "Colorectal"),
                          make_display_table(stad_result, "STAD"))
print(comparison_table[, c("dataset", "method", "p_value_display")],
      row.names = FALSE)

# 5. VERIFY AGAINST THE INCLUDED REFERENCE VALUES -------------------------
# A small numerical tolerance allows floating-point differences between systems.
# This verifies reproducibility of the implemented calculation, not the accuracy
# of the asymptotic/provisional null-tail approximations.
expected <- read.csv(file.path(bundle_dir, "EXPECTED_VALUES.csv"),
                     stringsAsFactors = FALSE)
calculated_key <- paste(comparison_table$dataset, comparison_table$method)
reference_key <- paste(expected$dataset, expected$method)
ii <- match(reference_key, calculated_key)
if (anyNA(ii)) stop("The calculated and reference method labels do not match.")
verification <- data.frame(dataset = expected$dataset, method = expected$method,
  expected_log10_p = expected$log10_p_value,
  calculated_log10_p = comparison_table$log10_p_value[ii])
verification$absolute_difference <- abs(verification$calculated_log10_p -
                                        verification$expected_log10_p)
verification$matches <- verification$absolute_difference <= 1e-6
print(verification, row.names = FALSE)
if (all(verification$matches)) {
  message("All 14 p-values match the reference within 1e-6 on the log10 scale.")
} else {
  warning("Some calculations differ from the reference; inspect verification.")
}
message("Finished. Results remain in R memory; no result files were saved.")

# Useful objects to inspect after running:
# comparison_table
# colorectal_result$qcq_details
# colorectal_result$zws_details
# stad_result$qcq_details
# stad_result$zws_details
# colorectal_result$pooled_summaries
# stad_result$pooled_summaries
# verification
# The exceptionally small colorectal ZWS p-value is represented by its finite
# log10_p_value and a p_value_display character string. p_value is NA because
# the ordinary probability is below the floating-point range, not because the
# logarithmic calculation failed. See README.md for the calibration caveats.
