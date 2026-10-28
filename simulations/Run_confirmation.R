### =======
# Confirmation Results
### =======
# 1. In RStudio select Session > Set Working Directory > To Source File Location.
# Locate this saved script when sourced/Rscript; selected RStudio lines use getwd().
manuscript_dir <- local({
  paths <- unlist(lapply(sys.frames(), function(f) f$ofile), use.names = FALSE)
  if (length(paths)) {
    dirname(normalizePath(tail(paths, 1L), mustWork = TRUE))
  } else {
    arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(arg)) {
      path <- sub("^--file=", "", arg[1L])
      if (!file.exists(path)) path <- gsub("~+~", " ", path, fixed = TRUE)
      dirname(normalizePath(path, mustWork = TRUE))
    } else getwd()
  }
})
engine_file <- normalizePath(file.path(manuscript_dir, "EA_manuscript_simulations.R"), mustWork = TRUE)
source(engine_file)

# 2. Choose where the fresh results will be saved.
outdir <- file.path(dirname(engine_file), "confirmation_verification_results")

# 3. Actually RUN the simulations.
result <- run_manuscript(
  which = "confirmation",
  budget = "full",
  B = 19999L,
  R = 10000L,
  workers = 1L,
  backend = "auto",
  outdir = outdir,
  resume = TRUE
)

# 4. Arrange the results with methods in separate columns.
generate_tables(result, outdir)

# 5. Open the power tables in RStudio's data viewer.
power_tables <- read.csv(file.path(outdir, "power_wide.csv"), check.names = FALSE)

# Table 4: confirmation designs, known covariance
if (interactive()) View(subset(power_tables,scenario != "null" & panel == "known"))

# Table 5: confirmation designs, estimated covariance
if (interactive()) View(subset(power_tables, scenario != "null" & panel == "unknown"))

# Rejection probabilities and their 95% Wilson confidence intervals.
rejection_intervals <- subset(result$power, scenario == "null", select = c(group, scenario, panel, method,
             R, rejections, power, low, high))

# Known covariance.
if (interactive()) View(subset(rejection_intervals, panel == "known"))

# Estimated covariance.
if (interactive()) View(subset(rejection_intervals, panel == "unknown"))

# Save the intervals beside the other confirmation results.
write.csv(
  rejection_intervals,
  file.path(outdir, "rejection_intervals.csv"),
  row.names = FALSE
)


# ---- FIVE-VARIABLE DESIGN: SELECT AND RUN THIS SECTION SEPARATELY ----
# Select this section through its CSV-path prints to run only the five-variable design.
# Running the entire file runs confirmation, five-variable, moving-average and Table S6.
# Separate result and output-directory variables preserve the confirmation objects.

# Table 2: five-variable design, both covariance panels.
# The tables can be checked manually against the manuscript.

# 1. Select and load the complete simulation engine.
# When the file-selection window opens, choose EA_manuscript_simulations.R.
# Do NOT choose this Run_confirmation.R script.
# Choose the current complete engine containing run_new_cluster().
# five_variable_engine_file <- normalizePath(file.choose(), mustWork = TRUE)
# if (basename(five_variable_engine_file) != "EA_manuscript_simulations.R")
#   stop("Select EA_manuscript_simulations.R in the file-selection window.")
# source(five_variable_engine_file)

# 2. Save results beside the selected engine file.
five_variable_outdir <- file.path(dirname(engine_file), "five_variable_results")

# 3. Run the five-variable design at a = 3, 3.5, and 4.
# d = 50; n = 499 Wishart degrees, corresponding to 500 raw observations.
# Sigma = diag{R_5(0.75), I_45}; mu = a*(1, 0.5, 0.5, 0.5, 0.5, 0, ..., 0).
# Master seed 20261026 gives calibration/evaluation seeds 20271026/20271027.
# run_new_cluster() always computes a fresh run; it has no resume argument.
five_variable_result <- run_new_cluster(
  strengths = c(3, 3.5, 4),
  B = 19999L,
  R = 10000L,
  workers = 1L,
  backend = "auto",
  seed = 20261026L,
  engine_path = engine_file,
  progress = TRUE
)

# 4. Save full-precision summaries, seeds, and the complete five_variable_result object.
write_result_tables(five_variable_result, five_variable_outdir)
generate_tables(five_variable_result, five_variable_outdir)
saveRDS(five_variable_result, file.path(five_variable_outdir, "five_variable_result.rds"))

# 5. Format the two manuscript-style tables to exactly three decimal places.
# Full-precision values remain in five_variable_result$known_table, five_variable_result$estimated_table,
# power.csv, power_wide.csv, and the saved five_variable_result object.
# With R = 10000, integer rejection counts give exact half-up rounding.
# NA represents unavailable estimated-covariance IHC, not a zero power.
format_five_variable_power_table <- function(tab, R = 10000L) {
  stopifnot(R == 10000L)
  formatted <- tab
  for (method in setdiff(names(tab), "a")) {
    x <- tab[[method]]
    available <- !is.na(x)
    counts <- as.integer(round(x[available] * R))
    formatted[[method]] <- rep(NA_character_, length(x))
    formatted[[method]][available] <- sprintf(
      "%.3f", ((counts + 5L) %/% 10L) / 1000
    )
  }
  formatted
}

known_table <- format_five_variable_power_table(five_variable_result$known_table, R = five_variable_result$config$R)
estimated_table <- format_five_variable_power_table(five_variable_result$estimated_table, R = five_variable_result$config$R)
known_csv <- file.path(five_variable_outdir, "five_variable_known_covariance.csv")
estimated_csv <- file.path(five_variable_outdir, "five_variable_estimated_covariance.csv")

write.csv(known_table, known_csv, row.names = FALSE, na = "NA")
write.csv(estimated_table, estimated_csv, row.names = FALSE, na = "NA")

# 6. View the two tables and print the locations of their CSV files.
if (interactive()) View(known_table)
if (interactive()) View(estimated_table)
cat("\nKnown-covariance CSV:", normalizePath(known_csv), "\n")
cat("Estimated-covariance CSV:", normalizePath(estimated_csv), "\n")

# ---- TABLE 3: MOVING-AVERAGE DESIGN ----
# Moving-average design: run the simulation and export the power tables to CSV.
# The tables can be checked manually against the manuscript.

# 1. Load the updated engine if it has not already been loaded.
# source(engine_file)

# 2. Save results beside the selected engine file.
moving_average_outdir <- file.path(dirname(engine_file), "moving_average_results")

# 3. Run the moving-average design at a = 3, 3.5, and 4.
# d = 50; n = 499 Wishart degrees, corresponding to 500 raw observations.
# Sigma has unit diagonal, adjacent entries 20/41, and zero elsewhere.
# mu[26] = a; mu[25] = mu[27] = a/3; all other means are zero.
# Master seed 20261028 gives calibration/evaluation seeds 20271028/20271029.
# resume = FALSE forces a fresh simulation using these fixed seeds.
moving_average_result <- run_moving_average(
  strengths = c(3, 3.5, 4),
  B = 19999L,
  R = 10000L,
  workers = 1L,
  backend = "auto",
  seed = 20261028L,
  outdir = moving_average_outdir,
  resume = FALSE,
  progress = TRUE
)

# 4. Full-precision results are saved automatically by run_moving_average().
# These include power.csv, power_wide.csv, settings.csv,
# known_covariance.csv, estimated_covariance.csv,
# and moving_average_result.rds.
# No additional saving commands are needed for those files.

# 5. Format the two manuscript-style tables to exactly three decimal places.
# The original full-precision results remain unchanged.
# With R = 10000, integer rejection counts give exact half-up rounding.
# NA represents unavailable estimated-covariance IHC, not a zero power.
format_moving_average_power_table <- function(tab, R = 10000L) {
  stopifnot(R == 10000L)
  formatted <- tab
  
  for (method in setdiff(names(tab), "a")) {
    x <- tab[[method]]
    available <- !is.na(x)
    counts <- as.integer(round(x[available] * R))
    
    formatted[[method]] <- rep(NA_character_, length(x))
    formatted[[method]][available] <- sprintf(
      "%.3f", ((counts + 5L) %/% 10L) / 1000
    )
  }
  
  formatted
}

moving_average_known_table <- format_moving_average_power_table(
  moving_average_result$known_table,
  R = moving_average_result$config$R
)

moving_average_estimated_table <- format_moving_average_power_table(
  moving_average_result$estimated_table,
  R = moving_average_result$config$R
)

moving_average_known_csv <- file.path(moving_average_outdir, "moving_average_known_covariance.csv")

moving_average_estimated_csv <- file.path(moving_average_outdir, "moving_average_estimated_covariance.csv")

write.csv(moving_average_known_table, moving_average_known_csv, row.names = FALSE, na = "NA")

write.csv(moving_average_estimated_table, moving_average_estimated_csv, row.names = FALSE, na = "NA")

# 6. View the two tables and print the locations of their CSV files.
if (interactive()) View(moving_average_known_table)
if (interactive()) View(moving_average_estimated_table)


# ============================================================
# Table S6: Type I error for five-variable and moving-average designs
# Extract existing null results; no simulations are run here.
# ============================================================

# 1. Use the completed result objects.
# If RStudio was restarted, load the saved objects instead.
if (!exists("five_variable_result")) {
  five_variable_result <- readRDS(
    file.path(
      dirname(engine_file),
      "five_variable_results",
      "five_variable_result.rds"
    )
  )
}

if (!exists("moving_average_result")) {
  moving_average_result <- readRDS(
    file.path(
      dirname(engine_file),
      "moving_average_results",
      "moving_average_result.rds"
    )
  )
}

# 2. Map the engine's method names to the manuscript headings.
type1_methods <- c(
  EA  = "EA",
  ZM  = "Z-Max",
  ZK  = "Z-TopK",
  IHC = "IHC",
  FC  = "Feng-COM",
  QCQ = "QCQ",
  ZWS = "ZWS-{2,6}",
  MP  = "Min-p",
  Q   = "Quadratic"
)

# Extract one setting and covariance panel.
extract_type1 <- function(result, setting_label, panel_name) {
  null_rows <- subset(
    result$power,
    scenario == "null" & panel == panel_name
  )
  
  rates <- vapply(type1_methods, function(method_name) {
    entry <- null_rows[null_rows$method == method_name, ]
    
    # Estimated-covariance IHC is unavailable.
    if (nrow(entry) == 0L &&
        panel_name == "unknown" &&
        method_name == "IHC") {
      return(NA_real_)
    }
    
    if (nrow(entry) != 1L)
      stop("Missing or duplicated null result: ",
           setting_label, ", ", panel_name, ", ", method_name)
    
    if (entry$R != 10000L || entry$B != 19999L)
      stop("Expected B = 19999 and R = 10000.")
    
    entry$rejections / entry$R
  }, numeric(1))
  
  data.frame(
    Panel = if (panel_name == "known")
      "Known covariance" else "Estimated covariance",
    Setting = setting_label,
    as.list(rates),
    check.names = FALSE
  )
}

# 3. Combine rows in the same order as the manuscript table.
type1_table <- rbind(
  extract_type1(
    five_variable_result, "Correlated module", "known"
  ),
  extract_type1(
    moving_average_result, "Moving average peak", "known"
  ),
  extract_type1(
    five_variable_result, "Correlated module", "unknown"
  ),
  extract_type1(
    moving_average_result, "Moving average peak", "unknown"
  )
)

# Keep numeric values in type1_table.
# Create a display copy with exactly four decimal places.
type1_display <- type1_table

for (method in names(type1_methods)) {
  x <- type1_table[[method]]
  type1_display[[method]] <- ifelse(
    is.na(x), "---", sprintf("%.4f", x)
  )
}

# 4. View and save the combined table.
type1_outdir <- file.path(dirname(engine_file), "type1_results")
dir.create(type1_outdir, recursive = TRUE, showWarnings = FALSE)

type1_csv <- file.path(type1_outdir, "type1_error_table.csv")

write.csv(
  type1_display,
  type1_csv,
  row.names = FALSE
)

if (interactive()) View(type1_display)

# 5. Generate the LaTeX tabular environment.
latex_rows <- function(tab) {
  apply(
    tab[, c("Setting", names(type1_methods)), drop = FALSE],
    1,
    function(row) paste0(paste(row, collapse = " & "), " \\\\")
  )
}

type1_latex <- c(
  "\\begin{tabular}{lrrrrrrrrr}",
  "\\toprule",
  "Setting & EA & ZM & ZK & IHC & FC & QCQ & ZWS & MP & Q \\\\",
  "\\midrule",
  "\\multicolumn{10}{l}{\\emph{Known covariance}}\\\\",
  latex_rows(type1_display[1:2, ]),
  "\\midrule",
  "\\multicolumn{10}{l}{\\emph{Estimated covariance}}\\\\",
  latex_rows(type1_display[3:4, ]),
  "\\bottomrule",
  "\\end{tabular}"
)

type1_tex <- file.path(type1_outdir, "type1_error_table.tex")
writeLines(type1_latex, type1_tex)

# Print the LaTeX table for copying into the manuscript.
cat(paste(type1_latex, collapse = "\n"), "\n")
cat("\nCSV saved to:", normalizePath(type1_csv), "\n")
cat("LaTeX saved to:", normalizePath(type1_tex), "\n")
