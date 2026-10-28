# EA_additional_analyses.R
# This file contains the settings, analysis functions, and executable run sections.
# In RStudio, first run all function definitions above A0, then select:
#   A0 -> A1 -> A2 for the two-stage bootstrap (A1 saves/reuses simulation scores).
#   A0 -> A3 for Table S7 from completed simulations.
# Sourcing this entire file runs BOTH analyses, including A1's simulations.
# A0 loads the engine beside this script; no file-selection dialog is required.
# Run_confirmation.R is not needed for these analyses.
#
# Analyses for the five-variable correlated module and moving-average design:
# 1. Two-stage paired bootstrap of Monte Carlo calibrated EA minus each baseline.
# 2. Table S7: analytic EA versus the strongest Monte Carlo calibrated baseline.
# No packages beyond the existing simulation engine are required.

additional_analysis_settings <- function() {
  strengths <- c(3, 3.5, 4)
  Sigma_ma <- diag(50)
  Sigma_ma[abs(row(Sigma_ma) - col(Sigma_ma)) == 1L] <- 20 / 41
  list(
    cluster5_r075 = list(
      Sigma = make_block_cov(c(5, 45), c(.75, 0)), n = 499L,
      alts = lapply(strengths, function(a) list(
        id = paste0("a", a), mu = c(a, rep(a / 2, 4), rep(0, 45)),
        design = "amplitude", value = a
      )),
      stage = "additional_setting", source_id = "cluster5_r075",
      master_seed = 20261026L, index = 1L
    ),
    moving_average = list(
      Sigma = Sigma_ma, n = 499L,
      alts = lapply(strengths, function(a) {
        mu <- numeric(50)
        mu[c(25, 26, 27)] <- c(a / 3, a, a / 3)
        list(id = paste0("a", a), mu = mu, design = "amplitude", value = a)
      }),
      stage = "moving_average", source_id = "moving_average",
      source_dir = NA_character_, master_seed = 20261028L, index = 1L
    )
  )
}

.ea_analysis_labels <- function() c(
  cluster5_r075 = "Correlated module", moving_average = "Moving average peak"
)

# Retain the original scores. The dedicated cache avoids overwriting earlier runs.
# First execution simulates; resume = TRUE subsequently reuses compatible caches.
# Defaults reproduce the registered design seeds and the original draw order.
run_additional_score_simulations <- function(outdir, B = 19999L, R = 10000L,
    workers = 1L, backend = "auto", resume = TRUE, progress = TRUE) {
  settings <- additional_analysis_settings()
  results <- setNames(vector("list", length(settings)), names(settings))
  for (id in names(settings)) {
    results[[id]] <- run_setting(
      setting = settings[[id]], id = id, B = B, R = R,
      outdir = outdir, backend = backend, alpha = .05,
      save_scores = TRUE, resume = resume, progress = progress, workers = workers
    )
  }
  results
}

# Obtain existing summaries for Table S7 without launching simulations.
# Explicit in-memory objects take precedence, then the score cache, then the
# files written by the earlier five-variable and moving-average runners.
load_additional_analysis_results <- function(engine_file,
    five_variable_result = NULL, moving_average_result = NULL,
    score_dir = file.path(dirname(engine_file), "additional_analysis_results", "scores")) {
  results <- list(cluster5_r075 = five_variable_result,
                  moving_average = moving_average_result)
  saved_files <- c(
    cluster5_r075 = file.path(dirname(engine_file), "five_variable_results",
                             "five_variable_result.rds"),
    moving_average = file.path(dirname(engine_file), "moving_average_results",
                               "moving_average_result.rds")
  )
  for (id in names(results)) {
    if (!is.null(results[[id]])) next
    candidates <- c(file.path(score_dir, "groups", paste0(id, "_result.rds")),
                    saved_files[[id]])
    available <- candidates[file.exists(candidates)]
    if (!length(available)) stop(
      "No saved result for ", id, ". Run its simulation section or section A1 ",
      "of EA_additional_analyses.R first. No simulation was started automatically."
    )
    results[[id]] <- readRDS(available[1L])
  }
  results
}

.ea_check_analysis_results <- function(results, require_scores = FALSE) {
  ids <- names(.ea_analysis_labels())
  if (!is.list(results) || !identical(names(results), ids))
    stop("Supply a named list in order: cluster5_r075, moving_average.")
  for (id in ids) {
    x <- results[[id]]
    if (is.null(x$config) || is.null(x$power) || is.null(x$analytic_EA))
      stop("Incomplete result for ", id)
    if (!identical(x$config$setting$source_id, id))
      stop("The result is assigned to the wrong design: ", id)
    alts <- x$config$setting$alts
    alt_ids <- vapply(alts, `[[`, character(1), "id")
    if (!identical(alt_ids, c("a3", "a3.5", "a4")))
      stop("Expected amplitudes 3, 3.5, and 4 for ", id)
    if (any(x$power$B != x$config$B | x$power$R != x$config$R) ||
        any(x$analytic_EA$R != x$config$R))
      stop("Inconsistent replicate counts for ", id)
    if (require_scores) {
      labels <- c(paste0("known:", manuscript_methods),
                  paste0("unknown:", setdiff(manuscript_methods, "IHC")))
      if (is.null(x$bank) || is.null(x$scores)) stop(
        "Individual scores are missing for ", id, ". Run section A1 with ",
        "save_scores = TRUE; summary tables cannot support this bootstrap."
      )
      if (!identical(dim(x$bank), c(as.integer(x$config$B), length(labels))) ||
          !identical(dim(x$scores), c(as.integer(x$config$R), 4L, length(labels))) ||
          !identical(colnames(x$bank), labels) ||
          !identical(dimnames(x$scores)[[2]], c("null", alt_ids)) ||
          !identical(dimnames(x$scores)[[3]], labels) ||
          anyNA(x$bank) || anyNA(x$scores))
        stop("Missing or inconsistent score dimensions/labels for ", id)
    }
  }
  invisible(TRUE)
}

# Recompute thresholds from resampled calibration rows, using the engine's
# k = B + 1 - floor(alpha * (B + 1)) and strict score > cutoff rejection rule.
# All methods/panels share calibration indices. All methods/panels/amplitudes
# share evaluation indices. The two index samples are independent.
# Evaluation row multiplicities are used as weights; this is exactly equivalent
# to copying repeated rows. Null evaluation rows are excluded from the 90 tests.
two_stage_paired_bootstrap <- function(results, nboot = 2000L,
    seed = 20261030L, family_alpha = .05, outdir = NULL, progress = TRUE) {
  .ea_check_analysis_results(results, require_scores = TRUE)
  stopifnot(length(nboot) == 1L, is.finite(nboot), nboot >= 2L,
            nboot == floor(nboot), length(seed) == 1L, is.finite(seed),
            seed >= 0, seed <= .Machine$integer.max, seed == floor(seed),
            length(family_alpha) == 1L, is.finite(family_alpha),
            family_alpha > 0, family_alpha < 1)

  # Keep the caller's simulation RNG state unchanged, including if an error occurs.
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  initial_rng <- .Random.seed
  estimates <- replicates <- list()
  design_labels <- .ea_analysis_labels()

  for (id in names(results)) {
    x <- results[[id]]
    B <- x$config$B; R <- x$config$R
    alts <- x$config$setting$alts
    alt_ids <- vapply(alts, `[[`, character(1), "id")
    amplitudes <- vapply(alts, `[[`, numeric(1), "value")
    K <- length(alt_ids)
    labels <- colnames(x$bank)
    # Matrix columns run through amplitudes within each method/panel.
    evaluation <- matrix(x$scores[, alt_ids, , drop = FALSE], nrow = R)
    cut <- rank_cutoff(x$bank, alpha = x$config$alpha)
    point_power <- colMeans(sweep(evaluation, 2, rep(cut, each = K), ">"))
    rows <- list(); ea_columns <- other_columns <- integer(0)
    for (panel in c("known", "unknown")) {
      ee <- match(paste0(panel, ":EA"), labels)
      competitors <- which(startsWith(labels, paste0(panel, ":")))
      competitors <- setdiff(competitors, ee)
      for (a in seq_len(K)) for (j in competitors) {
        ec <- (ee - 1L) * K + a; oc <- (j - 1L) * K + a
        ea_columns <- c(ea_columns, ec); other_columns <- c(other_columns, oc)
        rows[[length(rows) + 1L]] <- data.frame(
          setting = unname(design_labels[id]), group = id, a = amplitudes[a],
          scenario = alt_ids[a], panel = panel,
          comparator = sub("^[^:]+:", "", labels[j]),
          EA_power = point_power[ec], comparator_power = point_power[oc],
          delta = point_power[ec] - point_power[oc], B = B, R = R,
          stringsAsFactors = FALSE
        )
      }
    }
    tab <- do.call(rbind, rows)
    boot <- matrix(NA_real_, nrow = nboot, ncol = nrow(tab))
    if (progress) cat("Bootstrap", id, "0 /", nboot, "\n")
    for (b in seq_len(nboot)) {
      cal_indices <- sample.int(B, B, replace = TRUE)
      eval_indices <- sample.int(R, R, replace = TRUE)
      boot_cut <- rank_cutoff(x$bank[cal_indices, , drop = FALSE], x$config$alpha)
      weights <- tabulate(eval_indices, nbins = R)
      rejected <- sweep(evaluation, 2, rep(boot_cut, each = K), ">")
      boot_power <- colSums(rejected * weights) / R
      boot[b, ] <- boot_power[ea_columns] - boot_power[other_columns]
      if (progress && (b == nboot || b %% max(1L, floor(nboot / 10)) == 0L))
        cat("Bootstrap", id, b, "/", nboot, "\n")
    }
    colnames(boot) <- paste(id, tab$panel, tab$scenario, tab$comparator, sep = "|")
    estimates[[id]] <- tab
    replicates[[id]] <- boot
  }
  bounds <- do.call(rbind, estimates); rownames(bounds) <- NULL
  draws <- do.call(cbind, replicates)
  family_size <- nrow(bounds) # 2 settings * 3 amplitudes * (8 + 7) = 90.
  z <- qnorm(1 - family_alpha / family_size) # One-sided Bonferroni bound.
  bounds$bootstrap_se <- apply(draws, 2, sd)
  bounds$lower_pointwise <- bounds$delta - qnorm(1 - family_alpha) * bounds$bootstrap_se
  bounds$lower_simultaneous <- bounds$delta - z * bounds$bootstrap_se
  bounds$positive_lower_bound <- bounds$lower_simultaneous > 0
  bounds$family_size <- family_size
  weakest <- which.min(bounds$lower_simultaneous)
  summary <- data.frame(
    nboot = nboot, bootstrap_seed = seed, family_alpha = family_alpha,
    family_size = family_size, z = z,
    positive_bounds = sum(bounds$positive_lower_bound),
    all_positive = all(bounds$positive_lower_bound),
    minimum_lower_bound = bounds$lower_simultaneous[weakest],
    weakest_setting = bounds$setting[weakest], weakest_a = bounds$a[weakest],
    weakest_panel = bounds$panel[weakest], weakest_comparator = bounds$comparator[weakest]
  )
  answer <- list(bounds = bounds, summary = summary, bootstrap_differences = draws,
    config = list(nboot = nboot, seed = seed, family_alpha = family_alpha,
                  family_size = family_size, z = z,
                  EA_calibration = "Monte Carlo rank cutoff recomputed in each resample",
                  simulation_configs = lapply(results, `[[`, "config")),
    RNG = list(kind = RNGkind(), initial = initial_rng, final = .Random.seed))
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    write.csv(bounds, file.path(outdir, "bootstrap_power_differences.csv"), row.names = FALSE)
    write.csv(summary, file.path(outdir, "bootstrap_summary.csv"), row.names = FALSE)
    saveRDS(answer, file.path(outdir, "two_stage_bootstrap.rds"))
    writeLines(c(
      sprintf("Two-stage paired bootstrap: %d resamples; bootstrap seed %d.", nboot, seed),
      "Calibration and evaluation rows were resampled independently with replacement.",
      "Methods, covariance panels and evaluation amplitudes retained their row pairing.",
      "These differences concern Monte Carlo calibrated EA, not analytic EA.",
      sprintf("One-sided normal-approximation Bonferroni bounds: delta - %.8f * SE; family size %d.", z, family_size),
      sprintf("%d of %d lower bounds are positive; smallest lower bound = %.8f.",
              sum(bounds$positive_lower_bound), family_size, min(bounds$lower_simultaneous)),
      "The minimum and positivity are computed from this run, not assumed from manuscript text."
    ), file.path(outdir, "bootstrap_report.txt"))
  }
  answer
}

# Half-up rounding of a count/R estimate; numeric full precision is kept separately.
.ea_format_analysis_power <- function(x, R, digits = 3L) {
  counts <- round(x * R)
  rounded <- floor((2 * counts * 10^digits + R) / (2 * R)) / 10^digits
  sprintf(paste0("%.", digits, "f"), rounded)
}

# Summaries suffice: this function never simulates and never uses bootstrap draws.
# Baselines exclude both versions of EA. Ties are retained in the detailed CSV.
# Table S7's 'best' is selected from full-precision powers before any rounding.
make_analytic_EA_table <- function(results, outdir = NULL) {
  .ea_check_analysis_results(results)
  if (any(!vapply(results, function(x) isTRUE(all.equal(x$config$alpha, .05)), logical(1))))
    stop("Table S7 expects the analytic rule P_EA <= 0.05 and 5% calibrated baselines.")
  detail <- nulls <- list()
  design_labels <- .ea_analysis_labels()
  for (id in names(results)) {
    x <- results[[id]]
    for (panel in c("known", "unknown")) {
      null <- x$analytic_EA[x$analytic_EA$scenario == "null" &
                              x$analytic_EA$panel == panel, ]
      if (nrow(null) != 1L) stop("Missing/duplicated analytic EA null row for ", id)
      nulls[[length(nulls) + 1L]] <- data.frame(
        Setting = unname(design_labels[id]), panel = panel, R = null$R,
        null_rejection_rate = null$power, low = null$low, high = null$high
      )
      for (alt in x$config$setting$alts) {
        analytic <- x$analytic_EA[x$analytic_EA$scenario == alt$id &
                                   x$analytic_EA$panel == panel, ]
        methods <- setdiff(manuscript_methods, c("EA", if (panel == "unknown") "IHC"))
        competitors <- x$power[x$power$scenario == alt$id &
                                 x$power$panel == panel & x$power$method %in% methods, ]
        if (nrow(analytic) != 1L || anyDuplicated(competitors$method) ||
            !setequal(competitors$method, methods) || anyNA(competitors$power))
          stop("Missing/duplicated analytic or baseline results for ", id, " ", alt$id, " ", panel)
        best <- max(competitors$power)
        detail[[length(detail) + 1L]] <- data.frame(
          Setting = unname(design_labels[id]), group = id, a = alt$value,
          panel = panel, R = analytic$R, analytic_EA = analytic$power,
          best_baseline = best,
          best_baseline_method = paste(competitors$method[competitors$power == best], collapse = "; "),
          difference = analytic$power - best, stringsAsFactors = FALSE
        )
      }
    }
  }
  detail <- do.call(rbind, detail); nulls <- do.call(rbind, nulls)
  table <- formatted <- list()
  for (id in names(results)) for (a in c(3, 3.5, 4)) {
    kn <- detail[detail$group == id & detail$a == a & detail$panel == "known", ]
    un <- detail[detail$group == id & detail$a == a & detail$panel == "unknown", ]
    row <- data.frame(
      Setting = unname(design_labels[id]), a = a,
      known_analytic_EA = kn$analytic_EA, known_best_baseline = kn$best_baseline,
      estimated_analytic_EA = un$analytic_EA, estimated_best_baseline = un$best_baseline
    )
    table[[length(table) + 1L]] <- row
    row[3:6] <- lapply(row[3:6], .ea_format_analysis_power, R = results[[id]]$config$R)
    formatted[[length(formatted) + 1L]] <- row
  }
  table <- do.call(rbind, table); formatted <- do.call(rbind, formatted)
  formatted$a <- as.character(formatted$a)
  rows <- apply(formatted, 1, function(row) paste0(paste(row, collapse = " & "), " \\\\"))
  replicate_counts <- unique(detail$R)
  caption <- paste0("Power of analytically calibrated EA and the strongest Monte Carlo calibrated comparator. ",
    if (length(replicate_counts) == 1L)
      paste0("Each entry uses ", format(replicate_counts, big.mark = ",", scientific = FALSE, trim = TRUE), " power replicates.")
    else "Replicate counts are recorded in the accompanying detailed CSV.")
  latex <- c(
    "\\begin{table}\\centering", paste0("\\caption{", caption, "}\\label{tab:new-settings-analytic}"),
    "{\\small\\setlength{\\tabcolsep}{6pt}\\renewcommand{\\arraystretch}{1.12}",
    "\\begin{tabular}{lrrrrr}", "\\toprule",
    "Setting & $a$ & \\multicolumn{2}{c}{Known covariance} & \\multicolumn{2}{c}{Estimated covariance} \\\\",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
    " & & Analytic EA & Best baseline & Analytic EA & Best baseline \\\\",
    "\\midrule", rows[1:3], "\\midrule", rows[4:6],
    "\\bottomrule", "\\end{tabular}}", "\\end{table}"
  )
  answer <- list(table = table, formatted = formatted, detail = detail,
                 analytic_null_rates = nulls, latex = latex,
                 config = lapply(results, `[[`, "config"))
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    write.csv(table, file.path(outdir, "table_S7.csv"), row.names = FALSE)
    write.csv(formatted, file.path(outdir, "table_S7_formatted.csv"), row.names = FALSE)
    write.csv(detail, file.path(outdir, "table_S7_details.csv"), row.names = FALSE)
    write.csv(nulls, file.path(outdir, "analytic_EA_null_rates.csv"), row.names = FALSE)
    writeLines(latex, file.path(outdir, "table_S7.tex"))
    saveRDS(answer, file.path(outdir, "table_S7.rds"))
    writeLines(c(
      "Analytic EA uses P_EA <= 0.05; baselines retain their Monte Carlo calibrated thresholds.",
      sprintf("Analytic EA exceeds the strongest baseline in %d of %d rows/panels (unrounded estimates).",
              sum(detail$difference > 0), nrow(detail)),
      sprintf("The smallest point-estimate advantage is %.8f.", min(detail$difference)),
      "These are point-estimate comparisons, not simultaneous significance claims.",
      "Null rejection rates and Wilson intervals are in analytic_EA_null_rates.csv."
    ), file.path(outdir, "table_S7_report.txt"))
  }
  answer
}

# ============================================================================
# ADDITIONAL ANALYSES: TWO-STAGE BOOTSTRAP AND ANALYTIC EA TABLE S7
# The definitions above and the executable sections below are in this file.
# Run the function definitions first, then A0 and the desired analysis sections.
# Sourcing this entire file runs A0-A3, including simulations and the bootstrap.
# ============================================================================

# ---- A0. SETUP: run this before A1, A2, or A3 ----
# In RStudio select Session > Set Working Directory > To Source File Location.
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
.ea_runtime$source_file <- normalizePath(engine_file, mustWork = TRUE)

additional_analysis_outdir <- file.path(dirname(engine_file), "additional_analysis_results")
additional_score_outdir <- file.path(additional_analysis_outdir, "scores")

# ---- A1. PREPARE SCORES FOR THE BOOTSTRAP: first run simulates both designs ----
# Five-variable seed: 20261026 (calibration/evaluation 20271026/20271027).
# Moving-average seed: 20261028 (calibration/evaluation 20271028/20271029).
# Each design has B = 19999, R = 10000, and a = 3, 3.5, 4.
# The existing wrappers do not save individual scores. This call uses the same
# engine and designs with save_scores = TRUE and a separate output directory.
# Later executions reuse compatible saved score files when resume = TRUE.
# Keep backend/workers unchanged to reuse this cache. Set resume = FALSE to
# explicitly regenerate it with the same simulation seeds.
additional_score_results <- run_additional_score_simulations(
  outdir = additional_score_outdir,
  B = 19999L,
  R = 10000L,
  workers = 1L,
  backend = "auto",
  resume = TRUE,
  progress = TRUE
)

# ---- A2. RUN THE TWO-STAGE PAIRED BOOTSTRAP: 2000 resamples ----
# If reopening RStudio after A1 completed, run A0 and A2; the stored scores load.
# Resamples reuse existing scores; no new Gaussian datasets are generated here.
# This analysis concerns Monte Carlo calibrated EA minus each comparator.
# Bootstrap seed 20261030 is separate from the simulation seeds, and the
# function restores the caller's RNG state when it finishes.
if (!exists("additional_score_results")) {
  additional_score_results <- load_additional_analysis_results(
    engine_file, score_dir = additional_score_outdir
  )
}
additional_bootstrap <- two_stage_paired_bootstrap(
  results = additional_score_results,
  nboot = 2000L,
  seed = 20261030L,
  family_alpha = .05,
  outdir = file.path(additional_analysis_outdir, "bootstrap"),
  progress = TRUE
)

if (interactive()) View(additional_bootstrap$bounds)
print(additional_bootstrap$summary, row.names = FALSE, digits = 8)
cat("\nBootstrap results:",
    normalizePath(file.path(additional_analysis_outdir, "bootstrap")), "\n")
# Check positive_bounds and minimum_lower_bound. The earlier manuscript's
# claim is checked from this run; the reference minimum exceeds 0.0188.

# ---- A3. GENERATE TABLE S7: analytic EA and strongest calibrated baseline ----
# A0 + A3 can run independently of the bootstrap. No new simulation is launched.
# Use A1's results if present; otherwise use the earlier design result objects
# or their saved RDS files. Missing results produce an explanatory error.
if (exists("additional_score_results")) {
  analytic_table_results <- additional_score_results
} else {
  analytic_table_results <- load_additional_analysis_results(
    engine_file,
    five_variable_result = if (exists("five_variable_result")) five_variable_result else NULL,
    moving_average_result = if (exists("moving_average_result")) moving_average_result else NULL,
    score_dir = additional_score_outdir
  )
}

table_S7 <- make_analytic_EA_table(
  results = analytic_table_results,
  outdir = file.path(additional_analysis_outdir, "table_S7")
)

if (interactive()) View(table_S7$formatted)
if (interactive()) View(table_S7$analytic_null_rates)
cat(paste(table_S7$latex, collapse = "\n"), "\n")
cat("\nTable S7 CSV/LaTeX files:",
    normalizePath(file.path(additional_analysis_outdir, "table_S7")), "\n")
# table_S7.csv: numeric, full-precision powers.
# table_S7_formatted.csv and table_S7.tex: manuscript values, three decimals.
# table_S7_details.csv: identity of best baseline(s), differences, and replicate counts.
# analytic_EA_null_rates.csv: analytic EA's Type I error and Wilson intervals.
