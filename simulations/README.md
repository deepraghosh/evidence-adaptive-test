# Simulations and manuscript reproduction

This folder contains the four R scripts needed to reproduce the simulation analyses. Results are generated locally when you run the scripts; historical result files are not required or included. All manuscript runs use a nominal level of 0.05, **19,999 null calibration replicates** (`B`) and **10,000 evaluation replicates** (`R`) per simulation group.

The model is `X ~ N_d(mu, Sigma)` and an independent `A ~ Wishart_d(n, Sigma)`, with null hypothesis `mu = 0`. Here `n` is the Wishart degrees of freedom. The engine generates `n + 1` Gaussian rows to construct this experiment. Signal vectors describe the mean of `X`, rather than the mean of a single raw observation. Known-covariance panels use the true covariance. Estimated-covariance panels estimate their statistics from the scatter matrix; their Monte Carlo critical values use oracle null banks generated from the true simulation covariance.

## Files and execution order

| File | Purpose | Effect of sourcing the whole file |
|---|---|---|
| `EA_manuscript_simulations.R` | Shared statistics, design registry, simulation, table and plotting functions | Defines functions; does not start simulations |
| `Run_confirmation.R` | Confirmation designs, five-variable design, moving-average design and their null rejection table | Runs these simulations and writes their tables |
| `EA_explore_settings.R` | Editable correlation, dimension and amplitude sweeps | Runs six custom sweeps and displays twelve plots |
| `EA_additional_analyses.R` | Two-stage bootstrap and analytic EA Table S7 | Runs A0–A3, including score-generating simulations and the bootstrap |

Base R is sufficient. With `backend = "auto"`, an installed Rcpp package and working C++ compiler can accelerate the same calculations; the code falls back to base R when acceleration is unavailable. The scripts do not install packages. `workers = 1L` is the portable reproduction setting used throughout this guide. Larger runs can take substantial time, especially without acceleration.

In RStudio, open a script in this folder and choose **Session → Set Working Directory → To Source File Location**. Keep all four scripts together. Run the following setup before the code examples in this guide:

```r
engine_file <- normalizePath("EA_manuscript_simulations.R", mustWork = TRUE)
source(engine_file)
```

If a selection window opens in a driver, choose **EA_manuscript_simulations.R**, not the driver itself. Paths below are relative to this `simulations` folder.

For the main workflow:

1. Run `Run_confirmation.R` for Tables 2–5 and S6.
2. Run the fixed curve and original-design commands below for Figure 1, Figure S1, Table S2 and Tables S4/S5.
3. Run the desired sections of `EA_additional_analyses.R` for the 90 bootstrap comparisons and Table S7.

Use `EA_explore_settings.R` to explore changed input values; its defaults are deliberately separate from exact manuscript reproduction.

## Which output supplies each manuscript item?

| Manuscript item | Calculation | Main result location |
|---|---|---|
| Table 2 and 3 | Sixteen confirmation alternatives, known and estimated covariance | `confirmation_verification_results/` |
| Table 4 | Five-variable correlated module, both covariance panels | `five_variable_results/` |
| Tables 5 | Moving-average noise, both covariance panels | `moving_average_results/` |
| Figure 1 | Six main matched-pair curve panels | `figures/01_...png` through `06_...png` |
| Figure S1 | Six original sparse-design curve panels | `figures/07_...png` through `12_...png` |
| Table S1 and S2 | Wilson null rejection intervals for 21 groups, both covariance | `manuscript_tables/table_S2_...csv`, `table_S3.csv` |
| Tables S4 and S5 | Fourteen original sparse alternatives, known and estimated covariance | `manuscript_tables/table_S4.csv`, `table_S5.csv` |
| Table S6 | Null rejection estimates for the two additional designs | `type1_results/type1_error_table.csv` and `.tex` |
| Table S7 | Analytic EA and the strongest calibrated baseline | `additional_analysis_results/table_S7/` |
| Ninety simultaneous bounds | Two-stage paired bootstrap of calibrated EA minus each baseline | `additional_analysis_results/bootstrap/` |

## Tables 2–5 and S6: run the main driver

All paths below are relative to the repository root. Run the commands from that folder after completing the engine setup above.

Run the complete driver once:

```r
source("Run_confirmation.R")
confirmation_result <- result
```

Alternatively, select its sections separately in RStudio. The confirmation section produces `result`; the later sections produce `five_variable_result`, `moving_average_result` and `type1_table`. The later designs have separate output directories and do not overwrite the confirmation summaries. If running a design section alone in a new session, run the engine setup above first.

### Tables 2 and 3: sixteen confirmation alternatives

The CONFIRMATION section calls `run_manuscript(which = "confirmation", budget = "full", B = 19999L, R = 10000L, workers = 1L, backend = "auto", resume = TRUE, ...)`. It produces 16 alternatives across 12 calibration/evaluation groups: nine matched-pair alternatives M1–M9, four complementary alternatives B1–B4 and three controls C1–C3.

The full results are in `confirmation_verification_results/power.csv`; `power_wide.csv` places methods in columns. Non-null known-covariance rows supply **Table 2** and non-null estimated-covariance rows supply **Table 3**. The value `"unknown"` in the code identifies the estimated-covariance panel. The export instructions below put them in manuscript order and attach the M/B/C labels.

### Table 4: five-variable correlated module

The FIVE-VARIABLE section uses `d = 50`, `n = 499`, covariance `diag(R_5(0.75), I_45)`, and mean `a * (1, 0.5, 0.5, 0.5, 0.5, 0, ..., 0)` for `a = 3, 3.5, 4`. The registered master seed is `20261026`, giving calibration/evaluation seeds `20271026` and `20271027`.

The call is to `run_new_cluster()`. This wrapper always makes a fresh run and has no `resume` argument. The driver writes:

- `five_variable_results/five_variable_known_covariance.csv`
- `five_variable_results/five_variable_estimated_covariance.csv`
- Full-precision `power.csv`, `power_wide.csv`, `settings.csv` and `five_variable_result.rds`.

The first two CSVs provide the two panels of **Table 4**, with three decimal places. Missing estimated-covariance IHC is `NA`, not zero power.

### Table 5: moving-average noise

The MOVING-AVERAGE section uses `d = 50`, `n = 499`, covariance with diagonal 1, adjacent entries `20/41`, and all other entries zero. Its mean has `mu[26] = a`, `mu[25] = mu[27] = a/3`, and zero elsewhere, for `a = 3, 3.5, 4`. The registered master seed is `20261028`, giving calibration/evaluation seeds `20271028` and `20271029`.

The driver calls `run_moving_average()` with `resume = FALSE` for a fresh seeded run. It writes:

- `moving_average_results/moving_average_known_covariance.csv`
- `moving_average_results/moving_average_estimated_covariance.csv`
- Full-precision `power.csv`, `power_wide.csv`, `settings.csv` and `moving_average_result.rds`.

The first two CSVs provide the two panels of **Table 5**, with three decimal places. The engine uses a covariance square root to generate the stated moving-average distribution; replacing its generation algorithm with a different but distributionally equivalent one would change the observations produced by a fixed seed.

### Table S6: the two additional null rejection estimates

Run the final TYPE I ERROR section after completing the five-variable and moving-average sections. It extracts `scenario == "null"` from their results; it does not simulate again. In a restarted session, the section can read the two saved result RDS files after the engine setup.

It writes `type1_results/type1_error_table.csv` and `type1_results/type1_error_table.tex`. The table contains known- and estimated-covariance panels for both designs. Each estimate is the number of null rejections divided by 10,000. The amplitudes share a single null evaluation within each design; there are not three independent null rates per design.

## Fixed curve and original-design runs

Run these commands after the engine setup. They preserve the registered manuscript settings and seeds. `resume = TRUE` reuses compatible completed groups in the specified directories.

```r
curve_result <- run_manuscript(
  which = "curves", budget = "full", B = 19999L, R = 10000L,
  workers = 1L, backend = "auto", outdir = "curve_results", resume = TRUE
)
generate_tables(curve_result, "curve_results")

# The requested 14 sparse alternatives do not include the separate geometry study.
original_ids <- setdiff(names(manuscript_settings("original")), "geometry")
original_result <- run_manuscript(
  which = "original", budget = "full", group_ids = original_ids,
  B = 19999L, R = 10000L, workers = 1L, backend = "auto",
  outdir = "original_results", resume = TRUE
)
generate_tables(original_result, "original_results")
```

The curve run uses nine calibration groups. The original run above uses ten groups and 14 distinct alternatives. Excluding the geometry group does not renumber or change the remaining seeds.

If you restart R after completing these runs, load their summaries without rerunning simulations:

```r
read_result_summaries <- function(directory) {
  list(
    power = read.csv(file.path(directory, "power.csv"), check.names = FALSE),
    settings = read.csv(file.path(directory, "settings.csv"), check.names = FALSE)
  )
}
confirmation_result <- read_result_summaries("confirmation_verification_results")
curve_result <- read_result_summaries("curve_results")
original_result <- read_result_summaries("original_results")
```

These summaries suffice for the exports and plots below. The two-stage bootstrap additionally needs individual scores, which A1 saves separately.

## Export Tables 2 and 3, and Tables S4 and S5, in manuscript order

Run this helper and the two mapping blocks after obtaining `confirmation_result` and `original_result` above. It writes formatted four-decimal CSVs; full-precision summaries remain in the corresponding result folders.

```r
dir.create("manuscript_tables", showWarnings = FALSE)
method_columns <- c(
  EA = "EA", ZM = "Z-Max", ZK = "Z-TopK", IHC = "IHC",
  FC = "Feng-COM", QCQ = "QCQ", ZWS = "ZWS-{2,6}",
  MP = "Min-p", Q = "Quadratic"
)

export_power_table <- function(result, row_map, panel, filename) {
  methods <- method_columns
  if (panel == "unknown") methods <- methods[names(methods) != "IHC"]
  answer <- data.frame(Setting = row_map$Setting)
  for (column in names(methods)) {
    values <- vapply(seq_len(nrow(row_map)), function(i) {
      p <- result$power
      entry <- p[p$group == row_map$group[i] &
                   p$scenario == row_map$scenario[i] &
                   p$panel == panel & p$method == methods[[column]], ]
      if (nrow(entry) != 1L) stop("Missing or duplicated table cell: ",
                                 row_map$Setting[i], ", ", column)
      if (entry$R != 10000L || entry$B != 19999L)
        stop("These manuscript tables require B = 19999 and R = 10000.")
      entry$rejections / entry$R
    }, numeric(1))
    answer[[column]] <- sprintf("%.4f", values)
  }
  write.csv(answer, filename, row.names = FALSE)
  answer
}

confirmation_map <- data.frame(
  Setting = c(paste0("M", 1:9), paste0("B", 1:4), paste0("C", 1:3)),
  group = c(
    "matched_d20_r0.9_n99", "matched_d20_r0.95_n99",
    "matched_d50_r0.9_n249", "matched_d50_r0.95_n99",
    rep("matched_d50_r0.95_n249", 3), "matched_d50_r0.95_n499",
    "matched_d100_r0.95", rep("complementary_d50", 2),
    rep("complementary_d100", 2), "coherent_d50_r0.95",
    "identity_d50", "isolated_d50_r0.95"
  ),
  scenario = c(
    rep("matched", 5), "matched_t_offset_-0.05", "matched_t_offset_0.05",
    rep("matched", 2), rep(c("coherent_contrast_0.6", "coherent_contrast_0.8"), 2),
    "coherent", "isolated", "isolated"
  )
)
# Table 2: sixteen confirmation alternatives, known covariance.
table_2 <- export_power_table(confirmation_result, confirmation_map, "known",
                              "manuscript_tables/table_2.csv")
# Table 3: the same sixteen alternatives, estimated covariance.
table_3 <- export_power_table(confirmation_result, confirmation_map, "unknown",
                              "manuscript_tables/table_3.csv")

original_map <- data.frame(
  Setting = paste0("O", 1:14),
  group = c(paste0("cor_", c(.1, .3, .5, .7, .9)), rep("cor_0.9", 4),
            paste0("dim_", c(50, 100, 150, 200, 250))),
  scenario = c(paste0("rho=", c(.1, .3, .5, .7, .9)),
               paste0("a=", c(1, 1.5, 2.5, 3)),
               paste0("d=", c(50, 100, 150, 200, 250)))
)
table_S4 <- export_power_table(original_result, original_map, "known",
                               "manuscript_tables/table_S4.csv")
table_S5 <- export_power_table(original_result, original_map, "unknown",
                               "manuscript_tables/table_S5.csv")
if (interactive()) {
  View(table_2); View(table_3); View(table_S4); View(table_S5)
}
```

O1–O5 vary correlation at `d = 50`, `n = 249` and leading mean 2. O6–O9 vary leading mean at correlation 0.9, with O5 also serving as the amplitude-2 point. O10–O14 vary dimension from 50 to 250 with `n = 3*d - 1`; the active mean vector has leading value 3 and nine entries equal to 0.5. The dimension study changes both dimension and Wishart degrees of freedom. There are 14 distinct alternatives but 15 plotted points because O5 is shared by two sweeps.

## Tables S1 and S2: null rejection intervals for 21 groups

These tables combine the **12 confirmation null groups** and **nine main-curve null groups**. **Table S1** reports known-covariance intervals, and **Table S2** reports estimated-covariance intervals for the same 21 groups. Both report pointwise 95% Wilson intervals conditional on the realized calibration banks. These are different from simultaneous bootstrap intervals for power differences.

Run after the confirmation and curve runs. The following code writes a long numeric CSV and two formatted panel CSVs in manuscript order:

```r
null_map <- data.frame(
  Setting = c("M1", "M2", "M3", "M4", "M5/M6/M7", "M8", "M9",
              "B1/B2", "B3/B4", "C1", "C2", "C3", paste0("F", 1:9)),
  group = c(
    "matched_d20_r0.9_n99", "matched_d20_r0.95_n99",
    "matched_d50_r0.9_n249", "matched_d50_r0.95_n99",
    "matched_d50_r0.95_n249", "matched_d50_r0.95_n499",
    "matched_d100_r0.95", "complementary_d50", "complementary_d100",
    "coherent_d50_r0.95", "identity_d50", "isolated_d50_r0.95",
    "main_rho50", "main_rho70", "main_rho80", "main_rho90", "main_rho95",
    "main_d20", "main_d100", "main_d150", "main_d200"
  )
)
null_intervals <- subset(
  rbind(confirmation_result$power, curve_result$power), scenario == "null",
  select = c(group, panel, method, B, R, rejections, power, low, high)
)
stopifnot(length(unique(null_intervals$group)) == 21L,
          setequal(unique(null_intervals$group), null_map$group))
null_intervals$Setting <- null_map$Setting[match(null_intervals$group, null_map$group)]
null_intervals <- null_intervals[order(match(null_intervals$group, null_map$group),
                                      null_intervals$panel, null_intervals$method), ]
dir.create("manuscript_tables", showWarnings = FALSE)
write.csv(null_intervals, "manuscript_tables/tables_S1_S2_null_intervals.csv",
          row.names = FALSE)

for (panel in c("known", "unknown")) {
  methods <- method_columns
  if (panel == "unknown") methods <- methods[names(methods) != "IHC"]
  display <- data.frame(Setting = null_map$Setting)
  for (column in names(methods)) {
    display[[column]] <- vapply(null_map$group, function(id) {
      z <- null_intervals[null_intervals$group == id &
                            null_intervals$panel == panel &
                            null_intervals$method == methods[[column]], ]
      if (nrow(z) != 1L) stop("Missing or duplicated null interval.")
      sprintf("[%.3f, %.3f]", z$low, z$high)
    }, character(1))
  }
  table_id <- if (panel == "known") "S1" else "S2"
  write.csv(display,
            paste0("manuscript_tables/table_", table_id, "_", panel, ".csv"),
            row.names = FALSE)
}
```

M5/M6/M7 share a null group; B1/B2 share another and B3/B4 share another. F1–F5 correspond to correlations 0.50, 0.70, 0.80, 0.90 and 0.95; F6–F9 correspond to dimensions 20, 100, 150 and 200. F5 also provides the dimension-50 point and the amplitude sweep. The full long table contains 357 rows: 21 groups times (nine known-covariance methods plus eight estimated-covariance methods).

The exports from this section are:

- `manuscript_tables/table_S1_known.csv`: Table S1.
- `manuscript_tables/table_S2_unknown.csv`: Table S2.
- `manuscript_tables/tables_S1_S2_null_intervals.csv`: both tables in long numeric form, retaining rejection counts and interval endpoints.

## Figure 1 and Figure S1: twelve individual panels

Use the fixed `curve_result` and `original_result` objects above. These plotting commands do not run new simulations:

```r
figure_result <- list(
  power = rbind(curve_result$power, original_result$power),
  settings = rbind(curve_result$settings, original_result$settings)
)
generate_individual_figures(figure_result, outdir = "figures", figures = 1:12)
# To inspect one panel in RStudio, without saving or simulating:
if (interactive()) plot_manuscript_figure(1, figure_result)
```

| Panel number | Figure | Covariance | Horizontal axis |
|---|---|---|---|
| 1 | Figure 1 | Known | Correlation |
| 2 | Figure 1 | Known | Dimension |
| 3 | Figure 1 | Known | Amplitude |
| 4 | Figure 1 | Estimated | Correlation |
| 5 | Figure 1 | Estimated | Dimension |
| 6 | Figure 1 | Estimated | Amplitude |
| 7 | Figure S1 | Known | Correlation |
| 8 | Figure S1 | Known | Dimension |
| 9 | Figure S1 | Known | Amplitude |
| 10 | Figure S1 | Estimated | Correlation |
| 11 | Figure S1 | Estimated | Dimension |
| 12 | Figure S1 | Estimated | Amplitude |

Each PNG has a descriptive filename, for example `01_main_known_rho.png`. Supporting plotted-data CSVs, style CSVs and caption text are also written. The whiskers are pointwise 95% Wilson intervals conditional on the calibration bank. Lines connect the five sampled settings; they are not fitted power functions.

The main correlation sweep is `0.5, 0.7, 0.8, 0.9, 0.95`, at `d = 50`, `n = 249` and `a = 3.5`. Its dimension sweep is `20, 50, 100, 150, 200`, with `n = 5*d - 1`, correlation 0.95 and `a = 3.5`. Its amplitude sweep is `2.5, 3, 3.5, 4, 4.5`, at `d = 50`, `n = 249` and correlation 0.95. The original sweeps use the O1–O14 settings described above.

## Table S7 and the two-stage bootstrap

Open `EA_additional_analyses.R`. First select and run **all function definitions above A0**. Then run the desired executable sections:

| Section | Action |
|---|---|
| A0 | Load the engine and establish the output directories |
| A1 | Run both additional designs with their registered seeds, retaining individual scores |
| A2 | Perform 2,000 paired two-stage bootstrap resamples and write all 90 comparisons |
| A3 | Extract analytic EA powers and the strongest calibrated baseline to generate Table S7 |

Sourcing the entire file runs A0–A3, including simulations and bootstrap resampling. To obtain only Table S7 after `Run_confirmation.R` has completed, evaluate the definitions, then **A0 → A3**. A3 can load the existing five-variable and moving-average RDS files. It does not launch simulations or use bootstrap resamples.

### Table S7

For analytic EA, the rejection rule is `P_EA <= 0.05`; its cutoff is not estimated from the null bank. Every baseline retains its Monte Carlo calibrated cutoff. For each design, amplitude and covariance panel, A3 selects the largest **unrounded** non-EA power as the best baseline. The identity can change between rows; ties are retained in the detailed CSV.

Outputs under `additional_analysis_results/table_S7/` are:

- `table_S7.csv`: full-precision powers.
- `table_S7_formatted.csv` and `table_S7.tex`: manuscript table, three decimals.
- `table_S7_details.csv`: best-baseline identity or ties, powers and differences.
- `analytic_EA_null_rates.csv`: analytic EA null rejection rates and Wilson intervals.
- `table_S7.rds`: the complete table object.

Three-decimal power tables use half-up rounding of integer rejection counts. The analytic EA comparisons are point estimates; the A2 bootstrap bounds concern Monte Carlo calibrated EA and do not establish simultaneous bounds for Table S7.

### All 90 simultaneous lower bounds

Run the definitions and **A0 → A1 → A2**. A1 uses the same two designs, seeds and budgets as Tables 2 and 3, but retains individual calibration and evaluation scores. Ordinary summary tables cannot reconstruct the joint decisions needed for the bootstrap. A1 therefore has a separate score cache and may repeat these two simulations on its first run. Later compatible A1 runs reuse that cache.

For each bootstrap resample, A2 samples calibration rows jointly across methods and panels, recomputes each rank threshold, and independently resamples evaluation rows jointly across methods, panels and amplitudes. Sampling is with replacement from the stored scores; it does not generate new Gaussian observations. There are 90 comparisons: two designs times three amplitudes times eight known-covariance and seven estimated-covariance comparators.

For every original power difference `delta = power_EA - power_comparator`, the bootstrap SE is the sample standard deviation of its 2,000 resampled differences. The simultaneous one-sided lower bound is `delta - qnorm(1 - 0.05/90) * bootstrap_se`, with multiplier approximately 3.260767. Bonferroni adjustment targets approximate 95% simultaneous coverage under the bootstrap and normal approximation.

Inspect the results after A2:

```r
print(additional_bootstrap$summary, row.names = FALSE, digits = 10)
bounds <- additional_bootstrap$bounds
stopifnot(nrow(bounds) == 90L)
sum(bounds$lower_simultaneous > 0)
all(bounds$lower_simultaneous > 0)
min(bounds$lower_simultaneous)
bounds[which.min(bounds$lower_simultaneous),
       c("setting", "a", "panel", "comparator", "delta", "bootstrap_se", "lower_simultaneous")]
```

The recorded seeded run has 90 positive bounds and minimum **0.0188522653536516**, for the moving-average design at `a = 3.5` under estimated covariance, comparing EA with Z-Max. Thus all 90 bounds exceed **0.0188**. Reporting “the smallest is 0.0188” truncates that value; conventional rounding to four decimal places gives 0.0189. The code computes the minimum and positivity from each run rather than assuming the manuscript statement.

The full 90-row table is `additional_analysis_results/bootstrap/bootstrap_power_differences.csv`. Its `lower_simultaneous` column supports the paragraph. The summary is `bootstrap_summary.csv`, and the full object, including seeds and bootstrap differences, is `two_stage_bootstrap.rds`. To inspect a completed run after restarting R:

```r
additional_bootstrap <- readRDS("additional_analysis_results/bootstrap/two_stage_bootstrap.rds")
print(additional_bootstrap$summary, row.names = FALSE, digits = 10)
```

## Exploring changed settings

`EA_explore_settings.R` provides SETUP, six simulation blocks S1–S6, and twelve plot blocks. Evaluate SETUP, edit one five-value vector, run its S block, then either corresponding plot block. Each simulation supplies both covariance panels. Plotting again only redraws the existing result.

Defaults are `budget = "quick"` (`B = 199`, `R = 200`), `workers = 4`, `backend = "base"`, and custom seed base `2026101500`. For larger exploratory runs set `budget <- "full"`, `B <- 19999L`, `R <- 10000L`, `workers <- 1L`, and `backend <- "auto"` before executing the simulation blocks. `pilot` is also supported by this custom-sweep workflow.

The exploratory main-correlation vector is `0.5, 0.6, 0.7, 0.8, 0.9`, and custom-sweep seed allocation differs from the fixed manuscript registry. Increasing its budget alone does not reproduce the exact manuscript curves. Use the fixed commands earlier in this guide for Figure 1/S1. Exploratory results remain in R memory by default; plots can be saved through RStudio's Plots → Export. Sourcing the entire explorer runs all six custom sweeps.

## Seeds, saved results and reproducibility

The additional design master seeds are `20261026` and `20261028`, and the bootstrap seed is `20261030`. Fixed original, confirmation and curve settings have their own registered master seeds and indices. Calibration uses `master_seed + 10000 * index`; evaluation uses the next integer. The engine sets the R RNG kinds to Mersenne-Twister, Inversion and Rejection. The bootstrap restores the caller's RNG state on completion.

Use `settings_table(manuscript_settings("all"))` to inspect fixed-setting seeds. Result RDS objects also record simulation configuration and RNG states. The bootstrap saves its seed and RNG states in its own RDS object. File checks and recorded settings support reproducibility; no saved target values influence the test statistics or estimated powers.

Completed-group caches are used when an output directory and `resume = TRUE` are supplied. They require a matching design, budget, scoring implementation and backend. If you intentionally change a setting or implementation, choose a new output directory or explicitly request a fresh run. Interrupted groups are recomputed; already completed compatible groups can be reused. Avoid sourcing a complete driver when you intend only to inspect or reformat existing results.

Using the same scripts, seeds, budgets and computational environment reproduces the seeded experiments. Different R versions, numerical libraries or base/C++ backends may introduce small floating-point differences. Never replace generated numbers with the recorded examples when they differ; retain the configuration and inspect the cause.
