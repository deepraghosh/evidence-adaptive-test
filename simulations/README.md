# Simulations and manuscript reproduction

The model is `X ~ N_d(mu, Sigma)` and an independent
`A ~ Wishart_d(n, Sigma)`, with null hypothesis `mu = 0`. Here `n` is the
Wishart degrees of freedom, and the signal vectors specify the mean of
`X`. Setting K uses the known covariance in its scores. Setting W uses
the scatter matrix and quantities fitted from it.

## Files and execution order

| File | Purpose |
|---|---|
| `EA_run.R` | RStudio setup, table and figure functions, and example commands |
| `EA_methods.R` | EA and comparator scores |
| `EA_simulations.R` | Fixed design registry, conditional reference sampler, estimates, and Wilson intervals |
| `EA_fast.cpp` | Optional C++ acceleration |
| `EA_component_benchmarks.R` | Analytic EA and component comparisons for Tables S8 and S9 |
| `results/` | Saved manuscript runs and the Figure 1 and Figure S1 PDFs |

Open `EA_run.R` in RStudio, select **Session → Set Working Directory →
To Source File Location**, then click **Source**, or run the following
from the `simulations` folder.

```r
source("EA_run.R")
```

Sourcing defines the functions and the full design registry. The run
commands below are explicit. Viewing the tables or plots does not write
files. The archived conditional runs used `backend = "cpp"`, which
requires Rcpp and a C++ compiler. Use `backend = "R"` for the base R
implementation. `workers = 1L` runs serially. Long conditional runs can
be skipped by loading the saved results.

## Conditional calibration and budgets

For each observed $(X,A)$, the code fixes $A^*=A+XX^\top$. It draws
independent $Z_b\sim N_d(0,I_d)$ and $V_b\sim\chi^2_{n-d+1}$, and forms
$X_b^0=(A^*)^{1/2}Z_b/\sqrt{Z_b^\top Z_b+V_b}$ and
$A_b^0=A^*-X_b^0(X_b^0)^\top$. The conditional reference distribution
does not need the true covariance. All methods share the observed
dataset and its reference draws within a replicate. Estimated quantities
are refit on each reference draw. Larger scores give stronger evidence,
so the rank p-value is `(1 + sum(T_ref >= T_obs))/(B + 1)`.

All results use nominal level `.05`. The conditional table runs use
`R = 1000` outer replicates per setting and `B = 99` conditional draws
per observation. The figure runs use `R = 200` and `B = 99`. Table S7
applies the analytic EA cutoff to the same observations as Tables 4
and 5, and uses their conditionally calibrated comparator powers.
Tables S8 and S9 apply the analytic procedures at level `.05` using
`R = 1000` outer replicates.

## Which output supplies each manuscript item?

| Manuscript item | Calculation | Saved result |
|---|---|---|
| Table 1 | Fixed M1–M9 design parameters | `simulation_designs()` |
| Tables 2 and 3 | Sixteen fixed alternatives, known and estimated covariance | `results/fixed_R1000_B99.rds` |
| Tables 4 and 5 | Five-variable correlated module and moving-average peak | `results/fixed_R1000_B99.rds` |
| Table 6 | The two TCGA applications | [Real-data guide](../real%20data%20analysis/README.md) |
| Tables S1 and S2 | Wilson null rejection intervals for the 12 fixed groups, known and estimated covariance | `results/fixed_R1000_B99.rds` |
| Table S3 | Fourteen sparse-design settings | `simulation_designs()` |
| Tables S4 and S5 | Sparse-design power, known and estimated covariance | `results/broader_R1000_B99.rds` |
| Table S6 | Wilson null rejection intervals for the module and moving-average designs | `results/fixed_R1000_B99.rds` |
| Table S7 | Analytic EA and the strongest conditionally calibrated comparator | Replay the observation streams in `results/fixed_R1000_B99.rds` |
| Tables S8 and S9 | D1–D5 analytic power and Wilson null rejection intervals | `results/components_R1000.rds` |
| Figures 1 and S1 | Matched-pair and sparse-design power curves | `results/figure_R200_B99.rds` and `results/figures/` |

## Tables 2–5 and S6: run the main driver

Choose one of these ways to assign `result`.

```r
# Read the saved manuscript run.
result <- readRDS(file.path(bundle_dir, "results", "fixed_R1000_B99.rds"))
```

Or rerun the fixed designs. Use backend = "R" for base R.

```r
result <- run_simulations(
  designs, families = c("matched", "module", "moving_average"),
  R = 1000L, B = 99L, alpha = .05, master_seed = 20261118L,
  workers = 4L, backend = "cpp", bundle_dir = bundle_dir
)
```

The `matched` family contains 16 alternatives in 12 groups. It includes
the nine matched-pair alternatives M1–M9, four complementary alternatives
B1–B4, and three controls C1–C3. The `module` and `moving_average`
families use `d = 50`, `n = 499`, and `a = 3, 3.5, 4`. Within each
family, the three alternatives share a null setting because they have
the same covariance and degrees of freedom. Observations at different
alternative amplitudes are drawn separately. The value `"unknown"`
identifies estimated covariance in the code.

```r
tables <- make_tables(result)
Table2 <- tables[["2"]]
Table3 <- subset(tables[["3"]], select = -IHC)
Table4 <- tables[["4"]]
Table5 <- tables[["5"]]
View(Table2)
View(Table3)
View(Table4)
View(Table5)
```

In a terminal, use `print(Table2, row.names = FALSE)`, and likewise
for the other tables. For Tables 4 and 5 the `setting` column is the
amplitude and the `panel` column identifies the covariance setting.
IHC is evaluated only in the known-covariance panel.

## Tables S1, S2 and S6: null rejection intervals

Use the same `result`. The object `result$type1` keeps point estimates
and the lower and upper pointwise 95% Wilson endpoints at full
precision. The following commands put those intervals in manuscript
order. Displayed endpoints use three decimal places.

```r
null_methods <- c(EA = "EA", ZM = "Z-Max", ZK = "Z-TopK", IHC = "IHC",
                  FC = "Feng-COM", QCQ = "QCQ", ZWS = "ZWS",
                  MP = "MinP", Q = "Quadratic")
null_table <- function(rows, settings, methods = null_methods) {
  out <- data.frame(Setting = settings)
  for (label in names(methods)) {
    z <- rows[rows$method == methods[[label]], ]
    i <- match(settings, z$setting)
    entries <- sprintf("[%.3f,%.3f]", z$lower[i], z$upper[i])
    out[[label]] <- ifelse(is.na(i), NA_character_, entries)
  }
  out
}

fixed_labels <- c("M1", "M2", "M3", "M4", "M5/M6/M7", "M8", "M9",
                  "B1/B2", "B3/B4", "C1", "C2", "C3")
TableS1 <- null_table(subset(result$type1, family == "matched" & panel == "known"),
                     fixed_labels)
TableS2 <- null_table(subset(result$type1, family == "matched" & panel == "unknown"),
                     fixed_labels, null_methods[names(null_methods) != "IHC"])

additional_labels <- c("Correlated module", "Moving average peak")
additional <- subset(result$type1, family %in% c("module", "moving_average"))
TableS6 <- rbind(
  cbind(Panel = "Known", null_table(subset(additional, panel == "known"),
                                   additional_labels)),
  cbind(Panel = "Estimated", null_table(subset(additional, panel == "unknown"),
                                       additional_labels))
)
View(TableS1)
View(TableS2)
View(TableS6)
```

## Tables 1 and S3: design parameters

The design registry `designs <- simulation_designs()` supplies all
means, covariance matrices, and degrees of freedom. To view the
parameters displayed in the two design tables, run the following commands.

```r
matched_ids <- paste0("M", 1:9)
matched_groups <- c("M1", "M2", "M3", "M4", "M5", "M5", "M5", "M8", "M9")
matched_offsets <- c(0, 0, 0, 0, 0, -.05, .05, 0, 0)
Table1 <- do.call(rbind, lapply(seq_along(matched_ids), function(i) {
  g <- designs$groups[[matched_groups[i]]]
  data.frame(Setting = matched_ids[i], d = nrow(g$Sigma), n = g$n,
             rho = g$Sigma[1, 2], a = g$means[[matched_ids[i]]][1],
             eta = matched_offsets[i])
}))

points <- subset(designs$display, study == "broader")
points <- points[match(paste0("O", 1:14), points$setting), ]
TableS3 <- do.call(rbind, lapply(seq_len(nrow(points)), function(i) {
  point <- points[i, ]
  g <- designs$groups[[point$group]]
  variable <- c(rho = "Correlation", dimension = "Dimension", amplitude = "Amplitude")
  data.frame(Setting = point$setting, Variable = unname(variable[point$sweep]),
             d = nrow(g$Sigma), n = g$n, rho = g$Sigma[1, 2],
             a = g$means[[point$setting]][1])
}))
View(Table1)
View(TableS3)
```

The alternative O5 is also the amplitude-two point in the figure.
It appears once in the design and power tables.

## Tables S4 and S5: broader-design power

```r
# Read the saved broader study.
broader <- readRDS(file.path(bundle_dir, "results", "broader_R1000_B99.rds"))
```

Or rerun it.

```r
broader <- run_simulations(
  designs, families = "broader",
  R = 1000L, B = 99L, alpha = .05, master_seed = 20261118L,
  workers = 4L, backend = "cpp", bundle_dir = bundle_dir
)
```

Then display the two panels.

```r
broader_tables <- make_tables(broader)
TableS4 <- broader_tables[["S4"]]
TableS5 <- subset(broader_tables[["S5"]], select = -IHC)
View(TableS4)
View(TableS5)
```

## Figure 1 and Figure S1: twelve individual panels

Use the saved `R = 200` estimates to match the plots in the manuscript.

```r
figure_result <- readRDS(file.path(bundle_dir, "results", "figure_R200_B99.rds"))

plot_figure(figure_result, 1)  # Choose 1 through 12.
plot_six_panels(figure_result, "main")     # Figure 1.
plot_six_panels(figure_result, "broader")  # Figure S1.
```

Panels 1–6 refer to the matched-pair comparison, and panels 7–12 refer
to the sparse-design comparison. In each group of six, known covariance
comes first. Each row varies correlation, dimension and signal strength,
in that order. Optional settings include `label_size`, `tick_size`,
`line_width`, `point_size`, and `show_legend`.

The following call regenerates the power estimates for the figures.

```r
figure_result <- run_simulations(
  designs, families = c("main_curves", "broader"),
  R = 200L, B = 99L, alpha = .05, master_seed = 20261118L,
  workers = 4L, backend = "cpp", bundle_dir = bundle_dir
)
```

These commands save vector versions of the two figures.

```r
pdf("Figure_1.pdf", width = 16, height = 10)
tryCatch(plot_six_panels(figure_result, "main"), finally = dev.off())

pdf("Figure_S1.pdf", width = 16, height = 10)
tryCatch(plot_six_panels(figure_result, "broader"), finally = dev.off())
```

The PDFs in `results/figures/` are the exact figure assets used in the manuscript.
Their saved run records the original driver and engine checksums.
Graphics devices and systems can differ in font metrics and PDF
metadata, so a regenerated PDF can differ in bytes.

## Table S7: analytic EA and the strongest comparator

Use the `result` from Tables 4 and 5. The code below replays the recorded
observation streams and applies `P_EA <= .05`. It reads the largest
non-EA power estimate in each covariance panel from the saved comparator
results. `B = 99` applies to those comparators.

```r
make_table_S7 <- function(result) {
  labels <- c(module = "Correlated module",
              moving_average = "Moving average peak")
  panels <- c("known", "unknown")
  if (!isTRUE(result$R == 1000L && result$B == 99L && result$alpha == .05) ||
      !all(names(labels) %in% result$power$family)) {
    stop("Use the completed Tables 4 and 5 result at R=1000, B=99, alpha=.05.")
  }

  old_kind <- RNGkind()
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  })
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")

  rows <- list()
  null_rows <- list()
  for (id in names(labels)) {
    g <- result$designs$groups[[id]]
    seed <- result$seeds[result$seeds$group == id, paste0("state", 1:7)]
    if (nrow(seed) != 1L) stop("Missing recorded seed for ", id)
    stream <- as.integer(unlist(seed, use.names = FALSE))
    root <- chol(g$Sigma)
    omega <- chol2inv(root)
    means <- c(list(null = numeric(nrow(g$Sigma))), g$means)
    counts <- matrix(0L, nrow = length(means), ncol = length(panels),
                     dimnames = list(names(means), panels))

    for (r in seq_len(result$R)) {
      for (setting in names(means)) {
        assign(".Random.seed", stream, envir = .GlobalEnv)
        Y <- draw_gaussian_rows(g$n + 1L, g$Sigma, root)
        A <- crossprod(sweep(Y, 2L, colMeans(Y), "-"))
        X <- drop(draw_gaussian_rows(1L, g$Sigma, root)) + means[[setting]]
        scores <- c(ea_statistic(X, Sigma = g$Sigma, inverse = omega),
                    ea_statistic(X, A = A, n = g$n))
        counts[setting, ] <- counts[setting, ] + (scores >= -log(result$alpha))

        # Reserve the conditional reference substream used by the comparator run.
        stream <- parallel::nextRNGSubStream(parallel::nextRNGSubStream(stream))
      }
    }

    null_rows[[id]] <- data.frame(
      Setting = labels[[id]], Known = counts["null", "known"] / result$R,
      Estimated = counts["null", "unknown"] / result$R
    )
    for (a in names(g$means)) {
      z <- subset(result$power, family == id & setting == a & method != "EA")
      known <- z$estimate[z$panel == "known"]
      unknown <- z$estimate[z$panel == "unknown"]
      if (!length(known) || !length(unknown)) {
        stop("Missing comparator powers for ", id, " at a = ", a)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        Setting = labels[[id]], a = as.numeric(a),
        Known_EA = counts[a, "known"] / result$R,
        Known_best = max(known),
        Estimated_EA = counts[a, "unknown"] / result$R,
        Estimated_best = max(unknown)
      )
    }
  }
  table <- do.call(rbind, rows)
  attr(table, "analytic_null") <- do.call(rbind, unname(null_rows))
  table
}

TableS7 <- make_table_S7(result)
View(TableS7)
View(attr(TableS7, "analytic_null"))
```

## Tables S8 and S9: analytic component comparisons

These five designs use the same marginal and precision-adjusted input
p-values for each comparator. The tables display EA, marginal and
precision-adjusted Bonferroni, pooled Bonferroni, and marginal Simes.

Read the saved `R = 1000` component analysis.

```r
components <- readRDS(file.path(bundle_dir, "results", "components_R1000.rds"))
```

Or rerun it with base R.

```r
source(file.path(bundle_dir, "EA_component_benchmarks.R"))
components <- run_component_benchmarks(R = 1000L, alpha = .05, seed = 20261119L)
```

Then prepare the two tables.

```r
methods <- c("EA", "Bonf_P", "Bonf_Q", "Pooled", "Simes_P")
columns <- c("design", "covariance", methods)
TableS8 <- components$power[, columns]
TableS8 <- TableS8[order(match(TableS8$covariance, c("known", "sample")),
                        TableS8$design), ]

TableS9 <- components$type1[, columns]
TableS9 <- TableS9[order(match(TableS9$covariance, c("known", "sample")),
                        TableS9$design), ]
z <- qnorm(.975)
total <- components$R
for (method in methods) {
  p <- TableS9[[method]]
  center <- (p + z^2 / (2 * total)) / (1 + z^2 / total)
  half <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) /
    (1 + z^2 / total)
  TableS9[[method]] <- sprintf("[%.4f,%.4f]",
                               pmax(0, center - half), pmin(1, center + half))
}
View(TableS8)
View(TableS9)
```

The value `"sample"` in the component result identifies estimated
covariance. Point estimates of the null rejection probabilities are
available in `components$type1`. To save any table to your chosen
directory, call `write.csv(TableS8, "Table_S8.csv", row.names = FALSE)`,
substituting the table object and file name as needed.

## Seeds, saved results and reproducibility

- The conditional runs use master seed `20261118` with
  L'Ecuyer-CMRG, Inversion and Rejection. A separate stream is assigned
  to each group in the full design registry. Within a group, the
  observed data and conditional reference draws use separate
  substreams. Filter via `families` to retain the original streams.
- The component study uses seed `20261119` with Mersenne-Twister,
  Inversion and Rejection.
- To inspect per-group seed states and runtime information, use
  `result$seeds` and `result$session` for the fixed run, and the
  corresponding fields of `broader` and `figure_result`. The figure run
  records the earlier driver and simulation-engine revisions that
  produced the archived estimates. Its score and C++ source checksums
  match the distributed source. The broader figure and tables use the same
  registered streams, with different values of `R`.
- `export_results(result, "my_results")` is an explicit export for a
  conditional run. Use the dedicated commands above for manuscript
  Wilson intervals and vector figures. In particular, export the
  manuscript plots from the `R = 200` figure result.
- Saved results provide reference estimates and do not enter score or
  cutoff calculations. Numerical libraries and graphics devices can
  introduce small cross-platform differences.
