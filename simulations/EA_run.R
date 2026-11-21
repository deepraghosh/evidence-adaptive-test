# Open this file in RStudio. Set the working directory to this file's folder.
# Source this file once, then run the examples at the bottom as needed.
bundle_dir <- local({
  paths <- unlist(lapply(sys.frames(), function(frame) frame$ofile), use.names = FALSE)
  if (length(paths)) dirname(normalizePath(tail(paths, 1L))) else getwd()
})
source(file.path(bundle_dir, "EA_methods.R"))
source(file.path(bundle_dir, "EA_simulations.R"))

table_titles <- c(
  "2" = "Matched-pair and control designs, known covariance",
  "3" = "Matched-pair and control designs, estimated covariance",
  "4" = "Unequal responses in a five-variable correlated module",
  "5" = "Localized peaks with moving average noise",
  "S1" = "Null rejection estimates, known covariance",
  "S2" = "Null rejection estimates, estimated covariance",
  "S4" = "Broader-study power, known covariance",
  "S5" = "Broader-study power, estimated covariance",
  "S6" = "Null rejection estimates for the two additional designs")

make_tables <- function(result) {
  estimates <- result$estimates
  groups <- result$designs$groups
  methods <- c(EA = "EA", ZM = "Z-Max", ZK = "Z-TopK", IHC = "IHC", FC = "Feng-COM",
               QCQ = "QCQ", ZWS = "ZWS", MP = "MinP", Q = "Quadratic")
  wide <- unique(estimates[, c("group", "family", "setting", "panel")])
  for (heading in names(methods)) {
    wide[[heading]] <- vapply(seq_len(nrow(wide)), function(i) {
      row <- wide[i, ]
      estimates$estimate[estimates$group == row$group & estimates$setting == row$setting &
                          estimates$panel == row$panel & estimates$method == methods[[heading]]]
    }, numeric(1))
  }
  tables <- list()
  add <- function(number, data) {
    if (!nrow(data)) return(invisible(NULL))
    tables[[number]] <<- data[, c("setting", "panel", names(methods)), drop = FALSE]
  }
  order_settings <- c(paste0("M", 1:9), paste0("B", 1:4), paste0("C", 1:3))
  for (panel in c("known", "unknown")) {
    z <- wide[wide$family == "matched" & wide$setting != "null" & wide$panel == panel, ]
    z <- z[order(match(z$setting, order_settings)), ]
    add(if (panel == "known") "2" else "3", z)
  }
  for (family in c("module", "moving_average")) {
    z <- wide[wide$family == family & wide$setting != "null", ]
    z <- z[order(match(z$panel, c("known", "unknown")), as.numeric(z$setting)), ]
    add(if (family == "module") "4" else "5", z)
  }
  labels <- c("M1", "M2", "M3", "M4", "M5/M6/M7", "M8", "M9", "B1/B2", "B3/B4",
              "C1", "C2", "C3", paste0("F", 1:9))
  for (panel in c("known", "unknown")) {
    z <- wide[wide$setting == "null" & wide$panel == panel, ]
    z$setting <- vapply(z$group, function(id) groups[[id]]$null_label, character(1))
    z <- z[z$setting %in% labels, ]
    z <- z[order(match(z$setting, labels)), ]
    add(if (panel == "known") "S1" else "S2", z)
    z <- wide[wide$family == "broader" & wide$setting != "null" & wide$panel == panel, ]
    z <- z[order(match(z$setting, paste0("O", 1:14))), ]
    add(if (panel == "known") "S4" else "S5", z)
  }
  z <- subset(wide, family %in% c("module", "moving_average") & setting == "null")
  z$setting <- vapply(z$group, function(id) groups[[id]]$null_label, character(1))
  z <- z[order(match(z$panel, c("known", "unknown")), match(z$group, c("module", "moving_average"))), ]
  add("S6", z)
  tables[names(table_titles)[names(table_titles) %in% names(tables)]]
}

curve_styles <- data.frame(
  method = method_names,
  color = c("#111111", "#0072B2", "#D55E00", "#6F4C9B", "#009E73", "#E69F00",
            "#CC79A7", "#56B4E9", "#666666"),
  point = c(18, 17, 25, 16, 15, 8, 7, 1, 4),
  line = c(1, 1, 1, 2, 2, 2, 3, 3, 4))

plot_curve <- function(result, study = "main", panel = "known", sweep = "rho",
                        metric = "power", label_size = 1.7, tick_size = 1.3,
                        line_width = 4, point_size = 1.3, show_legend = FALSE) {
  display <- result$designs$display
  display <- display[display$study == study & display$sweep == sweep, ]
  data <- result$estimates[result$estimates$panel == panel, ]
  if (metric == "type1") {
    display$setting <- "null"
  } else {
    stopifnot(metric == "power")
  }
  data <- merge(display, data, by = c("group", "setting"))
  if (!nrow(data)) stop("Run the requested design family first.")
  x <- sort(unique(display$x))
  span <- diff(range(x))
  xlabel <- switch(sweep, rho = "Correlations", dimension = "Number of tests", amplitude = "Signal strengths")
  ylabel <- if (metric == "power") "Simulated power" else "Simulated type-I error"
  plot(range(x) + c(-1, 1) * span * .035, c(-.03, 1.03), type = "n", axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  abline(h = seq(0, 1, .2), col = "#e7e9ec", lwd = .7)
  axis(1, at = x, labels = format(x, trim = TRUE), cex.axis = tick_size)
  axis(2, at = seq(0, 1, .2), las = 1, cex.axis = tick_size)
  mtext(xlabel, side = 1, line = 3.5, cex = label_size)
  mtext(ylabel, side = 2, line = 3.8, cex = label_size)
  box(bty = "l")
  available <- character()
  for (i in seq_len(nrow(curve_styles))) {
    style <- curve_styles[i, ]
    z <- data[data$method == style$method, ]
    z <- z[order(z$x), ]
    if (all(is.na(z$estimate))) next
    available <- c(available, style$method)
    lines(z$x, z$estimate, type = "b", col = style$color, bg = style$color,
          pch = style$point, lty = style$line, lwd = line_width, cex = point_size)
  }
  if (show_legend) {
    style <- curve_styles[curve_styles$method %in% available, ]
    legend("bottomright", legend = style$method, col = style$color, pch = style$point,
           lty = style$line, lwd = line_width, bty = "n", cex = .8)
  }
  invisible(data)
}

# Numbers 1--6 give the main curves, 7--12 the broader curves.
# Within each set, known correlation/dimension/amplitude come first,
# followed by estimated correlation/dimension/amplitude.
plot_figure <- function(result, number, ...) {
  stopifnot(length(number) == 1, number %in% 1:12)
  within_set <- (number - 1) %% 6 + 1
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  par(mar = c(6, 6.5, 1, 1), family = "sans")
  plot_curve(result, study = if (number <= 6) "main" else "broader",
    panel = if (within_set <= 3) "known" else "unknown",
    sweep = c("rho", "dimension", "amplitude")[(within_set - 1) %% 3 + 1], ...)
}

plot_six_panels <- function(result, study = "main", ...) {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  par(mfrow = c(2, 3), mar = c(6, 6.5, 1, 1), family = "sans", cex = 1)
  for (panel in c("known", "unknown")) {
    for (sweep in c("rho", "dimension", "amplitude")) {
      plot_curve(result, study = study, panel = panel, sweep = sweep, ...)
    }
  }
}

write_tables_tex <- function(tables, file) {
  output <- c("% Requires booktabs. Entries are simulated rejection proportions.",
              "% A dash denotes the original IHC being restricted to known covariance. No intervals are used.")
  for (number in names(tables)) {
    table <- tables[[number]]
    headings <- setdiff(names(table), c("setting", "panel"))
    if (all(table$panel == "unknown")) headings <- setdiff(headings, "IHC")
    label <- if (number %in% c("4", "5")) "$a$" else "Setting"
    power_table <- number %in% c("2", "3", "4", "5", "S4", "S5")
    caption <- paste0(table_titles[[number]], ".")
    if (power_table) caption <- paste(caption, "Bold denotes the largest point estimate among available methods.")
    caption <- paste(caption, "Conditional Monte Carlo calibration; IHC is retained only for known covariance.")
    output <- c(output, paste0("% Table ", number), "\\begin{table}[htbp]\\centering", "\\small",
      "\\setlength{\\tabcolsep}{3pt}\\renewcommand{\\arraystretch}{1.12}",
      paste0("\\caption{", caption, "}"),
      paste0("\\label{tab:conditional-", number, "}"),
      paste0("\\begin{tabular}{l", paste(rep("r", length(headings)), collapse = ""), "}"),
      "\\toprule", paste0(paste(c(label, headings), collapse = " & "), " \\\\"), "\\midrule")
    last_panel <- ""
    for (i in seq_len(nrow(table))) {
      if (length(unique(table$panel)) > 1 && table$panel[i] != last_panel) {
        panel <- if (table$panel[i] == "known") "Known covariance" else "Estimated covariance"
        output <- c(output, paste0("\\multicolumn{", length(headings) + 1,
                                   "}{l}{\\emph{", panel, "}}\\\\"))
        last_panel <- table$panel[i]
      }
      values <- vapply(table[i, headings], function(x) if (is.na(x)) "---" else sprintf("%.4f", x), "")
      if (power_table) {
        largest <- max(unlist(table[i, headings]), na.rm = TRUE)
        best <- which(!is.na(table[i, headings]) & table[i, headings] == largest)
        values[best] <- paste0("\\textbf{", values[best], "}")
      }
      output <- c(output, paste0(paste(c(table$setting[i], values), collapse = " & "), " \\\\"))
    }
    output <- c(output, "\\bottomrule", "\\end{tabular}", "\\end{table}", "")
  }
  writeLines(output, file)
}

# Explicit export is the only operation that saves results or figures.
export_results <- function(result, directory = "results") {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  tables <- make_tables(result)
  stacked <- do.call(rbind, lapply(names(tables), function(number) {
    data.frame(table = number, tables[[number]], check.names = FALSE)
  }))
  write.csv(stacked, file.path(directory, "tables.csv"), row.names = FALSE, na = "")
  write_tables_tex(tables, file.path(directory, "tables.tex"))
  write.csv(result$estimates, file.path(directory, "estimates.csv"), row.names = FALSE, na = "")
  write.csv(result$seeds, file.path(directory, "seeds.csv"), row.names = FALSE)
  saveRDS(result, file.path(directory, "results.rds"))
  capture.output(result$session, file = file.path(directory, "sessionInfo.txt"))
  if (!any(result$estimates$setting != "null")) return(invisible(tables))
  for (number in 1:12) {
    study <- if (number <= 6) "main" else "broader"
    if (!study %in% result$designs$display$study ||
        !any(result$estimates$family == if (study == "main") "main_curves" else "broader")) next
    png(file.path(directory, sprintf("figure_%02d.png", number)), width = 2100, height = 1600, res = 300)
    tryCatch(plot_figure(result, number), finally = dev.off())
  }
  for (study in c("main", "broader")) {
    family <- if (study == "main") "main_curves" else "broader"
    if (!any(result$estimates$family == family)) next
    name <- if (study == "main") "Figure_1.png" else "Figure_S1.png"
    png(file.path(directory, name), width = 4800, height = 3000, res = 300)
    tryCatch(plot_six_panels(result, study), finally = dev.off())
  }
  invisible(tables)
}

# ---- CHOOSE YOUR SETTINGS ----
# R controls precision of estimated power/size; B controls each reference rank.
# Start with R=20, B=19 to check your setup. Production runs are much slower.
R <- 200L                    # Pilot; increase for final power estimates.
B <- 99L                     # References per observed dataset.
alpha <- .05
master_seed <- 20261118L       # Fixed before looking at power results.
workers <- 4L                 # Four workers plus the coordinating R process.
backend <- "R"                # Optional "cpp" requires Rcpp and a C++ compiler.
designs <- simulation_designs()

# ---- RUN ALL REQUESTED DESIGNS ----
# result <- run_simulations(designs, R=R, B=B, alpha=alpha, workers=workers,
#   backend=backend, bundle_dir=bundle_dir, master_seed=master_seed)
# tables <- make_tables(result)
# tables[["2"]]               # Choose 2,3,4,5,S1,S2,S4,S5,S6.
# plot_figure(result, 1)       # Choose any individual graph from 1 through 12.
# plot_six_panels(result)      # Figure 1.
# plot_six_panels(result, "broader")  # Figure S1.

# ---- RUN ONE FAMILY ONLY ----
# Choose "matched", "module", "moving_average", "main_curves", or "broader".
# module <- run_simulations(designs, families="module", R=R, B=B, alpha=alpha,
#   workers=workers, backend=backend, bundle_dir=bundle_dir, master_seed=master_seed)
# make_tables(module)[["4"]]

# ---- POWER AND TYPE I ERROR FOR THE THREE TABLE DESIGNS ----
# result <- run_simulations(designs,
#   families = c("matched", "module", "moving_average"), null_only = FALSE,
#   R = 1000L, B = 99L, alpha = .05, master_seed = master_seed,
#   workers = 4L, backend = "cpp", bundle_dir = bundle_dir)
# If the matching null results are already in memory, include
# reuse_null = null_results in this call to compute only the alternatives.
# View(result$type1)
# View(result$power)

# ---- OPTIONAL EXPLICIT EXPORT ----
# Running, sourcing and plotting do not save results. Only this command saves:
# export_results(result, "my_results")
# Workers close automatically, including on errors; no cluster remains running.
