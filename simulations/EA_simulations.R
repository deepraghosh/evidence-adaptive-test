# Design definitions and simulations. All results stay in R memory.
# Source EA_methods.R first. Each observation receives B conditional references.

block_covariance <- function(sizes, correlations) {
  Sigma <- diag(sum(sizes))
  start <- 1L
  for (j in seq_along(sizes)) {
    index <- seq.int(start, length.out = sizes[j])
    Sigma[index, index] <- correlations[j]
    diag(Sigma)[index] <- 1
    start <- start + sizes[j]
  }
  attr(Sigma, "block_sizes") <- sizes
  attr(Sigma, "block_correlations") <- correlations
  Sigma
}

matched_mean <- function(d, rho, a = 3.5, offset = 0) {
  t_star <- rho / (1 + sqrt(1 - rho^2))
  c(a, a * (t_star + offset), rep(0, d - 2))
}

simulation_designs <- function() {
  groups <- list()
  display <- list()
  add <- function(id, family, Sigma, n, means, null_label = "") {
    groups[[id]] <<- list(id = id, family = family, Sigma = Sigma, n = n,
                         means = means, null_label = null_label)
  }
  point <- function(group, setting, study, sweep, x) {
    display[[length(display) + 1L]] <<- data.frame(group, setting, study, sweep, x)
  }

  # 1. MATCHED PAIRS. Tables 2 and 3 retain all M1--M9, B1--B4 and C1--C3.
  fixed <- data.frame(setting = c("M1", "M2", "M3", "M4", "M5", "M8", "M9"),
    d = c(20, 20, 50, 50, 50, 50, 100), n = c(99, 99, 249, 99, 249, 499, 499),
    rho = c(.9, .95, .9, .95, .95, .95, .95), a = c(rep(3.5, 6), 3.75))
  for (j in seq_len(nrow(fixed))) {
    g <- fixed[j, ]
    means <- setNames(list(matched_mean(g$d, g$rho, g$a)), g$setting)
    if (g$setting == "M5") {
      means$M6 <- matched_mean(g$d, g$rho, g$a, -.05)
      means$M7 <- matched_mean(g$d, g$rho, g$a, .05)
    }
    label <- if (g$setting == "M5") "M5/M6/M7" else g$setting
    add(g$setting, "matched", block_covariance(c(2, g$d - 2), c(g$rho, 0)),
        g$n, means, label)
  }
  for (d in c(50, 100)) {
    labels <- if (d == 50) c("B1", "B2") else c("B3", "B4")
    means <- setNames(lapply(c(.6, .8), function(a) c(3.5, 3.5, a, -a, rep(0, d - 4))), labels)
    add(labels[1], "matched", block_covariance(c(2, 2, d - 4), c(.95, .9, 0)),
        5 * d - 1, means, paste(labels, collapse = "/"))
  }
  pair_covariance <- block_covariance(c(2, 48), c(.95, 0))
  add("C1", "matched", pair_covariance, 249, list(C1 = c(3.5, 3.5, rep(0, 48))), "C1")
  add("C2", "matched", diag(50), 249, list(C2 = c(3.5, rep(0, 49))), "C2")
  add("C3", "matched", pair_covariance, 249, list(C3 = c(1.25, rep(0, 49))), "C3")

  # 2. CORRELATED MODULE. Table 4 keeps d=50, n=499, rho=.75 and a=3,3.5,4.
  strengths <- c(3, 3.5, 4)
  add("module", "module", block_covariance(c(5, 45), c(.75, 0)), 499,
      setNames(lapply(strengths, function(a) c(a, rep(a / 2, 4), rep(0, 45))), strengths),
      "Correlated module")

  # 3. MOVING AVERAGE. Table 5 keeps adjacent covariance 20/41 and peak 25:27.
  Sigma <- diag(50)
  Sigma[abs(row(Sigma) - col(Sigma)) == 1] <- 20 / 41
  means <- lapply(strengths, function(a) {
    mu <- numeric(50)
    mu[c(25, 26, 27)] <- c(a / 3, a, a / 3)
    mu
  })
  add("moving_average", "moving_average", Sigma, 499, setNames(means, strengths),
      "Moving average peak")

  # Both correlation figures now use .1,.3,.5,.7,.9, as requested.
  # S1/S2 retain their separately registered F1--F9 null groups.
  main_correlations <- c(.1, .3, .5, .7, .9)
  null_correlations <- c(.50, .70, .80, .90, .95)
  main_dimensions <- c(20, 50, 100, 150, 200)
  main_strengths <- c(2.5, 3, 3.5, 4, 4.5)
  for (j in seq_along(null_correlations)) {
    rho <- null_correlations[j]
    id <- paste0("F", j)
    amplitudes <- if (rho == .95) main_strengths else 3.5
    means <- setNames(lapply(amplitudes, function(a) matched_mean(50, rho, a)), amplitudes)
    add(id, "main_curves", block_covariance(c(2, 48), c(rho, 0)), 249, means, id)
  }
  for (rho in c(.1, .3)) {
    add(paste0("rho", 100 * rho), "main_curves", block_covariance(c(2, 48), c(rho, 0)),
        249, list("3.5" = matched_mean(50, rho)), "")
  }
  correlation_groups <- c("rho10", "rho30", "F1", "F2", "F4")
  for (j in seq_along(main_correlations)) {
    point(correlation_groups[j], "3.5", "main", "rho", main_correlations[j])
  }
  for (j in seq_along(main_dimensions)) {
    d <- main_dimensions[j]
    id <- if (d == 50) "F5" else paste0("F", match(d, c(20, 100, 150, 200)) + 5)
    if (d != 50) add(id, "main_curves", block_covariance(c(2, d - 2), c(.95, 0)),
                     5 * d - 1, list("3.5" = matched_mean(d, .95)), id)
    point(id, "3.5", "main", "dimension", d)
  }
  for (a in main_strengths) point("F5", as.character(a), "main", "amplitude", a)

  # BROADER STUDY. Tables S4/S5 and Figure S1 preserve the original sparse designs.
  broader_correlations <- c(.1, .3, .5, .7, .9)
  broader_dimensions <- c(50, 100, 150, 200, 250)
  broader_strengths <- c(1, 1.5, 2, 2.5, 3)
  for (j in seq_along(broader_correlations)) {
    rho <- broader_correlations[j]
    means <- setNames(list(c(2, rep(0, 49))), paste0("O", j))
    if (rho == .9) {
      for (k in seq_along(c(1, 1.5, 2.5, 3))) {
        a <- c(1, 1.5, 2.5, 3)[k]
        means[[paste0("O", k + 5)]] <- c(a, rep(0, 49))
      }
    }
    id <- paste0("original_rho", j)
    add(id, "broader", block_covariance(c(5, 45), c(rho, .1)), 249, means)
    point(id, paste0("O", j), "broader", "rho", rho)
  }
  amplitude_labels <- c("O6", "O7", "O5", "O8", "O9")
  for (j in seq_along(broader_strengths)) {
    point("original_rho5", amplitude_labels[j], "broader", "amplitude", broader_strengths[j])
  }
  for (j in seq_along(broader_dimensions)) {
    d <- broader_dimensions[j]
    id <- paste0("original_d", d)
    label <- paste0("O", j + 9)
    add(id, "broader", block_covariance(rep(10, d / 10), c(.7, rep(.1, d / 10 - 1))),
        3 * d - 1, setNames(list(c(3, rep(.5, 9), rep(0, d - 10))), label))
    point(id, label, "broader", "dimension", d)
  }
  list(groups = groups, display = do.call(rbind, display))
}

# Equicorrelation blocks have an inexpensive symmetric square root.
draw_gaussian_rows <- function(rows, Sigma, root) {
  noise <- matrix(rnorm(rows * nrow(Sigma)), rows, nrow(Sigma))
  sizes <- attr(Sigma, "block_sizes")
  if (is.null(sizes)) return(noise %*% root)
  correlations <- attr(Sigma, "block_correlations")
  start <- 1L
  for (j in seq_along(sizes)) {
    index <- seq.int(start, length.out = sizes[j])
    rho <- correlations[j]
    block <- noise[, index, drop = FALSE]
    noise[, index] <- sqrt(1 - rho) * block +
      (sqrt(1 + (sizes[j] - 1) * rho) - sqrt(1 - rho)) * rowMeans(block)
    start <- start + sizes[j]
  }
  noise
}

# Evaluate every complete score. Internal p-transforms in EA/Feng/ZWS are
# part of their statistics; the final rejection rule is conditional Monte Carlo.
score_dataset <- function(X, A, context, moments = NULL) {
  n <- context$n
  scores <- matrix(NA_real_, 2, length(method_names),
                   dimnames = list(c("known", "unknown"), method_names))
  safe <- function(expression) {
    value <- tryCatch(expression, error = function(e) NA_real_)
    if (length(value) != 1L || is.na(value)) return(NA_real_)
    unname(value)  # -Inf, including an empty HC scan, is a valid score.
  }
  chol_A <- tryCatch(chol(A), error = function(e) NULL)
  inverse <- if (is.null(chol_A)) NA_real_ else chol2inv(chol_A)
  precision <- tryCatch(nodewise_precision(A, n), error = function(e) NA_real_)
  if (is.null(moments)) moments <- tryCatch(prepare_rows(reconstruct_rows(X, A, n)),
                                          error = function(e) NA_real_)
  shared <- c("Feng-COM" = safe(feng_com_statistic(X, A, n)),
              QCQ = safe(qcq_statistic(X, A, n, moments)),
              ZWS = safe(zws_statistic(X, A, n, moments)))
  scores[2, names(shared)] <- shared
  if (!is.null(context$Sigma)) scores[1, names(shared)] <- shared
  if (!is.null(context$Sigma)) {
    scores[1, "EA"] <- safe(ea_statistic(X, Sigma = context$Sigma, inverse = context$Omega))
    scores[1, "Z-Max"] <- safe(zmax_statistic(X, precision = context$Omega))
    scores[1, "Z-TopK"] <- safe(ztopk_statistic(X, n = n, precision = context$Omega))
    scores[1, "IHC"] <- safe(ihc_statistic(X, transform = context$ihc))
    scores[1, "MinP"] <- safe(minp_statistic(X, Sigma = context$Sigma))
    scores[1, "Quadratic"] <- safe(quadratic_statistic(X, inverse = context$Omega,
                                                      Sigma = context$Sigma))
  }
  scores[2, "EA"] <- safe(ea_statistic(X, A, n, inverse = inverse))
  scores[2, "Z-Max"] <- safe(zmax_statistic(X, precision = precision))
  scores[2, "Z-TopK"] <- safe(ztopk_statistic(X, n = n, precision = precision))
  scores[2, "MinP"] <- safe(minp_statistic(X, A, n))
  scores[2, "Quadratic"] <- safe(quadratic_statistic(X, A, n, inverse = inverse))
  scores
}

# For one OBSERVED dataset, generate its own B references at fixed a=A+XX'.
# All methods share these references. Every A-dependent fit is recomputed.
conditional_test <- function(X, A, n, Sigma = NULL, B = 99L, alpha = .05,
                             context = NULL,
                             moments = NULL, keep_reference = FALSE) {
  d <- length(X)
  stopifnot(n >= d, n == floor(n), B >= 1, B == floor(B), alpha > 0, alpha < 1,
            is.matrix(A), all(dim(A) == d),
            is.null(Sigma) || (is.matrix(Sigma) && all(dim(Sigma) == d)))
  if (is.null(context)) context <- list(n = n, Sigma = Sigma,
    Omega = if (is.null(Sigma)) NULL else chol2inv(chol(Sigma)),
    ihc = if (is.null(Sigma)) NULL else prepare_ihc(Sigma))
  observed <- score_dataset(X, A, context, moments)
  a <- A + tcrossprod(X)
  root <- t(chol(a))
  reference <- array(NA_real_, c(B, 2, length(method_names)),
                     dimnames = list(NULL, c("known", "unknown"), method_names))
  for (b in seq_len(B)) {
    Z <- rnorm(d)
    V <- rchisq(1, n - d + 1)
    X_b <- drop(root %*% Z) / sqrt(sum(Z^2) + V)
    A_b <- a - tcrossprod(X_b)
    reference[b, , ] <- score_dataset(X_b, A_b, context)
  }
  available <- matrix(TRUE, 2, length(method_names), dimnames = dimnames(observed))
  available[1, ] <- !is.null(context$Sigma)
  available[2, "IHC"] <- FALSE
  rejected <- failed <- matrix(FALSE, 2, length(method_names), dimnames = dimnames(observed))
  critical <- matrix(NA_real_, 2, length(method_names), dimnames = dimnames(observed))
  k <- B + 1L - floor(alpha * (B + 1L))
  for (panel in 1:2) for (method in seq_along(method_names)) {
    if (!available[panel, method]) next
    null_scores <- reference[, panel, method]
    failed[panel, method] <- is.na(observed[panel, method]) || anyNA(null_scores)
    # If any evaluation is undefined, that method does not reject this dataset.
    if (failed[panel, method]) next
    critical[panel, method] <- if (k > B) Inf else sort(null_scores, partial = k)[k]
    rejected[panel, method] <- observed[panel, method] > critical[panel, method]
  }
  output <- list(statistic = observed, critical = critical, rejected = rejected,
                 failed = failed, available = available)
  if (keep_reference) output$reference <- reference
  output
}

# Separate L'Ecuyer streams for groups, substreams for observed data/reference
# banks. Changing B, R, worker count or chunk size does not change earlier draws.
simulation_seeds <- function(designs, master_seed) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(master_seed)
  state <- .Random.seed
  seeds <- list()
  for (g in designs$groups) {
    seeds[[g$id]] <- state
    state <- parallel::nextRNGStream(state)
  }
  seeds
}

simulate_chunk <- function(job, B, alpha, null_only = FALSE, skip_null = FALSE) {
  g <- job$group
  d <- nrow(g$Sigma)
  context <- list(Sigma = g$Sigma, Omega = chol2inv(chol(g$Sigma)), n = g$n,
                  ihc = prepare_ihc(g$Sigma))
  root <- chol(g$Sigma)
  means <- c(list(null = numeric(d)), g$means)
  if (null_only) means <- means["null"]
  if (skip_null) means <- g$means
  count <- failures <- array(0L, c(length(means), 2, length(method_names)),
    dimnames = list(names(means), c("known", "unknown"), method_names))
  state <- job$seed
  for (iteration in seq_len(job$R)) {
    # Reserve the seed positions for the reused null observation and its bank.
    if (skip_null) state <- parallel::nextRNGSubStream(parallel::nextRNGSubStream(state))
    for (setting in names(means)) {
      # Observed X is independent of the centered Gaussian scatter A.
      assign(".Random.seed", state, .GlobalEnv)
      Y <- draw_gaussian_rows(g$n + 1, g$Sigma, root)
      C <- sweep(Y, 2, colMeans(Y), "-")
      A <- crossprod(C)
      X <- drop(draw_gaussian_rows(1, g$Sigma, root)) + means[[setting]]
      moments <- tryCatch(prepare_rows(C), error = function(e) NA_real_)
      state <- parallel::nextRNGSubStream(state)
      assign(".Random.seed", state, .GlobalEnv)
      test <- conditional_test(X, A, g$n, g$Sigma, B, alpha,
                                context = context, moments = moments)
      count[setting, , ] <- count[setting, , ] + test$rejected
      failures[setting, , ] <- failures[setting, , ] + test$failed
      state <- parallel::nextRNGSubStream(state)
    }
    # Leave the unused alternatives' seed positions in place so each null
    # replicate matches the corresponding one in a full run.
    if (null_only) {
      for (step in seq_len(2L * length(g$means))) state <- parallel::nextRNGSubStream(state)
    }
  }
  list(group = g$id, count = count, R = job$R, failures = failures)
}

# Pointwise binomial Wilson intervals from the outer rejection counts.
wilson_intervals <- function(result, level = .95) {
  stopifnot(length(level) == 1L, is.finite(level), level > 0, level < 1)
  rows <- subset(result$estimates, !is.na(rejections))
  R <- rows$R
  p <- rows$rejections / R
  z <- qnorm((1 + level) / 2)
  denominator <- 1 + z^2 / R
  center <- (p + z^2 / (2 * R)) / denominator
  half_width <- z * sqrt(p * (1 - p) / R + z^2 / (4 * R^2)) / denominator

  under_null <- rows$setting == "null"
  groups <- result$designs$groups
  labels <- vapply(groups, function(g) {
    if (is.null(g$null_label)) "" else g$null_label
  }, character(1))
  names(labels) <- vapply(groups, function(g) g$id, character(1))
  null_labels <- unname(labels[rows$group[under_null]])
  unlabelled <- is.na(null_labels) | null_labels == ""
  null_labels[unlabelled] <- rows$group[under_null][unlabelled]
  rows$setting[under_null] <- null_labels
  rows$metric <- ifelse(under_null, "type1", "power")
  rows$lower <- pmax(0, center - half_width)
  rows$upper <- pmin(1, center + half_width)
  rows[, c("metric", "family", "setting", "panel", "method", "estimate", "lower", "upper")]
}

run_simulations <- function(designs = simulation_designs(), families = NULL,
                            R = 200L, B = 99L, alpha = .05, workers = 4L,
                            backend = "R", bundle_dir = getwd(),
                            master_seed = 20261118L,
                            chunk_size = 5L, progress = TRUE, null_only = FALSE,
                            reuse_null = NULL) {
  stopifnot(R >= 1, R == floor(R), B >= 1, B == floor(B),
            workers >= 1, workers <= 4, workers == floor(workers),
            chunk_size >= 1, chunk_size == floor(chunk_size),
            alpha > 0, alpha < 1, backend %in% c("R", "cpp"),
            is.logical(null_only), length(null_only) == 1L, !is.na(null_only))
  if (floor(alpha * (B + 1)) == 0) warning("B is too small to reject at this alpha.")
  old_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) .Random.seed else NULL
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, .GlobalEnv)
  }, add = TRUE)
  seeds <- simulation_seeds(designs, master_seed)
  thread_variables <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                        "VECLIB_MAXIMUM_THREADS")
  previous_threads <- Sys.getenv(thread_variables, unset = NA_character_)
  do.call(Sys.setenv, as.list(setNames(rep("1", length(thread_variables)), thread_variables)))
  on.exit({
    Sys.unsetenv(thread_variables[is.na(previous_threads)])
    saved <- previous_threads[!is.na(previous_threads)]
    if (length(saved)) do.call(Sys.setenv, as.list(saved))
  }, add = TRUE)
  groups <- designs$groups
  if (!is.null(families)) groups <- Filter(function(g) g$family %in% families, groups)
  if (!length(groups)) stop("No requested design groups.")
  seed_table <- do.call(rbind, lapply(groups, function(g) {
    state <- as.data.frame(t(seeds[[g$id]]))
    names(state) <- paste0("state", 1:7)
    data.frame(group = g$id, master_seed = master_seed,
               substreams_per_replicate = 2L * (length(g$means) + 1L), state)
  }))
  skip_null <- !is.null(reuse_null)
  if (skip_null) {
    if (null_only) stop("Use null_only = FALSE when reusing completed null results.")
    settings <- list(R = R, B = B, alpha = alpha, master_seed = master_seed, backend = backend)
    if (!isTRUE(all.equal(settings, reuse_null[names(settings)], tolerance = 0))) {
      stop("Reused null results must match R, B, alpha, master_seed and backend.")
    }
    same_designs <- vapply(groups, function(g) identical(g, reuse_null$designs$groups[[g$id]]), logical(1))
    saved_seeds <- reuse_null$seeds[match(seed_table$group, reuse_null$seeds$group), , drop = FALSE]
    if (!all(same_designs) || !isTRUE(all.equal(seed_table, saved_seeds,
                                               check.attributes = FALSE, tolerance = 0))) {
      stop("Reused null results must match the selected designs and their seed streams.")
    }
    saved_null <- subset(reuse_null$estimates, setting == "null" & group %in% seed_table$group)
    for (id in seed_table$group) {
      by_group <- saved_null[saved_null$group == id, ]
      methods_by_panel <- table(factor(by_group$method, levels = method_names),
                               factor(by_group$panel, levels = c("known", "unknown")))
      if (nrow(by_group) != 2L * length(method_names) || any(methods_by_panel != 1L)) {
        stop("Reused null results are incomplete for ", id, ".")
      }
    }
  }
  jobs <- list()
  for (g in groups) {
    if (skip_null && !length(g$means)) next
    state <- seeds[[g$id]]
    for (first in seq.int(1L, R, by = chunk_size)) {
      size <- min(chunk_size, R - first + 1L)
      jobs[[length(jobs) + 1L]] <- list(group = g, seed = state, R = size)
      for (j in seq_len(2L * (length(g$means) + 1L) * size)) {
        state <- parallel::nextRNGSubStream(state)
      }
    }
  }
  if (backend == "cpp") enable_fast(file.path(bundle_dir, "EA_fast.cpp")) else {
    rm(list = ls(.ea_fast), envir = .ea_fast)
  }
  pool <- NULL
  use_fork <- workers > 1 && .Platform$OS.type == "unix" && Sys.getenv("RSTUDIO") != "1"
  if (workers > 1 && !use_fork) {
    pool <- parallel::makePSOCKcluster(workers, outfile = "")
    on.exit(parallel::stopCluster(pool), add = TRUE)
    parallel::clusterCall(pool, function(path, backend) {
      source(file.path(path, "EA_methods.R"), local = .GlobalEnv)
      source(file.path(path, "EA_simulations.R"), local = .GlobalEnv)
      if (backend == "cpp") enable_fast(file.path(path, "EA_fast.cpp"))
      TRUE
    }, normalizePath(bundle_dir), backend)
  }
  started <- proc.time()[[3]]
  completed <- vector("list", length(jobs))
  batch_starts <- if (length(jobs)) seq.int(1L, length(jobs), by = workers) else integer()
  for (first in batch_starts) {
    index <- seq.int(first, min(length(jobs), first + workers - 1L))
    arguments <- list(B = B, alpha = alpha, null_only = null_only, skip_null = skip_null)
    batch <- if (use_fork) {
      do.call(parallel::mclapply, c(list(X = jobs[index], FUN = simulate_chunk,
        mc.cores = workers, mc.set.seed = FALSE), arguments))
    } else if (is.null(pool)) {
      do.call(lapply, c(list(X = jobs[index], FUN = simulate_chunk), arguments))
    } else do.call(parallel::parLapply, c(list(cl = pool, X = jobs[index], fun = simulate_chunk), arguments))
    if (any(vapply(batch, inherits, logical(1), "try-error"))) stop("A simulation worker failed.")
    completed[index] <- batch
    if (progress) message("Completed ", max(index), "/", length(jobs), " chunks; elapsed ",
                          round((proc.time()[[3]] - started) / 60, 1), " min")
  }
  rows <- list()
  for (g in groups) {
    pieces <- completed[vapply(completed, function(x) x$group == g$id, logical(1))]
    if (!length(pieces)) {
      rows[[length(rows) + 1L]] <- saved_null[saved_null$group == g$id, ]
      next
    }
    counts <- Reduce(`+`, lapply(pieces, `[[`, "count"))
    failures <- Reduce(`+`, lapply(pieces, `[[`, "failures"))
    z <- as.data.frame.table(counts, responseName = "rejections", stringsAsFactors = FALSE)
    names(z)[1:3] <- c("setting", "panel", "method")
    z$group <- g$id
    z$family <- g$family
    z$R <- R
    z$B <- B
    z$alpha <- alpha
    z$estimate <- z$rejections / R
    unavailable <- z$method == "IHC" & z$panel == "unknown"
    z$estimate[unavailable] <- NA_real_
    z$rejections[unavailable] <- NA_integer_
    z$master_seed <- master_seed
    z$n <- g$n
    z$d <- nrow(g$Sigma)
    z$failed_calibrations <- as.vector(failures)
    if (skip_null) {
      z <- rbind(saved_null[saved_null$group == g$id, ], z)
      z <- z[order(match(z$method, method_names), match(z$panel, c("known", "unknown")),
                   match(z$setting, c("null", names(g$means)))), ]
    }
    rows[[length(rows) + 1L]] <- z
  }
  estimates <- do.call(rbind, rows)
  row.names(estimates) <- NULL
  result <- list(estimates = estimates, designs = designs, seeds = seed_table,
                 R = R, B = B, alpha = alpha, master_seed = master_seed,
                 backend = backend, null_only = null_only,
                 implementation_md5 = setNames(unname(tools::md5sum(file.path(bundle_dir,
                   c("EA_methods.R", "EA_simulations.R", "EA_run.R", "EA_fast.cpp")))),
                   c("EA_methods.R", "EA_simulations.R", "EA_run.R", "EA_fast.cpp")),
                 elapsed_seconds = proc.time()[[3]] - started, session = sessionInfo())
  if (skip_null) result$reused_null <- reuse_null[c("implementation_md5", "session")]
  intervals <- wilson_intervals(result)
  result$type1 <- subset(intervals, metric == "type1", select = -metric)
  result$power <- subset(intervals, metric == "power", select = -metric)
  result
}
