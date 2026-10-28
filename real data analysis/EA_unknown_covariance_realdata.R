# EA_unknown_covariance_realdata.R -- ANALYTIC P-VALUES ONLY
# Source EA_manuscript_simulations.R first, then source this file.
# Nothing runs on sourcing; no numerical results or figures are saved.
# There are no permutations, bootstrap draws, Monte Carlo banks, or workers.
#
# Example:
#   data_dir <- "/path/to/frozen/TCGA/files"
#   tcga <- load_tcga_unknown(data_dir, "stad")
#   groups <- tcga_nodal_groups(tcga)
#   answer <- analytic_unknown_groups(groups$Y1, groups$Y2)
#   answer$table
#   answer$components
#   answer$pooled_summaries  # X, A and n used by the five summary methods
#   answer$qcq_details      # original-row CQ base and enhancement components
#   answer$zws_details      # original-row order-2/order-6 statistics and variances
# CIN versus GS instead:
#   tcga <- load_tcga_unknown(data_dir, "coadread")
#   groups <- tcga_cin_gs_groups(tcga)
#   answer <- analytic_unknown_groups(groups$Y1, groups$Y2)
#   answer$table
# p_value_display and log10_p_value retain extremely small tail probabilities.
# When a positive p-value underflows ordinary R numbers, p_value is NA rather
# than a misleading zero; its display and log10 columns remain available.
#
# All covariance information is estimated from the observations.
# EA, Bonferroni MinP and Quadratic use finite-sample Gaussian calibrations.
# Z-Max uses an ASYMPTOTIC Gumbel approximation, subject to dependence and
# precision-estimation conditions; this is not an exact finite-sample p-value.
# Z-TopK is NA: its retained adaptive scan needs simulation-based calibration.
# IHC is NA: the retained implementation needs known covariance.
# EA, Z-Max, Feng-COM, MinP and Quadratic use only the scaled GROUP-MEAN
# DIFFERENCE X and pooled WITHIN-GROUP scatter A, with n=m1+m2-2.
# The grand pooled mean is not the mean contrast being tested.
# QCQ and ZWS use the ORIGINAL observations Y1 and Y2, preserving both groups.
# No Helmert transformation or reconstructed observations are used here.
# Feng-COM retains the report's summary-based score and asymptotic transform.
# QCQ is a documented two-sample adaptation: Chen--Qin (2010) quadratic base
# plus the retained Qiu enhancement. Its normal tail remains PROVISIONAL:
# the full primary QCQ two-sample specification has not been verified.
# ZWS uses the full two-sample mean U-statistics for orders 2 and 6, with a
# centered original-row variance adaptation and asymptotic analytic tails.
# These data-analysis adaptations differ from the canonical simulation scores.
#
# Original log-scale relative coordinates are retained. Only the Z-Max Lasso
# defaults to pooled coordinate standardization to make its penalty comparable
# across proteins. The returned settings disclose this choice.
#
# STAD was identified by an exploratory screen. These are model-based p-values,
# not selection-adjusted findings; the Gaussian/common-covariance assumptions
# still require assessment. Marginal normality fails for many STAD proteins.

read_tcga_clinical <- function(path) {
  read.delim(path, comment.char = "#", check.names = FALSE,
             stringsAsFactors = FALSE)
}

load_tcga_unknown <- function(data_dir, study = "stad") {
  paths <- file.path(data_dir, paste0(study,
    c("_data_rppa.txt", "_data_clinical_patient.txt")))
  rppa <- read.delim(paths[1L], check.names = FALSE, stringsAsFactors = FALSE)
  values <- as.matrix(rppa[-1L]); storage.mode(values) <- "double"
  # Apply completeness to the whole cohort, without using group labels.
  keep <- apply(values, 1L, function(x) all(is.finite(x)))
  Y <- t(values[keep, , drop = FALSE]); colnames(Y) <- rppa[[1L]][keep]
  patient <- read_tcga_clinical(paths[2L])
  metadata <- patient[match(substr(rownames(Y), 1L, 12L), patient$PATIENT_ID), , drop = FALSE]
  metadata$SAMPLE_ID <- rownames(Y)
  ord <- order(metadata$SAMPLE_ID)
  Y <- Y[ord, , drop = FALSE]; metadata <- metadata[ord, , drop = FALSE]
  if (anyDuplicated(substr(rownames(Y), 1L, 12L)))
    stop("Repeated patients need a prespecified sample-selection rule.")
  list(Y = Y, metadata = metadata, dropped_features = rppa[[1L]][!keep],
       source_md5 = tools::md5sum(paths))
}

tcga_nodal_groups <- function(tcga) {
  lab <- trimws(as.character(tcga$metadata$PATH_N_STAGE))
  g1 <- !is.na(lab) & lab == "N0"
  g2 <- !is.na(lab) & grepl("^N[123]", lab)
  list(Y1 = tcga$Y[g1, , drop = FALSE], Y2 = tcga$Y[g2, , drop = FALSE],
       labels = c("N0", "N+"), group_counts = table(lab[g1 | g2]))
}

tcga_cin_gs_groups <- function(tcga) {
  lab <- trimws(as.character(tcga$metadata$SUBTYPE))
  g1 <- !is.na(lab) & lab == "COAD_CIN"
  g2 <- !is.na(lab) & lab == "COAD_GS"
  if (!any(g1) || !any(g2)) stop("Both COAD_CIN and COAD_GS are required.")
  list(Y1 = tcga$Y[g1, , drop = FALSE], Y2 = tcga$Y[g2, , drop = FALSE],
       labels = c("COAD_CIN", "COAD_GS"), group_counts = table(lab[g1 | g2]))
}

# Sufficient summaries of the original groups under the common-covariance
# Gaussian model. A is scatter, not covariance; the pooled covariance is A/n.
# No individual rows are reconstructed or returned.
canonical_group_data <- function(Y1, Y2, standardize = FALSE) {
  Y1 <- as.matrix(Y1); Y2 <- as.matrix(Y2)
  m1 <- nrow(Y1); m2 <- nrow(Y2); n <- m1 + m2 - 2L; d <- ncol(Y1)
  stopifnot(ncol(Y2) == d, m1 >= 2L, m2 >= 2L, n >= d, n >= 5L, d >= 3L,
            all(is.finite(Y1)), all(is.finite(Y2)))
  X <- sqrt(m1 * m2 / (m1 + m2)) * (colMeans(Y1) - colMeans(Y2))
  C1 <- sweep(Y1, 2L, colMeans(Y1), "-")
  C2 <- sweep(Y2, 2L, colMeans(Y2), "-")
  A <- crossprod(C1) + crossprod(C2)
  scl <- if (standardize) sqrt(diag(A) / n) else rep(1, d)
  X <- X / scl; A <- A / outer(scl, scl)
  list(X = X, A = A, n = n, d = d, m1 = m1, m2 = m2, scale = scl)
}

# The retained Feng-COM score is determined by X,A,n. These expressions are
# algebraically identical to the former raw_scores() route, without making rows.
feng_summary_scores <- function(X, A, n) {
  d <- length(X); diagA <- diag(A)
  if (n <= 2 || any(!is.finite(diagA)) || any(diagA <= 0))
    stop("Undefined Feng-COM summary inputs.")
  correlation <- A / sqrt(outer(diagA, diagA))
  variance <- 2 * (sum(correlation^2) - d*d/n)
  if (!is.finite(variance) || variance <= 0)
    stop("Feng-COM sum variance is non-positive; no p-value returned.")
  w <- (n + 1) * X^2 / diagA
  sum_score <- (sum(w) - n*d/(n - 2)) / sqrt(variance)
  max_score <- max(w)
  log_p <- .log_combine_two(pnorm(sum_score, lower.tail=FALSE, log.p=TRUE),
    .log_gumbel_survival(max_score - 2*log(d) + log(log(d))))
  list(scores=c("Feng-SUM"=sum_score, "Feng-MAX"=max_score,
                "Feng-COM"=-exp(log_p)), log_p_value=log_p)
}

unknown_native_pvalues <- function(Y1, Y2) {
  z <- canonical_group_data(Y1, Y2, standardize = FALSE)
  Ai <- chol2inv(chol(z$A)); u <- drop(Ai %*% z$X); q <- sum(z$X * u)
  lp <- log_two_tail(sqrt(z$n) * z$X / sqrt(diag(z$A)), z$n)
  lq <- log_two_tail(sqrt(z$n-z$d+1L) * u /
    sqrt(diag(Ai) * (1 + q) - u*u), z$n-z$d+1L)
  c(EA = min(1, exp(-ea_logscore(lp, lq))),
    `Min-p Bonferroni` = min(1, z$d * exp(min(lp))),
    `Quadratic exact F` = pf((z$n-z$d+1L)*q/z$d, z$d,
                             z$n-z$d+1L, lower.tail = FALSE))
}

# Stable analytic tails: no resampling and no rounding of a positive p to zero.
.log_gumbel_survival <- function(x) {
  ell <- -x/2 - log(pi)/2
  if (ell < -35) ell else log(-expm1(-exp(ell)))
}
.log_combine_two <- function(lp, lq) {
  lm <- min(lp, lq)
  lm + log(2 - exp(lm))  # log{1-(1-min(p,q))^2}
}
.display_log_p <- function(lp) {
  if (is.na(lp)) return(NA_character_)
  if (lp == -Inf) return("Below finite log range")
  if (lp == 0) return("1")
  exponent <- floor(lp/log(10))
  sprintf("%.6fe%+d", exp(lp-exponent*log(10)), as.integer(exponent))
}

# Original-observation two-sample adaptation of the retained QCQ specification.
# This is NOT claimed to reproduce the unavailable full Qiu--Chen--Qiu paper.
# Verified sources:
# Chen and Qin (2010), Secs 2--3: https://arxiv.org/pdf/1002.4547
# Qiu JSM 2026 presentation, slides 23--24:
# https://ww2.amstat.org/meetings/jsm/2026/onlineprogram/files.cfm?id=138b443a71f9
# The presentation supplies the enhancement in the ideal independent-normal
# model, and states a general extension, without its full two-sample formula.
# We therefore combine the published Chen--Qin two-sample base statistic and
# its leave-out trace estimators with the report's M=2, c=1, delta=1/4 gates.
# Analytic normal calibration remains ASYMPTOTIC / PROVISIONAL for this
# adaptation, especially with correlated coordinates and extreme tail p-values.

.cq_trace_within <- function(Y) {
  m <- nrow(Y)
  G <- tcrossprod(Y)
  # A[i,j] = Y_i' (Y_j - mean excluding i and j).
  A <- ((m - 1) * G + diag(G) - rowSums(G)) / (m - 2)
  diag(A) <- 0
  sum(A * t(A)) / (m * (m - 1))
}

.cq_trace_between <- function(Y1, Y2) {
  m1 <- nrow(Y1)
  m2 <- nrow(Y2)
  H <- tcrossprod(Y1, Y2)
  # L[i,j] = Y1_i' (Y2_j - group-2 mean excluding j).
  # R[i,j] = Y2_j' (Y1_i - group-1 mean excluding i).
  L <- (m2 * H - rowSums(H)) / (m2 - 1)
  R <- sweep(m1 * H, 2L, colSums(H), "-") / (m1 - 1)
  sum(L * R) / (m1 * m2)
}

qcq_original_rows <- function(Y1, Y2, M=2, hc_lower_multiplier=1,
                              hc_upper_exponent=1/4, lambda=1/2) {
  Y1 <- as.matrix(Y1)
  Y2 <- as.matrix(Y2)
  storage.mode(Y1) <- storage.mode(Y2) <- "double"
  m1 <- nrow(Y1)
  m2 <- nrow(Y2)
  d <- ncol(Y1)
  if (m1 < 3 || m2 < 3 || d < 3 || ncol(Y2) != d ||
      any(!is.finite(Y1)) || any(!is.finite(Y2))) {
    stop("QCQ needs two finite matrices with >=3 rows and the same >=3 columns.")
  }
  if (M <= 0 || lambda < 0 || lambda > 1) stop("Invalid QCQ tuning.")
  # One common scale prevents overflow without changing relative coordinates.
  v1 <- apply(Y1, 2L, var)
  v2 <- apply(Y2, 2L, var)
  scl <- sqrt(max(v1, v2))
  if (!is.finite(scl) || scl <= 0) stop("Undefined QCQ scale.")
  Y1 <- Y1 / scl
  Y2 <- Y2 / scl
  v1 <- v1 / scl^2
  v2 <- v2 / scl^2
  difference <- colMeans(Y1) - colMeans(Y2)
  # This is algebraically the ordered-pair U-statistic in CQ2010 Sec 2.
  numerator <- sum(difference^2) - sum(v1) / m1 - sum(v2) / m2
  tr11 <- .cq_trace_within(Y1)
  tr22 <- .cq_trace_within(Y2)
  tr12 <- .cq_trace_between(Y1, Y2)
  variance <- 2 * tr11 / (m1 * (m1 - 1)) +
    2 * tr22 / (m2 * (m2 - 1)) + 4 * tr12 / (m1 * m2)
  if (!is.finite(variance) || variance <= 0) {
    stop("QCQ Chen--Qin leave-out variance is non-positive; no p-value returned.")
  }
  T2 <- numerator / sqrt(variance)
  # Preserve the report's common-covariance studentization for enhancement.
  # Only the CQ base variance retains its native group-specific trace terms.
  pooled_variance <- ((m1 - 1) * v1 + (m2 - 1) * v2) / (m1 + m2 - 2)
  marginal_variance <- (1 / m1 + 1 / m2) * pooled_variance
  if (any(!is.finite(marginal_variance)) || any(marginal_variance <= 0)) {
    stop("Undefined QCQ marginal variance.")
  }
  z <- difference / sqrt(marginal_variance)
  p <- sort(2 * pnorm(-abs(z)))
  take <- p >= hc_lower_multiplier * log(d) / d &
    p <= d^(-hc_upper_exponent)
  H <- if (any(take)) {
    pp <- p[take]
    max(sqrt(d) * (seq_len(d)[take] / d - pp) / sqrt(pp * (1 - pp)))
  } else -Inf
  max_cut <- 2 * log(d) + (2 * M - 1) * log(log(d))
  hc_cut_squared <- 2 * (M + 1) * log(log(d)) + log(log(log(d)))
  if (!is.finite(hc_cut_squared) || hc_cut_squared < 0) {
    stop("QCQ HC threshold is undefined for this dimension.")
  }
  Emax <- as.integer(max(z^2) > max_cut)
  Ehc <- as.integer(H > sqrt(hc_cut_squared))
  statistic <- T2 + sqrt(d) * (lambda * Emax + (1 - lambda) * Ehc)
  list(statistic=statistic,
       log_p_value=pnorm(statistic, lower.tail=FALSE, log.p=TRUE),
       components=c(CQ=T2, numerator=numerator, variance=variance,
                    trace11=tr11, trace22=tr22, trace12=tr12,
                    HC=H, Emax=Emax, Ehc=Ehc),
       implementation="Original-row two-sample QCQ adaptation; CQ2010 base, retained enhancement",
       calibration="Provisional asymptotic normal tail; full Qiu two-sample specification unverified")
}

# Native two-sample mean U-statistics using the original two groups.
# Numerator: ZWS arXiv 2303.08197 Section 2.2, kernel x_j-y_j;
# equivalently He et al. (2021), equation (4.2), Theorem 4.3.
# The centered analytic variance below extends ZWS Section 4.3 to the two
# independent groups. It is an explicitly disclosed adaptation, not a claim
# that ZWS's two-sample Section 4.5 supplied this exact studentizer.
# For true population covariances, the exact null variance is
# sum_s choose(q,s)^2/[choose(m1,s)*choose(m2,q-s)] *
#       sum_jk Sigma1[j,k]^s * Sigma2[j,k]^(q-s).
# Sample-centered covariance products estimate these terms. Centering introduces
# finite-sample bias; using the exact combinatorial weights does NOT make the
# estimated variance or its normal-tail p-value exact. Mean-special-case
# asymptotic variance: He et al. (2021), Theorem 4.3 and Corollary 4.1,
# https://www.wuchong.org/assets/pdf/he2021asymptotically.pdf (pp. 25--27).

# e_k for each coordinate, divided by choose(m,k). Degree zero is one.
.zws_coordinate_products <- function(Y, max_order) {
  m <- nrow(Y); d <- ncol(Y)
  E <- matrix(0, max_order + 1L, d); E[1L, ] <- 1
  for (i in seq_len(m)) {
    for (k in seq.int(min(i, max_order), 1L))
      E[k + 1L, ] <- E[k + 1L, ] + Y[i, ] * E[k, ]
  }
  E / choose(m, 0:max_order)
}

# For each k,j,l: average product of k distinct centered rows' C_ij C_il.
# Newton's identities calculate e_k from crossprod(C^r), avoiding subsets.
.zws_covariance_products <- function(C, max_order) {
  m <- nrow(C); d <- ncol(C)
  powers <- vector("list", max_order)
  C_power <- C
  for (r in seq_len(max_order)) {
    powers[[r]] <- crossprod(C_power)
    if (r < max_order) C_power <- C_power * C
  }
  E <- vector("list", max_order + 1L)
  E[[1L]] <- matrix(1, d, d)
  for (k in seq_len(max_order)) {
    value <- matrix(0, d, d)
    for (r in seq_len(k))
      value <- value + (-1)^(r - 1L) * E[[k - r + 1L]] * powers[[r]]
    E[[k + 1L]] <- value / k
  }
  lapply(0:max_order, function(k) E[[k + 1L]] / choose(m, k))
}

zws_original_two_sample <- function(Y1, Y2, orders = c(2L, 6L)) {
  Y1 <- as.matrix(Y1); Y2 <- as.matrix(Y2)
  storage.mode(Y1) <- storage.mode(Y2) <- "double"
  m1 <- nrow(Y1); m2 <- nrow(Y2); d <- ncol(Y1)
  if (ncol(Y2) != d || any(!is.finite(Y1)) || any(!is.finite(Y2)))
    stop("The two groups need finite observations on the same coordinates.")
  if (!length(orders) || any(!is.finite(orders)) || anyDuplicated(orders) ||
      any(orders %% 2L != 0L) || any(orders < 2L) ||
      max(orders) > min(m1, m2))
    stop("Orders must be positive even integers no larger than either group.")

  C1 <- sweep(Y1, 2L, colMeans(Y1), "-")
  C2 <- sweep(Y2, 2L, colMeans(Y2), "-")
  common_scale <- sqrt(max((colSums(C1 * C1) + colSums(C2 * C2)) /
                             (m1 + m2 - 2L)))
  if (!is.finite(common_scale) || common_scale <= 0)
    stop("A positive pooled variance is required.")
  # One common coordinate shift and one scalar scale improve numerics and do
  # not change any studentized statistic; there is no observation mixing.
  shift <- (m1 * colMeans(Y1) + m2 * colMeans(Y2)) / (m1 + m2)
  P1 <- .zws_coordinate_products(sweep(Y1, 2L, shift, "-") / common_scale,
                                 max(orders))
  P2 <- .zws_coordinate_products(sweep(Y2, 2L, shift, "-") / common_scale,
                                 max(orders))
  V1 <- .zws_covariance_products(C1 / common_scale, max(orders))
  V2 <- .zws_covariance_products(C2 / common_scale, max(orders))
  U <- variance <- numeric(length(orders))
  for (a in seq_along(orders)) {
    q <- orders[a]
    for (s in 0:q) {
      U[a] <- U[a] + choose(q, s) * (-1)^(q - s) *
        sum(P1[s + 1L, ] * P2[q - s + 1L, ])
      variance[a] <- variance[a] + choose(q, s)^2 /
        (choose(m1, s) * choose(m2, q - s)) *
        sum(V1[[s + 1L]] * V2[[q - s + 1L]])
    }
  }
  if (any(!is.finite(variance)) || any(variance <= 0))
    stop("The centered analytic ZWS variance is undefined; no truncation used.")
  T <- U / sqrt(variance)
  log_p <- pnorm(T, lower.tail = FALSE, log.p = TRUE)
  lm <- min(log_p)
  combined_log_p <- if (lm < -35) lm + log(length(orders)) else
    log(-expm1(length(orders) * log1p(-exp(lm))))
  list(statistics = setNames(T, paste0("ZWS-", orders)),
       log_p = setNames(log_p, paste0("ZWS-", orders)),
       combined_log_p = combined_log_p,
       # These two quantities share the numerical common-scale convention.
       numerator = setNames(U, orders), variance = setNames(variance, orders),
       common_scale = common_scale,
       calibration = paste("Asymptotic normal plus independent-min combination;",
         "full native two-sample mean numerator, centered analytic variance adaptation"))
}

analytic_unknown_groups <- function(Y1, Y2, z_standardize = TRUE,
                                     backend = "base") {
  if (!exists("unknown_scores", mode = "function"))
    stop("Source EA_manuscript_simulations.R first.")
  if (nrow(Y1) < 6L || nrow(Y2) < 6L)
    stop("The full comparison needs at least six original observations per group for ZWS-6.")
  enable_acceleration(backend)
  # Summary methods receive X,A,n; row methods receive the ORIGINAL two groups.
  pooled <- canonical_group_data(Y1, Y2, standardize = FALSE)
  zdata <- pooled
  if (z_standardize) {
    scl <- sqrt(diag(pooled$A) / pooled$n)
    zdata$X <- pooled$X / scl
    zdata$A <- pooled$A / outer(scl, scl)
  }
  scores <- unknown_scores(zdata$X, unknown_prepare(zdata$A, zdata$n))
  native <- unknown_native_pvalues(Y1, Y2)
  feng <- feng_summary_scores(pooled$X, pooled$A, pooled$n)
  qcq <- qcq_original_rows(Y1, Y2)
  zws <- zws_original_two_sample(Y1, Y2)
  gumbel_x <- unname(scores["Z-Max"]) - 2*log(pooled$d) + log(log(pooled$d))
  # Analytic tails only, evaluated in log space to avoid underflow.
  logp <- c(min(0, -scores["EA"]), .log_gumbel_survival(gumbel_x), NA_real_,
    feng$log_p_value, qcq$log_p_value, zws$combined_log_p,
    min(0, log(pooled$d)-scores["Min-p"]),
    pf((pooled$n-pooled$d+1L)*scores["Quadratic"]/(pooled$d*pooled$n),
       pooled$d, pooled$n-pooled$d+1L, lower.tail=FALSE, log.p=TRUE), NA_real_)
  p <- exp(logp)
  # NA identifies a value that cannot be represented as an ordinary R number;
  # the finite log10 and display columns give its value without calling it 0.
  p[!is.na(p) & p == 0] <- NA_real_
  table <- data.frame(
    method = c("EA", "Z-Max", "Z-TopK", "Feng-COM (summaries)",
               "QCQ (two-sample adaptation)", "ZWS-{2,6} (two-sample adaptation)",
               "MinP (Bonferroni)", "Quadratic (Hotelling)", "IHC"),
    p_value = unname(p),
    p_value_display = unname(vapply(logp, .display_log_p, "")),
    log10_p_value = unname(logp/log(10)),
    calibration = c("Finite-sample Gaussian", "Asymptotic Gumbel",
      "Unavailable: simulated calibration required", "Asymptotic combination",
      "Provisional two-sample normal approximation", "Asymptotic two-sample combination",
      "Finite-sample Gaussian", "Finite-sample Gaussian",
      "Unavailable: known covariance required"),
    input = c(rep("Group-mean contrast and pooled scatter", 4L),
      rep("Original observations in their two groups", 2L),
      rep("Group-mean contrast and pooled scatter", 2L), "Unavailable"),
    stringsAsFactors = FALSE)
  stopifnot(all(is.na(p) | (is.finite(p) & p >= 0 & p <= 1)))
  list(table = table,
       components = c(scores, feng$scores, QCQ=qcq$statistic, zws$statistics,
                      "ZWS-{2,6}"=-exp(zws$combined_log_p)),
       native_pvalues = native, qcq_details = qcq, zws_details = zws,
       pooled_summaries = pooled,
       m1 = nrow(Y1), m2 = nrow(Y2), d = pooled$d, n = pooled$n,
       z_standardize = z_standardize,
       moment_coordinates = "Original coordinates and original groups; no reconstruction",
       calibration = "Analytic only; exact Gaussian versus asymptotic explicitly labeled",
       results_saved = FALSE)
}

# Compatibility name for the earlier example. This now performs analytic
# calculations only: B, seed, and workers are deliberately not arguments.
compare_unknown_groups <- analytic_unknown_groups

gaussian_fit_summary <- function(Y1, Y2) {
  summarize <- function(Y) {
    p <- apply(Y, 2L, function(x) shapiro.test(x)$p.value)
    C <- scale(Y, TRUE, FALSE); variance <- colMeans(C*C)
    skew <- colMeans(C^3)/variance^(3/2)
    kurtosis <- colMeans(C^4)/variance^2 - 3
    c(subjects = nrow(Y), dimension = ncol(Y),
      marginal_Shapiro_p_below_05 = sum(p < .05),
      BH_Shapiro_p_below_05 = sum(p.adjust(p, "BH") < .05),
      median_absolute_skewness = median(abs(skew)),
      maximum_absolute_skewness = max(abs(skew)),
      maximum_excess_kurtosis = max(kurtosis))
  }
  list(groups = rbind(group1 = summarize(Y1), group2 = summarize(Y2)),
    variance_ratio_range = range(apply(Y1, 2L, var)/apply(Y2, 2L, var)),
    interpretation = "Nonsignificant marginal tests cannot establish joint Gaussianity.")
}
