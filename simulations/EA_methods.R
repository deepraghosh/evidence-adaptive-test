# Upward test statistics for X ~ N_d(mu, Sigma), A ~ W_d(n, Sigma), X independent of A.
# n always denotes the Wishart degrees of freedom. Source this file to define
# functions. Sourcing it does not simulate data or write results.

method_names <- c("EA", "Z-Max", "Z-TopK", "IHC", "Feng-COM", "QCQ",
                  "ZWS", "MinP", "Quadratic")
.ea_fast <- new.env(parent = baseenv())

two_sided_log_p <- function(z, df = Inf) {
  if (is.infinite(df)) return(log(2) + pnorm(-abs(z), log.p = TRUE))
  log(2) + pt(-abs(z), df = df, log.p = TRUE)
}

ea_statistic <- function(X, A = NULL, n = NULL, Sigma = NULL, inverse = NULL) {
  d <- length(X)
  if (is.null(A)) {
    if (is.null(inverse)) inverse <- chol2inv(chol(Sigma))
    marginal <- X / sqrt(diag(Sigma))
    adjusted <- drop(inverse %*% X) / sqrt(diag(inverse))
    lp <- two_sided_log_p(marginal)
    lq <- two_sided_log_p(adjusted)
  } else {
    if (is.null(inverse)) inverse <- chol2inv(chol(A))
    nu <- n - d + 1
    u <- drop(inverse %*% X)
    q <- sum(X * u)
    denominator <- diag(inverse) * (1 + q) - u^2
    if (any(denominator <= 0)) stop("EA adjusted denominator is nonpositive.")
    lp <- two_sided_log_p(sqrt(n) * X / sqrt(diag(A)), n)
    lq <- two_sided_log_p(sqrt(nu) * u / sqrt(denominator), nu)
  }
  log_simes <- pmin(log(2) + pmin(lp, lq), pmax(lp, lq))
  -log(d) - min(log_simes)
}

minp_statistic <- function(X, A = NULL, n = NULL, Sigma = NULL) {
  if (is.null(A)) {
    lp <- two_sided_log_p(X / sqrt(diag(Sigma)))
  } else {
    lp <- two_sided_log_p(sqrt(n) * X / sqrt(diag(A)), n)
  }
  -min(lp)
}

quadratic_statistic <- function(X, A = NULL, n = NULL, Sigma = NULL, inverse = NULL) {
  if (is.null(A)) {
    if (is.null(inverse)) inverse <- chol2inv(chol(Sigma))
    return(sum(X * (inverse %*% X)))
  }
  if (is.null(inverse)) inverse <- chol2inv(chol(A))
  n * sum(X * (inverse %*% X))
}

gumbel_max_p <- function(maximum, d) {
  centered <- maximum - 2 * log(d) + log(log(d))
  -expm1(-exp(-centered / 2) / sqrt(pi))
}

# Same Lasso objective, tuning, residual variance and unsymmetrized estimate
# as the previous manuscript bundle. The optional C++ version speeds this up.
nodewise_precision <- function(A, n, tolerance = 1e-10, max_iterations = 100000L) {
  d <- nrow(A)
  S <- A / (n + 1)
  penalty <- sqrt(2 * log(d) / (n + 1))
  if (exists("nodewise_lasso_cpp", .ea_fast, inherits = FALSE)) {
    return(.ea_fast$nodewise_lasso_cpp(S, penalty, tolerance, max_iterations)$Gamma)
  }
  Gamma <- matrix(0, d, d)
  for (j in seq_len(d)) {
    others <- setdiff(seq_len(d), j)
    G <- S[others, others, drop = FALSE]
    target <- S[others, j]
    beta <- numeric(d - 1L)
    for (iteration in seq_len(max_iterations)) {
      for (k in seq_along(others)) {
        residual <- target[k] - sum(G[k, ] * beta) + G[k, k] * beta[k]
        beta[k] <- sign(residual) * max(abs(residual) - penalty, 0) / G[k, k]
      }
      gradient <- drop(G %*% beta) - target
      error <- ifelse(beta == 0, pmax(abs(gradient) - penalty, 0),
                      abs(gradient + penalty * sign(beta)))
      if (max(error) <= tolerance) break
    }
    if (max(error) > tolerance) stop("Nodewise Lasso did not converge.")
    variance <- S[j, j] - 2 * sum(beta * target) +
      sum(beta * (G %*% beta)) + penalty * sum(abs(beta))
    if (variance <= 0) stop("Nodewise residual variance is nonpositive.")
    Gamma[j, j] <- 1 / variance
    Gamma[j, others] <- -beta / variance
  }
  Gamma
}

zmax_statistic <- function(X, A = NULL, n = NULL, Sigma = NULL, precision = NULL) {
  if (is.null(precision)) {
    precision <- if (is.null(A)) chol2inv(chol(Sigma)) else nodewise_precision(A, n)
  }
  z <- drop(precision %*% X) / sqrt(diag(precision))
  max(z^2)
}

ztopk_statistic <- function(X, A = NULL, n, Sigma = NULL, precision = NULL) {
  if (is.null(precision)) {
    precision <- if (is.null(A)) chol2inv(chol(Sigma)) else nodewise_precision(A, n)
  }
  z <- drop(precision %*% X) / sqrt(diag(precision))
  squared <- sort(z^2, decreasing = TRUE)
  k <- seq_len(min(40L, length(X), n))
  max(sqrt((1 - k / (n + 1)) / (2 * k)) * (cumsum(squared)[k] - k))
}

higher_criticism <- function(p, lower, upper) {
  p <- sort(p)
  j <- seq_along(p)
  eligible <- p >= lower & p <= upper
  if (!any(eligible)) return(-Inf)
  max(sqrt(length(p)) * (j[eligible] / length(p) - p[eligible]) /
        sqrt(p[eligible] * (1 - p[eligible])))
}

prepare_ihc <- function(Sigma, chol_factor = NULL) {
  d <- nrow(Sigma)
  b <- max(1L, floor(log(d)))
  if (is.null(chol_factor)) chol_factor <- chol(Sigma)
  L <- t(chol_factor)
  if (exists("banded_inverse_cpp", .ea_fast, inherits = FALSE)) {
    U <- .ea_fast$banded_inverse_cpp(L, b)
  } else {
    U <- matrix(0, d, d)
    for (j in seq_len(d)) {
      index <- j:min(d, j + b - 1L)
      U[index, j] <- forwardsolve(L[index, index, drop = FALSE], c(1, rep(0, length(index) - 1L)))
    }
  }
  U <- sweep(U, 2, sqrt(colSums(U^2)), "/")
  list(L = L, U = U, b = b)
}

# The original known-covariance IHC is retained.
ihc_statistic <- function(X, Sigma = NULL, transform = NULL) {
  if (is.null(transform)) transform <- prepare_ihc(Sigma)
  y <- drop(crossprod(transform$U, forwardsolve(transform$L, X)))
  p <- 2 * pnorm(-abs(y))
  higher_criticism(p, 1 / length(X), .5) / sqrt(2 * transform$b - 1)
}

feng_com_statistic <- function(X, A, n) {
  d <- length(X)
  squared <- (n + 1) * X^2 / diag(A)
  correlation <- cov2cor(A)
  variance <- 2 * (sum(correlation^2) - d^2 / n)
  if (variance <= 0) stop("Feng-SUM variance is nonpositive.")
  sum_statistic <- (sum(squared) - n * d / (n - 2)) / sqrt(variance)
  p_sum <- pnorm(sum_statistic, lower.tail = FALSE)
  p_max <- gumbel_max_p(max(squared), d)
  smallest <- min(p_sum, p_max)
  -smallest * (2 - smallest)
}

# Recover Gaussian observations when only X and A are supplied. In simulations,
# we generate the centered rows directly and reuse their moments for QCQ/ZWS.
reconstruct_rows <- function(X, A, n) {
  noise <- matrix(rnorm((n + 1) * length(X)), nrow = n + 1)
  noise <- sweep(noise, 2, colMeans(noise), "-")
  fit <- qr(noise, tol = 0)
  frame <- sweep(qr.Q(fit), 2, sign(diag(qr.R(fit))), "*")
  sweep(frame %*% chol(A), 2, X / sqrt(n + 1), "+")
}

# Newton's identities evaluate elementary symmetric polynomials from power sums.
elementary2 <- function(p1, p2) (p1^2 - p2) / 2
elementary6 <- function(p1, p2, p3, p4, p5, p6) {
  (p1^6 - 15 * p1^4 * p2 + 45 * p1^2 * p2^2 - 15 * p2^3 +
     40 * p1^3 * p3 - 120 * p1 * p2 * p3 + 40 * p3^2 -
     90 * p1^2 * p4 + 90 * p2 * p4 + 144 * p1 * p5 - 120 * p6) / 720
}

prepare_rows <- function(Y) {
  rows <- nrow(Y)
  n <- rows - 1L
  C <- sweep(Y, 2, colMeans(Y), "-")
  scale <- sqrt(max(colSums(C^2) / n))
  if (!is.finite(scale) || scale <= 0) stop("Row scale is nonpositive.")
  C <- C / scale
  if (exists("raw_pair_moments_cpp", .ea_fast, inherits = FALSE)) {
    moments <- .ea_fast$raw_pair_moments_cpp(C)
  } else {
    powers <- lapply(1:6, function(r) C^r)
    pair_sums <- lapply(powers, crossprod)
    moments <- list(A = pair_sums[[1]],
      moment = rbind(rep(rows, ncol(C)), t(sapply(powers, colSums))),
      V2sum = sum(elementary2(pair_sums[[1]], pair_sums[[2]])),
      V6sum = sum(do.call(elementary6, pair_sums)))
  }
  A <- moments$A
  row_norms <- rowSums(C^2)
  moments$n <- n
  moments$scale <- scale
  moments$qcq_v <- colSums(C * row_norms)
  moments$qcq_J0 <- (n^2 * sum(A^2) - rows^2 * sum(row_norms^2) +
                         sum(row_norms)^2) / (n - 1)^2
  moments
}

qcq_statistic <- function(X, A, n, moments = NULL) {
  if (is.null(moments)) moments <- prepare_rows(reconstruct_rows(X, A, n))
  rows <- n + 1
  d <- length(X)
  m <- X / (sqrt(rows) * moments$scale)
  J <- moments$qcq_J0 - 2 * rows * sum(m * moments$qcq_v) / (n - 1)^2 +
    rows * sum(m * (moments$A %*% m)) / (n - 1)
  if (!is.finite(J) || J <= 0) stop("QCQ variance is nonpositive.")
  quadratic <- (rows * n * sum(m^2) - sum(diag(moments$A))) / sqrt(2 * J)
  squared <- n * X^2 / diag(A)
  p <- sort(2 * pnorm(-sqrt(squared)))
  eligible <- p >= log(d) / d & p <= d^(-1 / 4)
  hc <- -Inf
  if (any(eligible)) {
    j <- seq_len(d)[eligible]
    p_scan <- p[eligible]
    hc <- max(sqrt(d) * (j / d - p_scan) / sqrt(p_scan * (1 - p_scan)))
  }
  max_enhancement <- as.numeric(max(squared) > 2 * log(d) + 3 * log(log(d)))
  hc_enhancement <- as.numeric(hc > sqrt(6 * log(log(d)) + log(log(log(d)))))
  quadratic + sqrt(d) / 2 * (max_enhancement + hc_enhancement)
}

zws_statistic <- function(X, A, n, moments = NULL) {
  if (n < 5) stop("ZWS order six requires n+1 >= 6 rows.")
  if (is.null(moments)) moments <- prepare_rows(reconstruct_rows(X, A, n))
  rows <- n + 1
  m <- X / (sqrt(rows) * moments$scale)
  shifted <- matrix(0, 6, length(X))
  for (r in 1:6) {
    for (k in 0:r) {
      shifted[r, ] <- shifted[r, ] + choose(r, k) * m^(r - k) * moments$moment[k + 1L, ]
    }
  }
  if (moments$V2sum <= 0 || moments$V6sum <= 0) stop("ZWS variance is nonpositive.")
  u2 <- sum(elementary2(shifted[1, ], shifted[2, ]))
  u6 <- sum(do.call(elementary6, lapply(1:6, function(r) shifted[r, ])))
  t2 <- rows * u2 / sqrt(2 * choose(rows, 2) * moments$V2sum)
  t6 <- rows^3 * u6 / sqrt(720 * choose(rows, 6) * moments$V6sum)
  smallest <- min(pnorm(t2, lower.tail = FALSE), pnorm(t6, lower.tail = FALSE))
  -smallest * (2 - smallest)
}

# Optional acceleration uses the same formulas. It requires Rcpp and a working
# C++ compiler. Compilation files go to R's temporary directory.
enable_fast <- function(file = "EA_fast.cpp") {
  if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Install Rcpp or use backend = 'R'.")
  Rcpp::sourceCpp(file, env = .ea_fast, cacheDir = tempdir())
  invisible(TRUE)
}
