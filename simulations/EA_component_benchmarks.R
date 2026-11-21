# Analytic EA versus benchmarks built from the same marginal and adjusted p-values.
# Source this file, then call run_component_benchmarks(). No simulation runs on source.

component_designs <- function() {
  equicor <- function(b, rho) (1 - rho) * diag(b) + matrix(rho, b, b)
  block <- function(correlation, d) {
    Sigma <- diag(d)
    index <- seq_len(nrow(correlation))
    Sigma[index, index] <- correlation
    Sigma
  }
  design <- function(Sigma, mu) {
    d <- nrow(Sigma)
    list(Sigma = Sigma, mu = c(mu, rep(0, d - length(mu))), n = 10L * d - 1L)
  }

  rho <- .9
  matched <- design(block(equicor(2, rho), 50),
                    3.5 * c(1, rho / (1 + sqrt(1 - rho^2))))
  module <- design(block(equicor(5, .75), 50), 3.5 * c(1, rep(.5, 4)))

  ma <- diag(50)
  ma[abs(row(ma) - col(ma)) == 1] <- 20/41
  ma_mean <- numeric(50)
  ma_mean[25:27] <- 3.5 * c(1/3, 1, 1/3)

  pairs <- kronecker(diag(25), equicor(2, rho))
  separate <- design(pairs, c(3.2, 3.2, 3.4 * sqrt(1 - rho^2), 0))

  star <- diag(5)
  star[1, 2:5] <- star[2:5, 1] <- .45
  star_peak <- design(block(star, 100), 3.7 * c(1, rep(.32, 4)))

  list(D1 = matched, D2 = module, D3 = design(ma, ma_mean),
       D4 = separate, D5 = star_peak)
}

# P and Q are the marginal and precision-adjusted coordinatewise p-values.
component_pvalues <- function(P, Q) {
  d <- length(P)
  p_min <- min(P)
  q_min <- min(Q)
  raw_min <- min(p_min, q_min)
  ea_pairs <- pmin(2 * pmin(P, Q), pmax(P, Q))
  c(EA = min(1, d * min(ea_pairs)),
    Bonf_P = min(1, d * p_min),
    Bonf_Q = min(1, d * q_min),
    Pooled = min(1, 2 * d * raw_min),
    Raw_min = min(1, d * raw_min),
    Simes_P = min(1, min(d * sort(P) / seq_len(d))))
}

# X has two columns, one under the null and one under the alternative.
known_component_inputs <- function(X, Sigma, Omega) {
  marginal <- X / sqrt(diag(Sigma))
  adjusted <- (Omega %*% X) / sqrt(diag(Omega))
  list(P = 2 * pnorm(-abs(marginal)), Q = 2 * pnorm(-abs(adjusted)))
}

sample_component_inputs <- function(X, A, n) {
  nu <- n - nrow(A) + 1L
  inverse <- chol2inv(chol(A))
  u <- inverse %*% X
  q <- colSums(X * u)
  denominator <- outer(diag(inverse), 1 + q) - u^2
  marginal <- sqrt(n) * X / sqrt(diag(A))
  adjusted <- sqrt(nu) * u / sqrt(denominator)
  list(P = 2 * pt(-abs(marginal), df = n), Q = 2 * pt(-abs(adjusted), df = nu))
}

run_component_benchmarks <- function(R = 1000L, alpha = .05, seed = 20261119L) {
  stopifnot(R >= 2, R == floor(R), alpha > 0, alpha < 1)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  designs <- component_designs()
  methods <- names(component_pvalues(1, 1))
  panels <- c("known", "sample")
  comparisons <- c("Bonf_P", "Bonf_Q", "Pooled", "Simes_P")
  power <- type1 <- paired_gain <- list()

  for (id in names(designs)) {
    g <- designs[[id]]
    d <- nrow(g$Sigma)
    root <- chol(g$Sigma)
    Omega <- chol2inv(root)
    null_X <- crossprod(root, matrix(rnorm(d * R), nrow = d))
    alt_X <- sweep(null_X, 1L, g$mu, "+")
    decisions <- array(FALSE, c(R, length(methods), 2L, 2L),
      dimnames = list(NULL, methods, c("null", "alternative"), panels))

    for (r in seq_len(R)) {
      # Common observations across methods. Under the alternative, only the
      # Gaussian mean changes. Scatter is independent of the Gaussian draws.
      X <- cbind(null_X[, r], alt_X[, r])
      A <- rWishart(1L, df = g$n, Sigma = g$Sigma)[, , 1L]
      inputs <- list(known = known_component_inputs(X, g$Sigma, Omega),
                     sample = sample_component_inputs(X, A, g$n))
      for (panel in panels) {
        p <- inputs[[panel]]$P
        q <- inputs[[panel]]$Q
        decisions[r, , "null", panel] <- component_pvalues(p[, 1], q[, 1]) <= alpha
        decisions[r, , "alternative", panel] <- component_pvalues(p[, 2], q[, 2]) <= alpha
      }
    }

    # Check the pointwise rejection-region orderings on every simulated draw.
    stopifnot(all(!decisions[, "Pooled", , ] | decisions[, "EA", , ]),
              all(!decisions[, "EA", , ] | decisions[, "Raw_min", , ]),
              all(!decisions[, "Bonf_P", , ] | decisions[, "Simes_P", , ]))
    for (panel in panels) {
      key <- paste(id, panel)
      base <- data.frame(design = id, covariance = panel)
      power[[key]] <- cbind(base, as.list(colMeans(decisions[, , "alternative", panel])))
      type1[[key]] <- cbind(base, as.list(colMeans(decisions[, , "null", panel])))

      ea <- as.integer(decisions[, "EA", "alternative", panel])
      paired_gain[[key]] <- do.call(rbind, lapply(comparisons, function(method) {
        difference <- ea - decisions[, method, "alternative", panel]
        gain <- mean(difference)
        mcse <- sd(difference) / sqrt(R)
        data.frame(base, comparator = method, gain = gain, mcse = mcse,
                   lower95 = gain - qnorm(.975) * mcse,
                   upper95 = gain + qnorm(.975) * mcse)
      }))
    }
    message("Completed ", id, " (d = ", d, ", n = ", g$n, ")")
  }

  list(power = do.call(rbind, unname(power)), type1 = do.call(rbind, unname(type1)),
       paired_gain = do.call(rbind, unname(paired_gain)),
       designs = designs, R = R, alpha = alpha, seed = seed)
}

# In RStudio, source this file and run these lines in the Console:
# results <- run_component_benchmarks()
# print(results$power, row.names = FALSE)
# print(results$type1, row.names = FALSE)
