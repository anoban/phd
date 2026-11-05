#!/usr/bin/env Rscript
## ── Exactness checks of the analytic integration ────────────────────
##
##   Rscript R/validate.R              # all checks (compiles Stan)
##   Rscript R/validate.R --no-stan    # R-only checks, seconds
##
## The Stan programs replace the brms likelihood with derived
## likelihoods in which the species effects are integrated out.  These
## checks compare the derivations with direct computations:
##
##   1  compound-symmetry formula vs a dense 3n_g multivariate normal
##   2  the pooled-W model-block form vs the sum of per-species terms
##   3  the O(1) rank-one LOO factor vs O(n) recomputation
##   4  Felsenstein pruning vs the dense A (x) Sigma_p construction
##   5  the LOO predictive (LOO-PIT) vs the ratio-of-marginals log_lik
##   6  Stan's log_lik and loo_pit_draw vs the R reference, at real draws
##
## Every check should pass at machine precision (< 1e-8).

## The package root is the parent of this file's directory under
## Rscript; when sourced, the working directory should already be it.
.script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
BASE_DIR <- if (length(.script))
  dirname(dirname(normalizePath(sub("^--file=", "", .script[1])))) else getwd()
setwd(BASE_DIR)

suppressPackageStartupMessages(library(ape))
source("R/build_data.R")

RUN_STAN <- !("--no-stan" %in% commandArgs(TRUE))
set.seed(20260830)

bundle <- build_anoban_data("data/stanforum.csv", "data/FRED4_1301.tre",
                            log_response = TRUE)
sd0 <- bundle$standata
D <- sd0$D; G <- sd0$G; N <- sd0$N
Y <- sd0$Y; gidx <- sd0$gidx; ng <- sd0$ng

pass <- function(dev, tol = 1e-8) if (dev < tol) "PASS" else "FAIL"
ok <- TRUE
report <- function(label, dev, tol = 1e-8) {
  ok <<- ok && dev < tol
  cat(sprintf("%-52s %9.3g  (%s)\n", label, dev, pass(dev, tol)))
}

## ── R reference implementations ─────────────────────────────────────
dmvn <- function(x, mu, S) {
  U <- chol(S); z <- backsolve(U, x - mu, transpose = TRUE)
  -0.5 * (length(x) * log(2 * pi) + 2 * sum(log(diag(U))) + sum(z^2))
}
## dense: y_g stacked observation-major, cov = I_n (x) Sr + J_n (x) Sb
block_dense <- function(Yg, mu, Sr, Sb) {
  n <- nrow(Yg)
  V <- kronecker(diag(n), Sr) + kronecker(matrix(1, n, n), Sb)
  dmvn(as.vector(t(Yg)), rep(mu, n), V)
}
block_fast <- function(Yg, mu, Sr, Sb) {
  n <- nrow(Yg); d <- ncol(Yg); yb <- colMeans(Yg)
  lp <- dmvn(yb, mu, Sr / n + Sb)
  if (n > 1) {
    U <- chol(Sr); Z <- backsolve(U, t(Yg) - yb, transpose = TRUE)
    lp <- lp - 0.5 * (2 * (n - 1) * sum(log(diag(U))) + sum(Z^2) +
                      d * (n - 1) * log(2 * pi))
  }
  lp - 0.5 * d * log(n)
}
## the O(1) rank-one form actually used in Stan
loo_fast <- function(Yg, mu, Sr, Sb) {
  n <- nrow(Yg); d <- ncol(Yg); yb <- colMeans(Yg)
  A_n <- dmvn(yb, mu, Sr / n + Sb)
  if (n == 1) return(A_n)
  n1 <- n - 1; U <- chol(Sr); ldSr <- 2 * sum(log(diag(U)))
  vapply(seq_len(n), function(i) {
    di <- Yg[i, ] - yb; yb1 <- yb - di / n1
    q <- sum(backsolve(U, di, transpose = TRUE)^2)
    A_n - dmvn(yb1, mu, Sr / n1 + Sb) -
      0.5 * (ldSr + (n / n1) * q + d * log(2 * pi)) -
      0.5 * d * (log(n) - log(n1))
  }, numeric(1))
}
rand_spd <- function(d, s = 1) { L <- matrix(rnorm(d * d), d, d); s * (L %*% t(L) / d + diag(d)) }

## Felsenstein pruning, upward pass (mirrors anoban_pruning.stan)
prune_up <- function(ta, ybar, ng, mu, Sr, St, Sp) {
  nn <- ta$n_node; D <- ncol(ybar); ld2pi <- D * log(2 * pi)
  Pacc <- array(0, c(nn, D, D)); hacc <- matrix(0, nn, D); gacc <- numeric(nn)
  for (v in ta$post_order) {
    if (v == ta$root) break
    if (ta$is_leaf[v] == 1L) {
      g <- ta$species_of[v]
      V <- Sr / ng[g] + St; m <- ybar[g, ] - mu[g, ]; K <- 0
    } else {
      P <- Pacc[v, , ]; h <- hacc[v, ]
      V <- solve(P); m <- as.vector(V %*% h)
      K <- gacc[v] + 0.5 * sum(m * h) + 0.5 * (ld2pi + log(det(V)))
    }
    W <- V + ta$edge_len[v] * Sp
    Pw <- solve(W); hw <- as.vector(Pw %*% m)
    p <- ta$parent[v]
    Pacc[p, , ] <- Pacc[p, , ] + Pw
    hacc[p, ] <- hacc[p, ] + hw
    gacc[p] <- gacc[p] + K - 0.5 * (ld2pi + log(det(W)) + sum(m * hw))
  }
  gacc[ta$root]
}
dense_between <- function(A, ybar, ng, mu, Sr, St, Sp) {
  G <- nrow(ybar); D <- ncol(ybar)
  V <- kronecker(A, Sp)
  for (g in seq_len(G)) {
    ii <- ((g - 1) * D + 1):(g * D)
    V[ii, ii] <- V[ii, ii] + Sr / ng[g] + St
  }
  dmvn(as.vector(t(ybar)), as.vector(t(mu)), V)
}

Sr <- rand_spd(D, 0.7); St <- rand_spd(D, 1.3); Sp <- rand_spd(D, 0.4)
mu_all <- matrix(rnorm(G * D), G, D)

cat("\n── data completeness ──\n")
## Every data variable either Stan model declares must be present in the
## list build_anoban_data() returns (plus cor_t / cor_p, which
## fit_anoban() sets per model).  Cheap, and it catches the case where a
## new Stan data variable is added but not wired into the builder.
declared_data_vars <- function(path) {
  txt <- readLines(path)
  i <- grep("^data \\{", txt)[1]
  j <- grep("^\\}", txt); j <- j[j > i][1]
  d <- trimws(sub("//.*", "", txt[(i + 1):(j - 1)]))
  d <- d[grepl(";", d)]
  sub("^.*[[:space:]]", "", sub(";.*$", "", d))
}
decl <- unique(unlist(lapply(
  file.path("stan", c("anoban_marginal.stan", "anoban_pruning.stan")),
  declared_data_vars)))
have <- c(names(sd0), "cor_t", "cor_p")
missing_vars <- setdiff(decl, have)
cat(sprintf("%-52s %9d  (%s)\n", "0  Stan data variables missing from standata",
            length(missing_vars), if (!length(missing_vars)) "PASS" else "FAIL"))
if (length(missing_vars)) {
  ok <- FALSE
  cat("   missing: ", paste(missing_vars, collapse = ", "), "\n")
}

cat("\n── R-level checks ──\n")

## 1: compound-symmetry formula vs a dense 3n_g MVN, every species
d1 <- max(vapply(seq_len(G), function(g) {
  Yg <- Y[gidx == g, , drop = FALSE]
  abs(block_fast(Yg, mu_all[g, ], Sr, St) - block_dense(Yg, mu_all[g, ], Sr, St))
}, numeric(1)))
report("1  formula vs dense MVN, all species", d1)

## 2: pooled-W form (what the model block uses) vs per-species sum
sum_fast <- sum(vapply(seq_len(G), function(g)
  block_fast(Y[gidx == g, , drop = FALSE], mu_all[g, ], Sr, St), numeric(1)))
between <- sum(vapply(seq_len(G), function(g)
  dmvn(sd0$ybar[g, ], mu_all[g, ], Sr / ng[g] + St), numeric(1)))
Ur <- chol(Sr); Zw <- backsolve(Ur, t(chol(sd0$W)), transpose = TRUE)
within <- -0.5 * (2 * sd0$N_minus_G * sum(log(diag(Ur))) + sum(Zw^2) +
                  D * sd0$N_minus_G * log(2 * pi)) - sd0$sum_half_D_log_ng
report("2  pooled-W vs per-species sum", abs(sum_fast - (between + within)), 1e-6)

## 3: O(1) rank-one LOO factor vs O(n) recomputation
gsel <- c(which(ng == 1)[1], which(ng == 2)[1], which(ng == 5)[1], which.max(ng))
d3 <- max(vapply(gsel, function(g) {
  Yg <- Y[gidx == g, , drop = FALSE]; mu <- mu_all[g, ]
  ref <- vapply(seq_len(nrow(Yg)), function(i)
    block_fast(Yg, mu, Sr, St) -
      if (nrow(Yg) == 1) 0 else block_fast(Yg[-i, , drop = FALSE], mu, Sr, St),
    numeric(1))
  max(abs(ref - loo_fast(Yg, mu, Sr, St)))
}, numeric(1)))
report("3  O(1) LOO factor vs O(n) recomputation", d3)

## 4: pruning vs the dense A (x) Sigma_p construction, full tree
tree <- read.tree("data/FRED4_1301.tre")
A <- vcv.phylo(tree, corr = TRUE)[bundle$species, bundle$species]
ta <- build_tree_arrays(tree, bundle$species)
up <- prune_up(ta, sd0$ybar, ng, mu_all, Sr, St, Sp)
dn <- dense_between(A, sd0$ybar, ng, mu_all, Sr, St, Sp)
report("4  pruning vs dense A (x) Sigma_p", abs(up - dn) / abs(dn), 1e-10)

## 5: the LOO predictive density (used for LOO-PIT) comes from
## conditioning on the species effect w given the other observations of
## the species, while log_lik is a ratio of species marginals.  These
## are two independent derivations of the same density, so they must
## agree.
loo_pred <- function(Yg, mu, Sr, Sb, want = c("logdens", "pit")) {
  want <- match.arg(want)
  n <- nrow(Yg); d <- ncol(Yg); yb <- colMeans(Yg)
  n1 <- if (n > 1) n - 1 else 0
  Sr_inv <- solve(Sr)
  Sw <- if (n1 > 0) solve(solve(Sb) + n1 * Sr_inv) else Sb
  Spred <- Sr + Sw; sd_pred <- sqrt(diag(Spred))
  out <- lapply(seq_len(n), function(i) {
    mu_pred <- mu
    if (n1 > 0) {
      yb1 <- yb - (Yg[i, ] - yb) / n1
      mu_pred <- as.vector(mu + Sw %*% (n1 * (Sr_inv %*% (yb1 - mu))))
    }
    if (want == "logdens") dmvn(Yg[i, ], mu_pred, Spred)
    else {
      ## same clamp as loo_pit_vec() in anoban_functions.stanfunctions
      eps <- .Machine$double.eps
      pmin(pmax(stats::pnorm((Yg[i, ] - mu_pred) / sd_pred), eps), 1 - eps)
    }
  })
  if (want == "logdens") unlist(out) else do.call(rbind, out)
}
d5 <- max(vapply(gsel, function(g) {
  Yg <- Y[gidx == g, , drop = FALSE]; mu <- mu_all[g, ]
  max(abs(loo_fast(Yg, mu, Sr, St) - loo_pred(Yg, mu, Sr, St, "logdens")))
}, numeric(1)))
report("5  LOO predictive vs ratio-of-marginals log_lik", d5)

if (!RUN_STAN) {
  cat(sprintf("\n%s (R-level checks only)\n", if (ok) "ALL PASS" else "FAILURES"))
  quit(status = if (ok) 0 else 1)
}

## 6: Stan's log_lik and loo_pit_draw vs the R reference, at real draws
cat("\n── Stan checks (short sampling runs) ──\n")
source("R/fit_models.R")
bundle$standata$compute_pit <- 1L
for (m in c("M4", "M3")) {
  fit <- fit_anoban(m, bundle, stan_dir = "stan", chains = 1,
                    iter_warmup = 150, iter_sampling = 10, sig_figs = 18,
                    refresh = 0, show_messages = FALSE)
  pick <- function(v) as.numeric(fit$draws(variables = v, format = "draws_matrix")[1, ])
  sigma <- pick("sigma")
  LSr <- diag(sigma) %*% matrix(pick("Lrescor"), D, D); Sr_d <- LSr %*% t(LSr)
  LSt <- matrix(pick("LSigma_t"), D, D); St_d <- LSt %*% t(LSt)
  b_d <- matrix(pick("b"), sd0$Kc, D)
  mu_d <- matrix(pick("Intercept"), G, D, byrow = TRUE) + sd0$Xc %*% b_d

  if (m == "M3") {
    ## In the pruning models the species mean and the between-species
    ## covariance include the downward message, which Stan does not
    ## output.  So this check compares the marginal likelihood instead.
    LSp <- matrix(pick("LSigma_p"), D, D); Sp_d <- LSp %*% t(LSp)
    ref <- prune_up(ta, sd0$ybar, ng, mu_d, Sr_d, St_d, Sp_d)
    report(sprintf("6  %s Stan lp_between vs R pruning", m),
           abs(pick("lp_between") - ref) / abs(ref), 1e-10)
  } else {
    idx_of <- function(g) which(gidx == g)
    ll_R <- numeric(N); pit_R <- matrix(NA_real_, N, D)
    for (g in seq_len(G)) {
      Yg <- Y[gidx == g, , drop = FALSE]; ii <- idx_of(g)
      ll_R[ii] <- loo_fast(Yg, mu_d[g, ], Sr_d, St_d)
      pit_R[ii, ] <- loo_pred(Yg, mu_d[g, ], Sr_d, St_d, "pit")
    }
    report(sprintf("6a %s Stan log_lik vs R reference", m),
           max(abs(pick("log_lik") - ll_R)))
    pit_stan <- matrix(pick("loo_pit_draw"), N, D)   # column-major (i, d)
    report(sprintf("6b %s Stan loo_pit_draw vs R reference", m),
           max(abs(pit_stan - pit_R)))
  }
}

cat(sprintf("\n%s\n", if (ok) "ALL CHECKS PASS" else "FAILURES: do not trust the output"))
quit(status = if (ok) 0 else 1)
