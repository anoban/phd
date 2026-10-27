## cv_core.R --- cross-validation functions copied from loobench
##
## PROVENANCE.  Each definition below is the unchanged source text of
## the same definition in loobench.  The comment above
## each definition names its source file and lines.  The comments that
## precede a definition in loobench are not copied; see METHODS.md of
## loobench for the full description of each method.
##
## Do not edit these definitions here.  To update them, copy them again
## from loobench (the extraction keeps the text exactly).
##
## The functions, with the sections of loobench METHODS.md:
##
##   * the CV bridge, plain and mean-MM components: .lb_bb_core and the
##     .lb_bridge_* functions (section 6).  The bridge trust rule
##     .lb_bridge_suspect (section 10.5) flags a fold when the plain
##     value is not finite, OR its MCSE >= 0.3, OR the cross-density gap
##     > 50 nats.
##   * bridge sampling of one log marginal likelihood against a fitted
##     normal proposal: .lb_normal_bridge_logml (section 8.11).  It
##     replaces the bridgesampling package.
##   * PSIS for one fold or a matrix of folds: psis_elpd_fold, pareto_k
##     (sections 3 and 10.1).  The Pareto k uses r_eff = 1 and the right
##     tail.
##   * brute force with its truncated-GPD MCSE and Jensen bias:
##     .lb_bf_vector_stats, truncgpd_mcse (sections 7.1 and 10.2).
##   * the compile arguments for cmdstanr < 1.0 and >= 1.0:
##     .lb_compile_args.
##
## An adapter `ad` is a list of two functions on a fit bundle (a list
## with element `fit`, a CmdStanMCMC object): `unconstrain(bundle)`
## returns the S x P matrix of unconstrained draws, and
## `log_density(bundle, upars)` the log density with Jacobian at the
## rows of `upars`.  R/logo.R builds it.

## loobench R/utils.R, lines 14-14
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

## loobench R/est_bridge.R, lines 126-129
LB_BRIDGE_DEFAULTS <- list(tol = 1e-7, maxiter = 1000L,
                           gap_wide = 50, gap_disjoint = 500,
                           iter_drop_t_max = 3L,
                           mcse_gate = 0.3, ij_mcse = FALSE)

## loobench R/est_bridge.R, lines 146-147
.lb_bridge_config <- function(config)
  utils::modifyList(LB_BRIDGE_DEFAULTS, config %||% list())

## loobench R/est_bridge.R, lines 178-185
.lb_require_fit <- function(bundle, what, who, spec_id) {
  if (!is.null(bundle$fit)) return(invisible(TRUE))
  stop(sprintf(paste0("%s(%s): the %s bundle has no fit object. This method ",
                      "evaluates log_prob() at draws the sampler never saw, ",
                      "so a draws-only bundle (fixture, restored checkpoint) cannot ",
                      "run it. Refit the cell, or run the test on the cluster."),
               who, spec_id, what))
}

## loobench R/est_bridge.R, lines 191-195
.lb_bridge_err <- function(msg, gap = NA_real_, tier_gap = NA_character_)
  list(elpd = NA_real_, mcse = NA_real_, gap = gap, tier_gap = tier_gap,
       log_r = NA_real_, n_iter = 0L, converged = FALSE, error = msg,
       ess_numi = NA_real_, ess_deni = NA_real_, influence = NULL,
       n_drop_numi = NA_integer_, n_drop_deni = NA_integer_)

## loobench R/est_bridge.R, lines 237-353
.lb_bb_core <- function(l1, l2, cfg = list(), return_influence = FALSE,
                        n_chains1 = 1L, iid2 = FALSE) {
  cfg <- .lb_bridge_config(cfg)
  l1 <- as.numeric(l1)
  l2 <- as.numeric(l2)
  s_loo <- length(l1)
  s_full <- length(l2)
  if (!s_loo || !s_full) return(.lb_bridge_err("empty log-ratio vector"))
  if (any(!is.finite(l1)) || anyNA(l2) || any(l2 == Inf))
    return(.lb_bridge_err("non-finite log-ratio"))

  ## Cross-density overlap. A diagnostic only: the log-scale iteration
  ## is defined at any gap. `.lb_bridge_suspect()` flags `gap > gap_wide`.
  gap <- abs(stats::median(l1) - stats::median(l2))
  tier_gap <- if (gap > cfg$gap_disjoint) "disjoint"
              else if (gap > cfg$gap_wide) "wide" else "ok"

  n_ch1 <- as.integer(n_chains1)
  ess1_of <- function(x) .lb_bb_ess(x, n_ch1)
  ess2_of <- function(x) if (isTRUE(iid2)) length(x) else .lb_bb_ess(x, 1L)
  lstar <- stats::median(l1)
  ess_loo <- ess1_of(l1)
  ess_full <- ess2_of(l2)
  s1 <- ess_loo / (ess_loo + ess_full)
  s2 <- ess_full / (ess_loo + ess_full)

  ## The iteration on `lr = log r` (relative to `lstar`), from r = 1.
  log_s1 <- log(s1)
  log_s2 <- log(s2)
  l1c <- l1 - lstar
  l2c <- l2 - lstar
  lae <- function(a, b) matrixStats::rowLogSumExps(cbind(a, b))
  lmean <- function(x) matrixStats::logSumExp(x) - log(length(x))
  lr <- 0
  conv <- FALSE
  step <- NA_real_
  log_numi <- NULL
  log_deni <- NULL
  iter_used <- 0L
  n_drop_numi <- 0L
  n_drop_deni <- 0L
  for (t in seq_len(as.integer(cfg$maxiter))) {
    iter_used <- t
    log_numi <- l2c - lae(log_s1 + l2c, log_s2 + lr)
    log_deni <- -lae(log_s1 + l1c, log_s2 + lr)
    ok_numi <- !is.na(log_numi)
    ok_deni <- !is.na(log_deni)
    if (t <= as.integer(cfg$iter_drop_t_max)) {
      n_drop_numi <- n_drop_numi + sum(!ok_numi)
      n_drop_deni <- n_drop_deni + sum(!ok_deni)
      if (!any(ok_numi))
        return(.lb_bridge_err(sprintf("all numi NA at iter %d (n=%d)",
                                      t, length(log_numi)), gap, tier_gap))
      if (!any(ok_deni))
        return(.lb_bridge_err(sprintf("all deni NA at iter %d (n=%d)",
                                      t, length(log_deni)), gap, tier_gap))
      lr_new <- lmean(log_numi[ok_numi]) - lmean(log_deni[ok_deni])
    } else {
      if (any(!ok_numi) || any(!ok_deni))
        return(.lb_bridge_err(sprintf("NA numi/deni at iter %d", t),
                              gap, tier_gap))
      lr_new <- lmean(log_numi) - lmean(log_deni)
    }
    if (!is.finite(lr_new))
      return(.lb_bridge_err(sprintf("non-finite log r = %g at iter %d",
                                    lr_new, t), gap, tier_gap))
    step <- abs(lr_new - lr)
    lr <- lr_new
    if (step < cfg$tol) {
      conv <- TRUE
      break
    }
  }
  ## Without convergence the value is not an estimate. This occurs when
  ## the two draw sets do not overlap (a disjoint gap): the iteration
  ## then moves `log r` by a constant step at each iteration.
  if (!conv)
    return(.lb_bridge_err(sprintf(
      "not converged in %d iterations (last |d log r| = %.3g, gap %.1f nats)",
      as.integer(cfg$maxiter), step, gap), gap, tier_gap))
  log_r <- lr + lstar

  ## MCSE: delta method on r = mean(numi) / mean(deni), with the ESS
  ## correction applied to both sides, then the log-normal map:
  ##   var_r / r^2 = (mean(numi) / mean(deni) / r)^2 *
  ##                 (rel_var(numi) + rel_var(deni)) / n
  ## The first factor is 1 at the fixed point.
  numi_s <- exp(log_numi - max(log_numi))
  deni_s <- exp(log_deni - max(log_deni))
  ess_numi <- ess2_of(numi_s)
  ess_deni <- ess1_of(deni_s)
  rel_numi <- stats::var(numi_s) * length(numi_s) / ess_numi / mean(numi_s)^2
  rel_deni <- stats::var(deni_s) * length(deni_s) / ess_deni / mean(deni_s)^2
  log_ab <- lmean(log_numi) - lmean(log_deni)
  var_logml <- log1p(exp(2 * (log_ab - lr)) * (rel_numi + rel_deni) /
                       length(numi_s))
  if (!is.finite(var_logml) || var_logml < 0) var_logml <- NA_real_

  ## The influence terms need `numi`, `deni` and `r` on the linear scale.
  ## A side that underflows there (a very wide gap) gives no influence.
  influence <- NULL
  if (isTRUE(return_influence)) {
    numi <- exp(log_numi)
    deni <- exp(log_deni)
    if (mean(numi) > 0 && mean(deni) > 0)
      influence <- .lb_bridge_influence(numi, deni, l2, lstar, s1, s2,
                                        exp(lr), s_full, s_loo)
  }

  list(elpd = -log_r, mcse = sqrt(var_logml), gap = gap, tier_gap = tier_gap,
       log_r = log_r, n_iter = iter_used,
       converged = conv,
       error = NA_character_, ess_numi = ess_numi, ess_deni = ess_deni,
       influence = influence,
       n_drop_numi = as.integer(n_drop_numi),
       n_drop_deni = as.integer(n_drop_deni))
}

## loobench R/est_bridge.R, lines 359-365
.lb_bb_ess <- function(x, n_chains = 1L) {
  n <- length(x)
  arg <- if (n_chains > 1L && n %% n_chains == 0L) matrix(x, ncol = n_chains)
         else x
  v <- tryCatch(as.numeric(posterior::ess_mean(arg)), error = function(e) n)
  if (length(v) != 1L || !is.finite(v) || v <= 0) n else v
}

## loobench R/est_bridge.R, lines 373-385
.lb_bridge_influence <- function(numi, deni, l2, lstar, s1, s2, r,
                                 s_full, s_loo) {
  tryCatch({
    if (!all(is.finite(numi)) || !all(is.finite(deni))) return(NULL)
    a_bar <- mean(numi)
    b_bar <- mean(deni)
    g <- 1 / (s1 * exp(l2 - lstar) + s2 * r)
    d <- -s2 * mean(numi * g) - b_bar + r * s2 * mean(deni * deni)
    if (!all(is.finite(g)) || !is.finite(d) || d >= 0) return(NULL)
    list(L2 = (numi - a_bar) / (r * d * s_full),
         L1 = -(deni - b_bar) / (d * s_loo))
  }, error = function(e) NULL)
}

## loobench R/est_bridge.R, lines 449-460
.lb_mm_mean_transform <- function(upars_full, upars_loo) {
  if (ncol(upars_full) != ncol(upars_loo))
    stop(sprintf(paste0(".lb_mm_mean_transform: parameter dim mismatch ",
                        "(full=%d, fold=%d). The fold refit must keep every ",
                        "parameter of the full fit (the weights = 0 trick)."),
                 ncol(upars_full), ncol(upars_loo)))
  delta <- colMeans(upars_loo) - colMeans(upars_full)
  list(forward = function(u) sweep(u, 2L, delta, "+"),
       backward = function(u) sweep(u, 2L, delta, "-"),
       delta = delta, log_jac = 0,
       delta_norm = sqrt(sum(delta^2)), delta_max = max(abs(delta)))
}

## loobench R/est_bridge.R, lines 484-499
.lb_bridge_qs <- function(ad, full, fold_refit, who = "est_cvbridge",
                          spec_id = "?") {
  .lb_require_fit(full, "full-data", who, spec_id)
  .lb_require_fit(fold_refit, "fold", who, spec_id)
  upars_full <- ad$unconstrain(full)
  upars_loo <- ad$unconstrain(fold_refit)
  if (ncol(upars_full) != ncol(upars_loo))
    return(list(error = sprintf(
      "parameter dim mismatch: full=%d, fold=%d (use the weights = 0 trick)",
      ncol(upars_full), ncol(upars_loo))))
  lp_full <- function(u) as.numeric(ad$log_density(full, u))
  lp_loo <- function(u) as.numeric(ad$log_density(fold_refit, u))
  list(upars_full = upars_full, upars_loo = upars_loo,
       lp_full = lp_full, lp_loo = lp_loo,
       q22 = lp_full(upars_full), q11 = lp_loo(upars_loo), error = NULL)
}

## loobench R/est_bridge.R, lines 503-509
.lb_bridge_arm_plain <- function(qs, cfg, return_influence = FALSE) {
  if (!is.null(qs$error)) return(.lb_bridge_err(qs$error))
  q21 <- qs$lp_loo(qs$upars_full)
  q12 <- qs$lp_full(qs$upars_loo)
  .lb_bb_core(qs$q11 - q12, q21 - qs$q22, cfg,
              return_influence = return_influence)
}

## loobench R/est_bridge.R, lines 515-527
.lb_bridge_arm_mm <- function(qs, cfg, return_influence = FALSE) {
  if (!is.null(qs$error)) return(.lb_bridge_err(qs$error))
  tr <- tryCatch(.lb_mm_mean_transform(qs$upars_full, qs$upars_loo),
                 error = function(e) NULL)
  if (is.null(tr)) return(.lb_bridge_err("MM transform (mean) failed"))
  q21_mm <- qs$lp_loo(tr$forward(qs$upars_full))
  q12_mm <- qs$lp_full(tr$backward(qs$upars_loo))
  r <- .lb_bb_core(qs$q11 - q12_mm, q21_mm - qs$q22, cfg,
                   return_influence = return_influence)
  r$delta_norm <- tr$delta_norm
  r$delta_max <- tr$delta_max
  r
}

## loobench R/est_bridge.R, lines 581-587
.lb_bridge_suspect <- function(bb, cfg) {
  e <- suppressWarnings(as.numeric(bb$elpd %||% NA_real_))
  m <- suppressWarnings(as.numeric(bb$mcse %||% NA_real_))
  g <- suppressWarnings(as.numeric(bb$gap %||% NA_real_))
  !isTRUE(is.finite(e)) || !isTRUE(is.finite(m)) ||
    isTRUE(m >= cfg$mcse_gate) || isTRUE(g > cfg$gap_wide)
}

## loobench R/est_bridge.R, lines 614-627
.lb_bf_truncgpd <- function(log_lik, n_chains = 1L, y_max = Inf,
                            y_max_alt = NULL) {
  na <- list(mcse = NA_real_, bias = NA_real_, k = NA_real_,
             mcse_hybrid = NA_real_, bias_hybrid = NA_real_)
  d <- tryCatch(truncgpd_mcse(log_lik, y_max = y_max, y_max_alt = y_max_alt,
                              n_chains = n_chains),
                error = function(e) NULL)
  if (is.null(d)) return(na)
  list(mcse = as.numeric(d$mcse_truncgpd),
       bias = as.numeric(d$bias_jensen_truncgpd),
       k = as.numeric(d$k_truncgpd),
       mcse_hybrid = as.numeric(d$mcse_truncgpd_hybrid),
       bias_hybrid = as.numeric(d$bias_jensen_truncgpd_hybrid))
}

## loobench R/est_bridge.R, lines 719-732
.lb_bf_vector_stats <- function(ll, n_chains = 1L, y_max = Inf,
                                y_max_full = NA_real_) {
  ll <- as.numeric(ll)
  k_bf <- if (all(is.finite(ll))) pareto_k(exp(ll - max(ll))) else NA_real_
  tg <- .lb_bf_truncgpd(ll, n_chains, y_max = y_max,
                        y_max_alt = .lb_hybrid_bound(y_max, y_max_full))
  list(elpd = matrixStats::logSumExp(ll) - log(length(ll)),
       mcse = .lb_mcse_log_mean_exp(ll, n_chains),
       k = as.numeric(k_bf),
       mcse_truncgpd = tg$mcse, bias_truncgpd = tg$bias,
       k_truncgpd = tg$k,
       mcse_truncgpd_hybrid = tg$mcse_hybrid,
       bias_truncgpd_hybrid = tg$bias_hybrid)
}

## loobench R/bridge_normal.R, lines 38-46
.lb_pd_repair <- function(m) {
  ev <- eigen((m + t(m)) / 2, symmetric = TRUE)
  top <- max(ev$values)
  if (!is.finite(top) || top <= 0) return(NULL)
  floor_v <- 1e-10 * top
  vals <- pmax(ev$values, floor_v)
  out <- ev$vectors %*% (vals * t(ev$vectors))
  list(m = (out + t(out)) / 2, n_floored = sum(ev$values < floor_v))
}

## loobench R/bridge_normal.R, lines 51-55
.lb_psd_power <- function(m, power) {
  ev <- eigen((m + t(m)) / 2, symmetric = TRUE)
  out <- ev$vectors %*% (pmax(ev$values, 0)^power * t(ev$vectors))
  (out + t(out)) / 2
}

## loobench R/bridge_normal.R, lines 64-76
.lb_geometric_mean_psd <- function(a, b) {
  out <- tryCatch({
    ah <- .lb_psd_power(a, 0.5)
    aih <- .lb_psd_power(a, -0.5)
    g <- ah %*% .lb_psd_power(aih %*% b %*% aih, 0.5) %*% ah
    (g + t(g)) / 2
  }, error = function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) {
    out <- (a + b) / 2
    attr(out, "fallback") <- TRUE
  }
  out
}

## loobench R/bridge_normal.R, lines 83-91
.lb_score_covariance <- function(grads) {
  g <- as.matrix(grads)
  ok <- rowSums(!is.finite(g)) == 0L
  if (sum(ok) < ncol(g) + 1L) return(NULL)
  g <- sweep(g[ok, , drop = FALSE], 2L, colMeans(g[ok, , drop = FALSE]))
  prec <- .lb_pd_repair(crossprod(g) / nrow(g))
  if (is.null(prec)) return(NULL)
  list(sigma = .lb_psd_power(prec$m, -1), n_used = sum(ok))
}

## loobench R/bridge_normal.R, lines 102-133
.lb_normal_proposal <- function(upars_fit, seed, n, grads = NULL,
                                proposal_fit = c("sample", "hybrid")) {
  proposal_fit <- match.arg(proposal_fit)
  u <- as.matrix(upars_fit)
  p <- ncol(u)
  mu <- colMeans(u)
  smp <- .lb_pd_repair(stats::cov(u))
  if (is.null(smp))
    stop(".lb_normal_proposal: the draw covariance has no positive eigenvalue")
  sigma <- smp$m
  combiner <- "none"
  n_grad <- NA_integer_
  if (identical(proposal_fit, "hybrid") && !is.null(grads)) {
    sc <- .lb_score_covariance(grads)
    if (!is.null(sc)) {
      gm <- .lb_geometric_mean_psd(smp$m, sc$sigma)
      combiner <- if (isTRUE(attr(gm, "fallback"))) "geometric_fallback_arithmetic"
                  else "geometric"
      attr(gm, "fallback") <- NULL
      sigma <- .lb_pd_repair(gm)$m
      n_grad <- as.integer(sc$n_used)
    }
  }
  ch <- chol(sigma)
  z <- .lb_with_seed(seed, matrix(stats::rnorm(n * p), n, p))
  draws <- sweep(z %*% ch, 2L, mu, "+")
  colnames(draws) <- colnames(u)
  list(draws = draws, mean = mu, sigma = sigma,
       n_floored = smp$n_floored, proposal_fit = proposal_fit,
       combiner = combiner, n_gradients_used = n_grad,
       log_density = function(x) .lb_mvn_log_density(x, mu, sigma))
}

## loobench R/bridge_normal.R, lines 139-146
.lb_log_density_rows <- function(ad, bundle, upars) {
  out <- tryCatch(as.numeric(ad$log_density(bundle, upars)),
                  error = function(e) NULL)
  if (!is.null(out) && length(out) == nrow(upars)) return(out)
  vapply(seq_len(nrow(upars)), function(i)
    tryCatch(as.numeric(ad$log_density(bundle, upars[i, , drop = FALSE]))[1L],
             error = function(e) NA_real_), numeric(1))
}

## loobench R/bridge_normal.R, lines 165-201
.lb_normal_bridge_logml <- function(bundle, ad, n_chains, seed, cfg = list(),
                                    proposal_fit = c("sample", "hybrid")) {
  proposal_fit <- match.arg(proposal_fit)
  u <- ad$unconstrain(bundle)
  s <- nrow(u)
  nc <- as.integer(n_chains)
  if (!is.na(nc) && nc >= 2L && s %% nc == 0L) {
    per <- s %/% nc
    n_fit_ch <- nc %/% 2L
    fit_idx <- seq_len(per * n_fit_ch)
    n_iter_ch <- nc - n_fit_ch
  } else {
    fit_idx <- seq_len(round(s / 2))
    n_iter_ch <- 1L
  }
  u_fit <- u[fit_idx, , drop = FALSE]
  u_it <- u[-fit_idx, , drop = FALSE]
  grads <- if (identical(proposal_fit, "hybrid") && is.function(ad$grad_log_density))
             tryCatch(ad$grad_log_density(bundle, u_fit), error = function(e) NULL)
           else NULL
  prop <- .lb_normal_proposal(u_fit, seed, nrow(u_it), grads = grads,
                              proposal_fit = proposal_fit)
  q11 <- .lb_log_density_rows(ad, bundle, u_it)
  q21 <- .lb_log_density_rows(ad, bundle, prop$draws)
  bad21 <- is.na(q21)
  q21[bad21] <- -Inf
  q12 <- prop$log_density(u_it)
  q22 <- prop$log_density(prop$draws)
  bb <- .lb_bb_core(q11 - q12, q21 - q22, cfg, n_chains1 = n_iter_ch,
                    iid2 = TRUE)
  list(logml = as.numeric(bb$log_r), mcse_logml = as.numeric(bb$mcse),
       n_iter = as.integer(bb$n_iter), converged = isTRUE(bb$converged),
       gap = as.numeric(bb$gap), error = bb$error,
       n_neg_inf = sum(bad21 | q21 == -Inf), n_floored = prop$n_floored,
       proposal_fit = proposal_fit, combiner = prop$combiner,
       n_gradients_used = prop$n_gradients_used)
}

## loobench R/est_integrated.R, lines 495-502
.lb_mvn_log_density <- function(B, mu, Sigma) {
  B <- if (is.null(dim(B))) matrix(B, 1L, length(mu)) else as.matrix(B)
  m <- length(mu)
  ch <- chol(Sigma)
  d <- sweep(B, 2L, mu, "-")
  z <- backsolve(ch, t(d), transpose = TRUE)
  -0.5 * (m * log(2 * pi) + 2 * sum(log(diag(ch))) + colSums(z * z))
}

## loobench R/fit.R, lines 148-148
LB_RNG_KIND <- c("L'Ecuyer-CMRG", "Inversion", "Rejection")

## loobench R/fit.R, lines 151-164
.lb_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  old_kind <- RNGkind()
  on.exit({
    RNGkind(old_kind[1L], old_kind[2L], old_kind[3L])
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  set.seed(as.integer(seed) %% .Machine$integer.max, kind = LB_RNG_KIND[1L],
           normal.kind = LB_RNG_KIND[2L], sample.kind = LB_RNG_KIND[3L])
  expr
}

## loobench R/fit.R, lines 405-413
.lb_compile_args <- function(local_stan,
                             pre_v1 = "compile" %in%
                               names(formals(cmdstanr::cmdstan_model))) {
  if (isTRUE(pre_v1))
    list(local_stan, force_recompile = TRUE,
         compile_model_methods = TRUE, stanc_options = list("O1"))
  else
    list(local_stan, stanc_options = list("O1"))
}

## loobench R/diagnostics.R, lines 54-65
pareto_k <- function(x, r_eff = 1, tail = "right", are_log_weights = FALSE) {
  x <- as.numeric(x)
  if (!length(x) || any(!is.finite(x))) return(NA_real_)
  if (!are_log_weights && all(x == 0)) return(NA_real_)
  kh <- tryCatch(
    suppressWarnings(posterior::pareto_khat(
      x, r_eff = r_eff, tail = tail, are_log_weights = are_log_weights)),
    error = function(e) NA_real_)
  ## `posterior::pareto_khat` returns a bare numeric in some versions and
  ## a list with `$khat` in others.
  if (is.list(kh)) as.numeric(kh$khat) else as.numeric(kh)
}

## loobench R/diagnostics.R, lines 81-101
reff_per_obs <- function(ll, n_chains = 1L) {
  ll <- as.matrix(ll)
  s_full <- nrow(ll)
  n <- ncol(ll)
  n_ch <- if (is.null(n_chains) || !is.finite(n_chains) || n_chains < 1L)
            1L else as.integer(n_chains)
  iter_per_chain <- s_full %/% n_ch
  vapply(seq_len(n), function(i) {
    ll_i <- ll[, i]
    if (!all(is.finite(ll_i))) return(1.0)
    lik_i <- exp(ll_i - max(ll_i))
    ess <- if (n_ch > 1L && s_full %% n_ch == 0L) {
      arr <- matrix(lik_i, nrow = iter_per_chain, ncol = n_ch)
      tryCatch(as.numeric(posterior::ess_mean(arr)), error = function(e) NA_real_)
    } else {
      tryCatch(as.numeric(posterior::ess_mean(lik_i)), error = function(e) NA_real_)
    }
    r_eff <- if (is.finite(ess) && ess > 0) ess / s_full else 1.0
    max(min(r_eff, 1.0), 1e-3)
  }, numeric(1))
}

## loobench R/diagnostics.R, lines 121-134
reff_correct_mcse <- function(mcse_r1, r_eff, na_on_fail = TRUE) {
  mcse_r1 <- as.numeric(mcse_r1)
  r_eff <- rep_len(as.numeric(r_eff), length(mcse_r1))
  vapply(seq_along(mcse_r1), function(i) {
    m <- mcse_r1[i]
    re <- r_eff[i]
    fail <- if (na_on_fail) NA_real_ else m
    if (!is.finite(m) || m <= 0) return(fail)
    if (!is.finite(re) || re <= 0) return(fail)
    v <- log(1 + (exp(m^2) - 1) / re)
    if (!is.finite(v) || v < 0) return(fail)
    sqrt(v)
  }, numeric(1))
}

## loobench R/diagnostics.R, lines 156-179
.lb_mcse_log_mean_exp <- function(log_lik, n_chains = 1L) {
  log_lik <- as.numeric(log_lik)
  s_draws <- length(log_lik)
  if (!s_draws || any(!is.finite(log_lik))) return(NA_real_)
  n_ch <- as.integer(n_chains)
  if (is.na(n_ch) || n_ch < 1L || s_draws %% n_ch != 0L) n_ch <- 1L
  elpd <- matrixStats::logSumExp(log_lik) - log(s_draws)
  e_epd <- exp(elpd)
  if (!is.finite(e_epd) || e_epd <= 0) return(NA_real_)
  lik <- exp(log_lik)
  if (any(!is.finite(lik))) return(NA_real_)
  r_eff <- tryCatch({
    if (n_ch > 1L)
      as.numeric(posterior::ess_mean(matrix(log_lik, s_draws %/% n_ch, n_ch))) /
        s_draws
    else as.numeric(posterior::ess_mean(log_lik)) / s_draws
  }, error = function(e) 1.0)
  if (!is.finite(r_eff) || r_eff <= 0) r_eff <- 1.0
  var_epd <- sum((lik - e_epd)^2) / (s_draws^2 * r_eff)
  if (!is.finite(var_epd) || var_epd < 0) return(NA_real_)
  var_log <- log(1 + var_epd / e_epd^2)
  if (!is.finite(var_log) || var_log < 0) return(NA_real_)
  sqrt(var_log)
}

## loobench R/diagnostics.R, lines 208-221
psis_elpd_fold <- function(log_lik, n_chains = 1L) {
  was_vector <- is.null(dim(log_lik))
  ll <- if (was_vector) matrix(as.numeric(log_lik), ncol = 1L) else as.matrix(log_lik)
  ## The Pareto-k warning is suppressed: `k` is returned, and a per-fold
  ## call would otherwise emit one warning per fold.
  lo <- suppressWarnings(loo::loo(ll, r_eff = 1, save_psis = TRUE))
  elpd <- as.numeric(lo$pointwise[, "elpd_loo"])
  mcse_reff1 <- as.numeric(lo$pointwise[, "mcse_elpd_loo"])
  k <- as.numeric(lo$diagnostics$pareto_k)
  r_eff <- reff_per_obs(ll, n_chains = n_chains)
  mcse <- reff_correct_mcse(mcse_reff1, r_eff, na_on_fail = TRUE)
  list(elpd = elpd, mcse = mcse, k = k, mcse_reff1 = mcse_reff1,
       r_eff = r_eff, psis_loo = lo)
}

## loobench R/diagnostics.R, lines 227-233
trunc_pareto_mean <- function(k, y) {
  if (abs(k - 1) < 1e-3) {
    log(y) / (1 - 1 / y)
  } else {
    (1 / (k - 1)) * (y^(1 - 1 / k) - 1) / (1 - y^(-1 / k))
  }
}

## loobench R/diagnostics.R, lines 237-244
trunc_pareto_var <- function(k, y) {
  mu <- trunc_pareto_mean(k, y)
  if (abs(k - 0.5) < 1e-3) {
    2 * log(y) / (1 - 1 / y^2) - mu^2
  } else {
    (1 / (2 * k - 1)) * (y^(2 - 1 / k) - 1) / (1 - y^(-1 / k)) - mu^2
  }
}

## loobench R/diagnostics.R, lines 252-267
trunc_pareto_raw_moments <- function(k, y, m_max = 4L) {
  if (!is.finite(k) || !is.finite(y) || y <= 1) return(rep(NA_real_, m_max))
  denom <- 1 - y^(-1 / k)
  if (!is.finite(denom) || denom <= 0) return(rep(NA_real_, m_max))
  out <- numeric(m_max)
  for (m in seq_len(m_max)) {
    if (abs(m * k - 1) < 1e-6) {
      out[m] <- m * log(y) / (1 - y^(-m))
    } else {
      num <- y^(m - 1 / k) - 1
      out[m] <- num / ((m * k - 1) * denom)
    }
    if (!is.finite(out[m])) return(rep(NA_real_, m_max))
  }
  out
}

## loobench R/diagnostics.R, lines 282-337
trunc_gpd_mean_var <- function(k, sigma, t_physical, S) {
  if (!is.finite(k) || !is.finite(sigma) || !is.finite(S) ||
        sigma <= 0 || S <= 1)
    return(list(mu = NA_real_, var = NA_real_, tau_eff = NA_real_,
                truncation_source = NA_character_))
  no_phys <- !is.finite(t_physical) || t_physical <= 0
  ## Exponential limit at small k.
  if (abs(k) < 1e-3) {
    if (no_phys) {
      return(list(mu = sigma, var = sigma^2, tau_eff = NA_real_,
                  truncation_source = "S^k (no physical bound)"))
    }
    a <- t_physical / sigma
    if (a > 700) {
      mu <- sigma
      e2 <- 2 * sigma^2
    } else {
      mu <- sigma * (1 - (1 + a) * exp(-a)) / (1 - exp(-a))
      e2 <- sigma^2 * (2 - (2 + 2 * a + a^2) * exp(-a)) / (1 - exp(-a))
    }
    return(list(mu = mu, var = e2 - mu^2, tau_eff = NA_real_,
                truncation_source = "exponential_limit"))
  }
  if (k < 0) {
    return(list(mu = sigma / (1 - k),
                var = sigma^2 / ((1 - k)^2 * (1 - 2 * k)),
                tau_eff = Inf,
                truncation_source = "none (bounded tail)"))
  }
  tau_sk <- S^k
  if (no_phys) {
    tau_phys <- Inf
    src_base <- "S^k (no physical bound)"
  } else {
    tau_phys <- 1 + k * t_physical / sigma
    src_base <- NA_character_
  }
  if (!is.finite(tau_phys) || tau_phys <= 1) {
    tau_eff <- tau_sk
    src <- if (no_phys) src_base else "S^k (tau_phys not usable)"
  } else if (tau_phys <= tau_sk) {
    tau_eff <- tau_phys
    src <- "physical (y_max informative)"
  } else {
    tau_eff <- tau_sk
    src <- "S^k (y_max too loose)"
  }
  mu_te <- trunc_pareto_mean(k, tau_eff)
  var_te <- trunc_pareto_var(k, tau_eff)
  if (!is.finite(mu_te) || !is.finite(var_te))
    return(list(mu = NA_real_, var = NA_real_, tau_eff = tau_eff,
                truncation_source = src))
  list(mu = (sigma / k) * (mu_te - 1),
       var = (sigma / k)^2 * var_te,
       tau_eff = tau_eff, truncation_source = src)
}

## loobench R/diagnostics.R, lines 347-403
trunc_gpd_raw_moments <- function(k, sigma, t_physical, S) {
  na_ret <- list(mu = NA_real_, m2 = NA_real_, m3 = NA_real_, m4 = NA_real_,
                 tau_eff = NA_real_, truncation_source = NA_character_)
  if (!is.finite(k) || !is.finite(sigma) || !is.finite(S) ||
        sigma <= 0 || S <= 1)
    return(na_ret)
  no_phys <- !is.finite(t_physical) || t_physical <= 0
  if (abs(k) < 1e-3) {
    if (no_phys) {
      return(list(mu = sigma, m2 = 2 * sigma^2, m3 = 6 * sigma^3,
                  m4 = 24 * sigma^4, tau_eff = NA_real_,
                  truncation_source = "S^k (no physical bound)"))
    }
    a <- t_physical / sigma
    denom_exp <- 1 - exp(-a)
    if (!is.finite(denom_exp) || denom_exp <= 0)
      return(c(na_ret, list(truncation_source = "exponential_limit (degenerate)")))
    e_y_j <- function(j) {
      sigma^j * gamma(j + 1) * stats::pgamma(a, j + 1, lower.tail = TRUE) / denom_exp
    }
    return(list(mu = e_y_j(1), m2 = e_y_j(2), m3 = e_y_j(3), m4 = e_y_j(4),
                tau_eff = NA_real_, truncation_source = "exponential_limit"))
  }
  tau_sk <- S^k
  if (no_phys) {
    tau_phys <- Inf
    src_base <- "S^k (no physical bound)"
  } else {
    tau_phys <- 1 + k * t_physical / sigma
    src_base <- NA_character_
  }
  if (!is.finite(tau_phys) || tau_phys <= 1) {
    tau_eff <- tau_sk
    src <- if (no_phys) src_base else "S^k (tau_phys not usable)"
  } else if (tau_phys <= tau_sk) {
    tau_eff <- tau_phys
    src <- "physical (y_max informative)"
  } else {
    tau_eff <- tau_sk
    src <- "S^k (y_max too loose)"
  }
  e_tau <- trunc_pareto_raw_moments(k, tau_eff, m_max = 4L)
  if (any(!is.finite(e_tau)))
    return(list(mu = NA_real_, m2 = NA_real_, m3 = NA_real_, m4 = NA_real_,
                tau_eff = tau_eff, truncation_source = src))
  s_over_k <- sigma / k
  e_tau_full <- c(1, e_tau)
  e_y <- numeric(4)
  for (j in 1:4) {
    coef <- vapply(0:j,
                   function(i) choose(j, i) * (-1)^(j - i) * e_tau_full[i + 1],
                   numeric(1))
    e_y[j] <- s_over_k^j * sum(coef)
  }
  list(mu = e_y[1], m2 = e_y[2], m3 = e_y[3], m4 = e_y[4],
       tau_eff = tau_eff, truncation_source = src)
}

## loobench R/diagnostics.R, lines 412-428
higher_order_jensen <- function(m1, m2, m3, m4, S) {
  na_ret <- list(bias_higher = NA_real_, leading = NA_real_, skew = NA_real_,
                 kurt = NA_real_, var_X = NA_real_, mu3_X = NA_real_)
  if (!is.finite(m1) || !is.finite(m2) || !is.finite(m3) ||
        !is.finite(m4) || !is.finite(S) || m1 <= 0 || S <= 1)
    return(na_ret)
  var_x <- m2 - m1^2
  mu3_x <- m3 - 3 * m1 * m2 + 2 * m1^3
  if (!is.finite(var_x) || var_x <= 0) return(na_ret)
  tau2 <- var_x / (S * m1^2)
  if (!is.finite(tau2) || tau2 <= 0) return(na_ret)
  leading <- -tau2 / 2
  skew <- mu3_x / (3 * S^2 * m1^3)
  kurt <- -3 * tau2^2 / 4
  list(bias_higher = leading + skew + kurt, leading = leading,
       skew = skew, kurt = kurt, var_X = var_x, mu3_X = mu3_x)
}

## loobench R/diagnostics.R, lines 470-587
truncgpd_mcse <- function(log_lik, y_max = Inf, y_max_alt = NULL,
                          n_chains = 1L) {
  na_row <- list(
    mcse_truncgpd = NA_real_,
    bias_jensen_truncgpd = NA_real_,
    truncation_source = NA_character_,
    k_truncgpd = NA_real_,
    sigma_truncgpd = NA_real_,
    threshold_u = NA_real_,
    tau_eff = NA_real_,
    mcse_truncgpd_hybrid = NA_real_,
    bias_jensen_truncgpd_hybrid = NA_real_,
    truncation_source_hybrid = NA_character_,
    tau_eff_hybrid = NA_real_,
    bias_jensen_higher_truncgpd = NA_real_,
    bias_jensen_higher_truncgpd_hybrid = NA_real_,
    skew_corr_truncgpd = NA_real_,
    kurt_corr_truncgpd = NA_real_,
    skew_corr_truncgpd_hybrid = NA_real_,
    kurt_corr_truncgpd_hybrid = NA_real_)
  log_lik <- as.numeric(log_lik)
  if (length(log_lik) < 50L || any(!is.finite(log_lik))) return(na_row)
  x <- exp(log_lik - max(log_lik))              # numerical stability
  scale <- exp(max(log_lik))
  s_draws <- length(x)
  ## Fit the GPD to the top-M exceedances with loo's wip prior, for
  ## consistency with the internals of `posterior::pareto_khat`.
  m <- as.integer(min(floor(0.2 * s_draws), floor(3 * sqrt(s_draws))))
  m <- max(m, 5L)
  x_sorted <- sort(x)
  u_thr <- x_sorted[s_draws - m]
  exc <- x_sorted[(s_draws - m + 1):s_draws] - u_thr
  exc <- exc[exc > 0]
  if (length(exc) < 5L) return(na_row)
  fit <- tryCatch(loo::gpdfit(exc, wip = TRUE), error = function(e) NULL)
  if (is.null(fit) || !is.finite(fit$k) || !is.finite(fit$sigma))
    return(na_row)
  ## Both bounds share the GPD fit; only `t_physical` differs.  The S^k
  ## cap is applied per bound inside `trunc_gpd_mean_var()`.
  y_max_norm <- y_max / scale
  t_physical <- y_max_norm - u_thr
  x_bulk <- x_sorted[1:(s_draws - m)]
  mu_bulk <- mean(x_bulk)
  v_bulk <- stats::var(x_bulk)
  p_tail <- m / s_draws
  m1_bulk <- mean(x_bulk)
  m2_bulk <- mean(x_bulk^2)
  m3_bulk <- mean(x_bulk^3)
  m4_bulk <- mean(x_bulk^4)
  combine_to_mcse <- function(mv) {
    if (!is.finite(mv$mu) || !is.finite(mv$var))
      return(list(mcse_log = NA_real_, bias_pred = NA_real_, ok = FALSE))
    mu_tail_x <- u_thr + mv$mu
    v_tail_x <- mv$var
    mu_total <- (1 - p_tail) * mu_bulk + p_tail * mu_tail_x
    between <- (1 - p_tail) * p_tail * (mu_bulk - mu_tail_x)^2
    v_total <- (1 - p_tail) * v_bulk + p_tail * v_tail_x + between
    if (!is.finite(v_total) || v_total <= 0 || mu_total <= 0)
      return(list(mcse_log = NA_real_, bias_pred = NA_real_, ok = FALSE))
    mcse_log <- sqrt(v_total / s_draws) / mu_total
    list(mcse_log = mcse_log, bias_pred = -mcse_log^2 / 2, ok = TRUE)
  }
  combine_to_higher_jensen <- function(mv_raw) {
    if (!is.finite(mv_raw$mu) || !is.finite(mv_raw$m2) ||
          !is.finite(mv_raw$m3) || !is.finite(mv_raw$m4))
      return(higher_order_jensen(NA, NA, NA, NA, s_draws))
    ## E[X_tail^j] = sum_i C(j,i) u_thr^(j-i) E[Y^i].
    e_y_full <- c(1, mv_raw$mu, mv_raw$m2, mv_raw$m3, mv_raw$m4)
    m_tail <- numeric(4)
    for (j in 1:4) {
      coef <- vapply(0:j,
                     function(i) choose(j, i) * u_thr^(j - i) * e_y_full[i + 1],
                     numeric(1))
      m_tail[j] <- sum(coef)
    }
    m1_x <- (1 - p_tail) * m1_bulk + p_tail * m_tail[1]
    m2_x <- (1 - p_tail) * m2_bulk + p_tail * m_tail[2]
    m3_x <- (1 - p_tail) * m3_bulk + p_tail * m_tail[3]
    m4_x <- (1 - p_tail) * m4_bulk + p_tail * m_tail[4]
    higher_order_jensen(m1_x, m2_x, m3_x, m4_x, s_draws)
  }
  ## ---- Physical bound branch ----
  mv_phys <- trunc_gpd_mean_var(fit$k, fit$sigma, t_physical, S = s_draws)
  m_phys <- combine_to_mcse(mv_phys)
  mv_phys_raw <- trunc_gpd_raw_moments(fit$k, fit$sigma, t_physical, S = s_draws)
  h_phys <- combine_to_higher_jensen(mv_phys_raw)
  ## ---- Hybrid bound branch (optional, shares the GPD fit) ----
  use_alt <- !is.null(y_max_alt) && is.finite(y_max_alt) && y_max_alt > 0
  if (use_alt) {
    t_alt <- y_max_alt / scale - u_thr
    mv_alt <- trunc_gpd_mean_var(fit$k, fit$sigma, t_alt, S = s_draws)
    m_alt <- combine_to_mcse(mv_alt)
    mv_alt_raw <- trunc_gpd_raw_moments(fit$k, fit$sigma, t_alt, S = s_draws)
    h_alt <- combine_to_higher_jensen(mv_alt_raw)
  } else {
    mv_alt <- list(tau_eff = NA_real_, truncation_source = NA_character_)
    m_alt <- list(mcse_log = NA_real_, bias_pred = NA_real_, ok = FALSE)
    h_alt <- higher_order_jensen(NA, NA, NA, NA, s_draws)
  }
  list(
    mcse_truncgpd = m_phys$mcse_log,
    bias_jensen_truncgpd = m_phys$bias_pred,
    truncation_source = mv_phys$truncation_source,
    k_truncgpd = fit$k,
    sigma_truncgpd = fit$sigma,
    threshold_u = u_thr * scale,                 # de-stabilized
    tau_eff = mv_phys$tau_eff,
    mcse_truncgpd_hybrid = m_alt$mcse_log,
    bias_jensen_truncgpd_hybrid = m_alt$bias_pred,
    truncation_source_hybrid = mv_alt$truncation_source,
    tau_eff_hybrid = mv_alt$tau_eff,
    bias_jensen_higher_truncgpd = h_phys$bias_higher,
    bias_jensen_higher_truncgpd_hybrid = h_alt$bias_higher,
    skew_corr_truncgpd = h_phys$skew,
    kurt_corr_truncgpd = h_phys$kurt,
    skew_corr_truncgpd_hybrid = h_alt$skew,
    kurt_corr_truncgpd_hybrid = h_alt$kurt)
}

## loobench R/diagnostics.R, lines 649-654
.lb_hybrid_bound <- function(y_max, y_max_full) {
  y_max_full <- suppressWarnings(as.numeric(y_max_full))
  if (length(y_max_full) == 1L && is.finite(y_max_full) && y_max_full > 0)
    min(as.numeric(y_max), y_max_full)
  else NA_real_
}

