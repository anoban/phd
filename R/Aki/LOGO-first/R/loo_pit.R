## ── Analytic LOO-PIT and its ECDF plots ────────────────────────────
##
## The probability integral transform of each observation under its own
## leave-one-out predictive distribution.  For a calibrated model the
## LOO-PIT values are uniform on (0, 1); a systematic departure shows
## how the model is wrong (too wide, too narrow, shifted).
##
## Here the LOO predictive distribution is a finite mixture of normal
## distributions over the posterior draws, so the LOO-PIT is computed in
## closed form,
##
##   PIT_id = E[ Phi((y_id - mu_id) / sigma_id) ]
##
## over the leave-one-out posterior.  Stan emits the per-draw Phi(.)
## term (`loo_pit_draw`) and loo::E_loo() takes that expectation with
## the PSIS weights.  No posterior predictive draws are simulated, so
## there is no y_rep Monte Carlo noise in the diagnostic.
##

library("cmdstanr")
library("loo")
library("bayesplot")
library("ggplot2")
library("patchwork")

#' Analytic LOO-PIT values for a fitted model.
#'
#' Runs the generated quantities again on the saved draws, without new
#' sampling.  `loo_pit_draw` is N x D per draw.
#'
#' @param model one of "M1".."M4", "N1".."N4"
#' @param bundle output of build_anoban_data()
#' @param fit the CmdStanMCMC object from fit_anoban()
#' @param stan_dir directory of the .stan files
#' @return list(pit = N x D matrix, pareto_k, responses)

anoban_loo_pit <- function(model, bundle, fit, stan_dir = "stan", responses = c("F00727", "F00679", "F00709")) {
  spec <- ANOBAN_MODELS[[model]]
  if (is.null(spec)) stop("unknown model: ", model)

  standata <- bundle$standata
  standata$cor_t <- spec$cor_t
  if (!is.na(spec$cor_p)) standata$cor_p <- spec$cor_p
  standata$compute_loo <- 1L
  standata$compute_pit <- 1L

  mod <- anoban_stan_model(file.path(stan_dir, spec$stan), stan_dir)

  dr <- fit$draws()
  gq <- mod$generate_quantities(fitted_params = dr, data = standata, seed = 1)

  N <- standata$N; D <- standata$D
  ll <- as.matrix(gq$draws(variables = "log_lik", format = "draws_matrix"))
  pd <- as.matrix(gq$draws(variables = "loo_pit_draw", format = "draws_matrix"))

  ## LOO-PIT_id is the expectation of the predictive CDF value under the
  ## leave-one-out posterior, so it is exactly what E_loo() computes
  ## from the PSIS weights.  loo::psis() takes log ratios, i.e. -log_lik.
  psis_obj <- loo::psis(-ll, r_eff = 1)
  pit <- vapply(seq_len(D), function(d) {
    cols <- sprintf("loo_pit_draw[%d,%d]", seq_len(N), d)
    loo::E_loo(pd[, cols, drop = FALSE], psis_obj, type = "mean", log_ratios = -ll)$value
  }, numeric(N))
  colnames(pit) <- responses

  list(pit = pit, pareto_k = psis_obj$diagnostics$pareto_k, responses = responses, model = model)
}



#' LOO-PIT ECDF plots, one panel per response.
#'
#' Uses bayesplot's `method = "correlated"` uniformity test (Tesso &
#' Vehtari 2026), which accounts for the LOO-PIT values being
#' correlated.  `plot_diff = TRUE` plots the ECDF minus its uniform
#' expectation.
#'
#' All panels share a y-axis so the magnitudes are comparable.
#'
#' `POT` and `PIET` require LOO-PIT values strictly inside (0, 1); the
#' values arrive that way, because `loo_pit_vec()` clamps each draw's
#' CDF value where `Phi()` would underflow to 0 or round to 1.
#'
#' @param pit_obj output of anoban_loo_pit()
#' @param test "POT" (default), "PRIT" or "PIET"
#' @return a patchwork object (a ggplot)

plot_anoban_loo_pit <- function(pit_obj, test = "POT", title_prefix = NULL) {

  prefix <- if (is.null(title_prefix)) pit_obj$model else title_prefix

  plots <- lapply(seq_along(pit_obj$responses), function(d) {
    ppc_loo_pit_ecdf(pit = pit_obj$pit[, d], method = "correlated", test = test, plot_diff = TRUE) +
      ggtitle(sprintf("%s - %s", prefix, pit_obj$responses[d]))
  })

  ## Shared y range so the panels are comparable.  `&` in patchwork
  ## applies to every panel at once, and coord_cartesian zooms rather
  ## than filters, so the simultaneous bands are not clipped.
  rng <- range(vapply(plots, function(g) layer_scales(g)$y$get_limits(), numeric(2)))
  wrap_plots(plots, nrow = 1) & coord_cartesian(ylim = rng)
}

