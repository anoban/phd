## ── Compile and fit the marginal models ─────────────────────────────

suppressPackageStartupMessages(library(cmdstanr))

## M1/M3 carry the phylogenetic effect (pruning model); M2/M4 do not.
## cor_* = 1 means the species effect is correlated ACROSS responses,
## i.e. brms `(1|q|...)` rather than `(1|...)`.
## The N models are the M models with the `state` fixed effect dropped:
## `~ 1` instead of `~ state`, everything else identical.  Because
## `state` is species-constant, removing it does not simply delete a
## term — that variation moves into the species-level effects, so the
## M-vs-N elpd difference measures what `state` contributes beyond the
## phylogenetic and iid species structure.
ANOBAN_MODELS <- local({

  spec <- function(stan, cor_t, cor_p, predictor, rhs) {
    list(stan = stan, cor_t = cor_t, cor_p = cor_p, predictor = predictor,
         brms = paste0("mvbind(...) ~ ", predictor, " + ", rhs))
  }

  g1 <- "(1 | gr(binominal, cov = corrmat)) + (1 | taxa)"
  g2 <- "(1 | taxa)"
  g3 <- "(1|p| gr(binominal, cov = corrmat)) + (1|q| taxa)"
  g4 <- "(1|q| taxa)"

  list(
    M1 = spec("anoban_pruning.stan",  0L, 0L, "state", g1),
    M2 = spec("anoban_marginal.stan", 0L, NA, "state", g2),
    M3 = spec("anoban_pruning.stan",  1L, 1L, "state", g3),
    M4 = spec("anoban_marginal.stan", 1L, NA, "state", g4),

    N1 = spec("anoban_pruning.stan",  0L, 0L, "1", g1),
    N2 = spec("anoban_marginal.stan", 0L, NA, "1", g2),
    N3 = spec("anoban_pruning.stan",  1L, 1L, "1", g3),
    N4 = spec("anoban_marginal.stan", 1L, NA, "1", g4)
  )
})

#' Compile a model, recompiling when an #include'd file has changed.
#'
#' cmdstanr decides whether the binary is stale from the .stan file
#' alone.  Both models get their functions from
#' `anoban_functions.stanfunctions` via #include, so editing that file
#' leaves the old executable in place and the change silently does
#' nothing.  Compare mtimes against every source the model depends on.
anoban_stan_model <- function(stan_file, stan_dir) {
  ## Let cmdstanr do its own up-to-date check first; this is cheap when
  ## the binary is current, and it gives the real executable path
  ## (exe_file() is empty before a model is compiled).
  mod <- cmdstan_model(stan_file, include_paths = stan_dir)
  includes <- list.files(stan_dir, pattern = "[.]stanfunctions$", full.names = TRUE)

  if (length(includes) &&
      any(file.mtime(includes) > file.mtime(mod$exe_file()))) {
    message("an #include is newer than the binary; recompiling")
    mod <- cmdstan_model(stan_file, include_paths = stan_dir, force_recompile = TRUE)
  }

  mod
}

#' Fit one model.
#'
#' @param model one of "M1".."M4"
#' @param bundle output of build_anoban_data()
#' @param stan_dir directory holding the .stan files
#' @param chains,iter_warmup,iter_sampling sampler settings
#' @param compute_loo emit log_lik (needed for `fit$loo()`)
#' @return a CmdStanMCMC object
fit_anoban <- function(model, bundle, stan_dir = "./stan", chains = 8, iter_warmup = 1000, iter_sampling = 1000, seed = 20260830, compute_loo = TRUE, ...) {
  spec <- ANOBAN_MODELS[[model]]
  if (is.null(spec)) stop("unknown model: ", model)

  standata <- bundle$standata
  standata$cor_t <- spec$cor_t
  standata$compute_loo <- as.integer(compute_loo)
  if (!is.na(spec$cor_p)) standata$cor_p <- spec$cor_p

  mod <- anoban_stan_model(file.path(stan_dir, spec$stan), stan_dir)

  mod$sample(data = standata, chains = chains, parallel_chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling, seed = seed, init = 0.3, ...)
}
