## ── Fit M1-M4 and N1-N4 and compute PSIS-LOO ────────────────────────
##
##   Rscript run_all.R                 # log responses (default)
##   Rscript run_all.R --raw           # original response scale
##   Rscript run_all.R --models M3     # subset
##   Rscript run_all.R --quick         # short chains, for a smoke test
##
## Outputs, under output/:
##   <model>_fit.Rds     CmdStanMCMC object
##   <model>_loo.Rds     loo object
##   loo_compare.txt     model comparison
##
## Per-model elpd_loo, p_loo and the Pareto-k table are printed by
## print(loo); sampler diagnostics are printed alongside.  Both are
## recoverable from the saved objects, so nothing is tabulated to disk.

source("./R/build_data.R")
source("./R/fit_models.R")
source("./R/loo_compare_scales.R")

options(mc.cores = 8)
options(loo.cores = 8)

LOG_RESPONSE <- TRUE
MODELS <- c("M1", "M2", "M3", "M4", # M models with a fixed effect
            "N1", "N2", "N3", "N4") # N models with no fixed effects
WARMUP <- 1000
SAMPLING <- 1000

dir.create("./output", showWarnings = FALSE)

cat(sprintf("response scale: %s\n", if (LOG_RESPONSE) "log" else "original"))

## Build data once per distinct right-hand side: M* use `state`, N* use `1` (intercept only).
anoban_data_for <- local({
  cache <- list()
  function(predictor) {
    key <- paste0("rhs_", predictor)
    if (is.null(cache[[key]])) {
      cache[[key]] <<- build_anoban_data("./data/stan.csv", "./data/FRED4_1301.tre", log_response = LOG_RESPONSE, predictor = predictor)
    }
    cache[[key]]
  }
})

anoban_data <- anoban_data_for(ANOBAN_MODELS[[MODELS[1]]]$predictor)
sd0 <- anoban_data$standata
cat(sprintf("N = %d observations, G = %d species, D = %d responses\n", sd0$N, sd0$G, sd0$D))
cat(sprintf("observations per species: min %d, median %d, max %d (%d singletons)\n", min(sd0$ng), stats::median(sd0$ng), max(sd0$ng), sum(sd0$ng == 1)))
cat(sprintf("predictors: %s\n", paste(anoban_data$predictor_names, collapse = ", ")))

## The original responses, needed for the cross-scale comparison.
Y_orig <- as.matrix(read.csv("data/stan.csv", stringsAsFactors = FALSE)[, c("F00727", "F00679", "F00709")])


for (m in MODELS) {
  cat(sprintf("\n── %s : %s ──\n", m, ANOBAN_MODELS[[m]]$brms))
  md <- anoban_data_for(ANOBAN_MODELS[[m]]$predictor)
  t0 <- proc.time()[["elapsed"]]
  ## Fit one model
  fit <- fit_anoban(m, md, stan_dir = "stan", iter_warmup = WARMUP, iter_sampling = SAMPLING, refresh = max(1, WARMUP %/% 4)) # this does not have the extra arg `sig_figs` for the .csv output decimal places????
  wall <- proc.time()[["elapsed"]] - t0

  ## Convergence diagnostics over the model parameters only.  `b` is
  ## zero-size in the intercept-only models, so leave it out there.
  dvars <- c(if (md$standata$Kc > 0) "b", "Intercept", "sigma", "sd_t", "sd_p")
  s <- fit$summary(variables = intersect(dvars, fit$metadata()$stan_variables))
  d <- fit$diagnostic_summary(quiet = TRUE)
  cat(sprintf("%s: %.0f s | min ess_bulk %.0f | max rhat %.3f | divergent %d | treedepth %d\n", m, wall, min(s$ess_bulk, na.rm = TRUE), max(s$rhat, na.rm = TRUE), sum(d$num_divergent), sum(d$num_max_treedepth)))

  ## PSIS-LOO
  lo <- fit$loo()

  print(lo)

  ## Save fit and loo objects.  Tag the loo with the response scale so
  ## the comparison below can refuse to mix log- and raw-scale fits.
  fit$save_object(file.path("output", sprintf("%s_fit.Rds", m)))
  attr(lo, "log_response") <- LOG_RESPONSE
  saveRDS(lo, file.path("output", sprintf("%s_loo.Rds", m)))

}

## Compare EVERY model with a saved loo in output/, not just the ones
## fitted in this run, so fitting a few more models updates the whole
## comparison without refitting the rest.
loo_files <- list.files("output", pattern = "_loo[.]Rds$", full.names = TRUE)
loos <- lapply(loo_files, readRDS)
names(loos) <- sub("_loo[.]Rds$", "", basename(loo_files))

scale_tag <- vapply(loos, function(x) {
  s <- attr(x, "log_response"); if (is.null(s)) NA else isTRUE(s)
}, logical(1))
if (any(!is.na(scale_tag) & scale_tag != LOG_RESPONSE)) {
  stop("output/ mixes log- and raw-response fits; elpd values are not comparable across scales. Clear output/ or refit on one scale.")
}
if (any(is.na(scale_tag))) {
  warning("some saved loo objects predate the response-scale tag; assuming they match the current run", call. = FALSE)
}

if (length(loos) > 1) {
  cmp <- loo_compare_original_scale(loos, Y_orig, LOG_RESPONSE)
  cat(sprintf("\n── loo_compare over %d models (original response scale) ──\n", length(loos)))
  print(cmp)
  capture.output(print(cmp), file = "output/loo_compare.txt")
}
cat("\nwrote output/\n")
