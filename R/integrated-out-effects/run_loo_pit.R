#!/usr/bin/env Rscript
## ── LOO-PIT plots for fitted models ─────────────────────────────────
##
## Run AFTER run_all.R, which saves the fits this reads.
##
##   Rscript run_loo_pit.R                  # all models found in output/
##   Rscript run_loo_pit.R --models M3      # one model
##   Rscript run_loo_pit.R --raw            # if the fits used --raw
##
## Outputs, under output/:
##   loo_pit_<model>.pdf
##
## Sourced from an interactive session the plots are also drawn on the
## current device, and `plots[[m]]` keeps each one for re-display or
## further tweaking.  They are NOT printed under `Rscript`, which has
## no display and would otherwise dump them into a stray Rplots.pdf.
##
## Reading the plots: the line is the LOO-PIT ECDF minus its uniform
## expectation, so a calibrated model stays flat near zero and inside
## the band.  A rising-then-falling (S) shape means PIT values pile up
## near 0.5 — the predictive is too WIDE.  The mirror image means it is
## too narrow.  A shift up or down means the predictive is biased.

## Run from the package root under Rscript; when sourced from an R
## session, stay in the working directory (which should be that root).
## Guard the --file= lookup: normalizePath(NA) warns when there is none.

source("./R/build_data.R")
source("./R/fit_models.R")
source("./R/loo_pit.R")

options(mc.cores = 8) # we only used 8 chains in sampling
options(loo.cores = 8)

LOG_RESPONSE <- TRUE
MODELS <- strsplit("M1,M2,M3,M4,N1,N2,N3,N4", ",")[[1]]
MODELS <- MODELS[file.exists(file.path("./output", sprintf("%s_fit.Rds", MODELS)))]
if (!length(MODELS)) stop("no fits in output/ — run run_all.R first")

cat(sprintf("response scale: %s | models: %s\n", if (LOG_RESPONSE) "log" else "original", paste(MODELS, collapse = ", ")))

## data once per distinct right-hand side (M* use `state`, N* use `1`)
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

plots <- list()

for (m in MODELS) {

  cat(sprintf("\n── %s ──\n", m))
  ## saved fit
  fit <- readRDS(file.path("./output", sprintf("%s_fit.Rds", m)))
  ## LOO-PIT
  pit_obj <- anoban_loo_pit(m, anoban_data_for(ANOBAN_MODELS[[m]]$predictor), fit, stan_dir = "./stan")

  ## LOO-PIT plots via bayesplot
  plots[[m]] <- plot_anoban_loo_pit(pit_obj)
  out <- file.path("./output", sprintf("loo_pit_%s.png", m))

  png(filename = out, width = 15, height = 5, units = "in", res = 500, type = "cairo-png") # choosing to use grDevices::png()
  # ggplot2::ggsave(units = "in", device = "pdf", filename = out, plot = plots[[m]], width = 15, height = 5, dpi = 500) # try saving as pdf files
  # ggsave() causes the following error on the server with png files
  # Error in f(...) : Graphics API version mismatch
  # Calls: <Anonymous> -> dev -> png_dev -> f
  print(plots[[m]])
  dev.off()

  cat(sprintf("wrote %s\n", out))
}
