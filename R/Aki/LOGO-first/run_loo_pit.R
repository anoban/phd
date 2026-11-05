#!/usr/bin/env Rscript
## ── LOO-PIT plots for fitted models ─────────────────────────────────


#----------------------------------------------------------------------
# MAKE SURE THE TEMPDIR IS SET TO SOMEWHERE PERSISTENT IN .Renviron
# TMPDIR=/home/u90963425/ScratchData/Rtmp
# TMP=/home/u90963425/ScratchData/Rtmp
# TEMP=/home/u90963425/ScratchData/Rtmp
#----------------------------------------------------------------------


## Run run_all.R first: this script reads the fits that it saves.
##
##   Rscript run_loo_pit.R                  # all models found in output/
##   Rscript run_loo_pit.R --models M3      # one model
##   Rscript run_loo_pit.R --raw            # if the fits used --raw
##
## Outputs, under output/:
##   loo_pit_<model>.pdf
##
## In an interactive session the script also draws the plots on the
## current device, and `plots[[m]]` keeps each plot for later use.
## Under `Rscript` it does not print them: there is no display, and
## printing would write an extra Rplots.pdf.
##
## Reading the plots: the line is the LOO-PIT ECDF minus its uniform
## expectation, so for a calibrated model it stays near zero and inside
## the band.  An S-shape (first up, then down) means that the PIT values
## are too concentrated near 0.5: the predictive distribution is too
## wide.  The mirror image means that it is too narrow.  A shift up or
## down means that the predictive distribution is biased.  Points at 0
## or 1 are observations far in the tails of their predictive
## distribution; if there are more of them than expected, the tails are
## too thin.

options(mc.cores = 8) # we only used 8 chains in sampling
options(loo.cores = 8)

source("./R/build_data.R")
source("./R/fit_models.R")
source("./R/loo_pit.R")

library("ggplot2")

LOG_RESPONSE <- TRUE
MODELS <- c("M1", "M2", "M3", "M4", "N1", "N2", "N3", "N4")
MODELS <- MODELS[file.exists(file.path("output", sprintf("%s_fit.rds", MODELS)))]

if (!length(MODELS)) stop("no fits in output/; run run_all.R first")

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

  fit <- readRDS(file.path("./output", sprintf("%s_fit.rds", m))) # saved fit
  pit_obj <- anoban_loo_pit(m, anoban_data_for(ANOBAN_MODELS[[m]]$predictor), fit, stan_dir = "./stan") # LOO-PIT

  # LOO-PIT plots via bayesplot
  plots[[m]] <- plot_anoban_loo_pit(pit_obj)
  out <- file.path("./output", sprintf("loo_pit_%s.pdf", m))
  ggsave(out, plots[[m]], width = 15, height = 4.5)
  cat(sprintf("wrote %s\n", out))
}
