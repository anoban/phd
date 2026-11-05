## ── Comparing models fitted on different response scales ────────────
##
## Code for comparing a log-response fit against a raw-response one,
## because they are densities for different random variables.

suppressPackageStartupMessages({
  library(loo)
})

#' Compare models on the original response scale.
#'
#' log-response and raw-response fits model different random variables,
#' so their elpd values are not comparable as they stand.  With
#' z = log y elementwise (a bijection),
#'
#'   p_Y(y_i | y_-i) = p_Z(z_i | z_-i) * prod_d (1 / y_id)
#'   => elpd_i^Y = elpd_i^Z - sum_d log y_id
#'
#' The correction is applied pointwise so loo_compare() still sees
#' matched folds and reports a paired SE.  Within one scale it is a
#' no-op, so it is always safe to call.
#'
#' @param loo_list named list of loo objects
#' @param Y the original (untransformed) response matrix
#' @param log_response TRUE if the fits used log responses
loo_compare_original_scale <- function(loo_list, Y, log_response) {
  if (log_response) {
    jacobian <- rowSums(log(as.matrix(Y)))
    loo_list <- lapply(loo_list, function(lo) {
      stopifnot(nrow(lo$pointwise) == length(jacobian))
      lo$pointwise[, "elpd_loo"] <- lo$pointwise[, "elpd_loo"] - jacobian
      n <- nrow(lo$pointwise)
      lo$estimates["elpd_loo", "Estimate"] <- sum(lo$pointwise[, "elpd_loo"])
      lo$estimates["elpd_loo", "SE"] <-
        sqrt(n * stats::var(lo$pointwise[, "elpd_loo"]))
      lo
    })
  }
  loo::loo_compare(loo_list)
}
