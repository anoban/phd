## ── Leave-one-species-out cross-validation (LOGO) ───────────────────
##
## The workflow has three stages:
##
##   1. Fit each model to all data (the full fit).  Its generated
##      quantities give log_lik_g[s, g] = log p(y_g | y_-g, theta_s) for
##      every species g; the species effects are integrated out (see the
##      Stan files).  PSIS on these columns gives the LOGO elpd, its MCSE
##      and the Pareto k of every species.
##   2. The escalation set is the set of species with Pareto k > 0.7 in
##      at least one model.  For every model and every species in it, fit
##      the model without the species (the fold refit) and compute the
##      CV bridge between the full fit and the fold refit.  Where the
##      bridge trust rule flags the plain CV bridge, use the mean-MM CV
##      bridge; the result is the composite `cvbridge_mm`.  Two
##      cross-checks use the same fits: brute force, and bridge sampling
##      of the marginal likelihoods (`mlbridge`).  Only `cvbridge_mm`
##      enters the comparison.
##   3. For each model, the PSIS values with the escalated species
##      replaced by `cvbridge_mm`, and the model comparison.
##
## The fold refit holds the species out by zero weight: it uses the same
## data with the Stan data field `holdout`, so the full fit and the fold
## refit have the same parameters, as the CV bridge requires.
##
## The estimators and diagnostics are the loobench functions in
## R/cv_core.R (loobench METHODS.md sections 3, 6, 7, 8.11 and 10).  The
## column names follow loobench (`elpd_cvbridge`, `mcse_cvbridge`, ...).
##
## Requires R/build_data.R, R/fit_models.R and R/cv_core.R (run_logo.R
## sources them).

options(mc.cores = 8)
options(loo.cores = 8)

LOGO_SETTINGS <- list(chains = 8L, iter_warmup = 1000L, iter_sampling = 1000L, init = 0.3, seed = 20261011L, k_threshold = 0.7, sig_figs = 18) # Aki's code did not have the sig_figs member of this list

## ── Paths ───────────────────────────────────────────────────────────

logo_stage1_path <- function(out_dir, model) file.path(out_dir, "stage1", sprintf("%s.Rds", model))

logo_stage2_path <- function(out_dir, model, species) file.path(out_dir, "stage2", model, sprintf("%s.Rds", species))

logo_tasks_path <- function(out_dir) file.path(out_dir, "tasks.Rds")

## ── Compiled models, with model methods ─────────────────────────────
##
## The CV bridge and `mlbridge` evaluate the log density at draws that
## the sampler did not produce, through the cmdstanr model methods
## (`$log_prob()`, `$unconstrain_draws()`).  `.lb_compile_args()`
## (R/cv_core.R) gives the compile arguments for both cmdstanr versions:
## before 1.0 the model must be compiled with
## `compile_model_methods = TRUE` and `force_recompile = TRUE`; from 1.0
## the methods compile when they are first used.
##
## The Stan files are copied to a temporary directory of the process and
## compiled there.  Otherwise `force_recompile` from several processes
## at the same time would rebuild the one executable in stan/.
##
## The compiled model is kept for the life of the process.
## logo_run_local() compiles before it forks, so the worker processes
## have the compiled methods and do not compile again.
.logo_models <- new.env()

logo_model <- function(model, stan_dir = "stan") {
  spec <- ANOBAN_MODELS[[model]]
  if (is.null(spec)) stop("unknown model: ", model)
  key <- spec$stan
  hit <- get0(key, envir = .logo_models, inherits = FALSE)
  if (!is.null(hit)) return(hit)
  dest <- file.path(tempdir(), "anoban_logo_stan")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  src <- list.files(stan_dir, pattern = "[.](stan|stanfunctions)$", full.names = TRUE)
  file.copy(src, dest, overwrite = TRUE)
  args <- .lb_compile_args(file.path(dest, spec$stan))
  args$include_paths <- dest
  mod <- do.call(cmdstanr::cmdstan_model, args)
  assign(key, mod, envir = .logo_models)
  mod
}

## ── Data, fits and the adapter ──────────────────────────────────────

#' Stan data of one model, with species `holdout` held out (0 = none).
logo_standata <- function(model, bundle, holdout = 0L) {
  spec <- ANOBAN_MODELS[[model]]
  sd <- bundle$standata
  sd$cor_t <- spec$cor_t
  if (!is.na(spec$cor_p)) sd$cor_p <- spec$cor_p
  sd$compute_loo <- 0L
  sd$compute_pit <- 0L
  sd$compute_logo <- 1L
  sd$holdout <- as.integer(holdout)
  sd
}

#' Seed of one fit, different for each model and each held-out species.
logo_seed <- function(model, holdout = 0L, settings = LOGO_SETTINGS) {
  m <- match(model, names(ANOBAN_MODELS))
  as.integer(settings$seed + 10000L * m + as.integer(holdout))
}

#' Fit one model to all data (the full fit, holdout = 0) or without one
#' species (a fold refit).
#'
#' @return a bundle: list(fit, model, holdout, seed, diag, wall)
logo_fit <- function(model, bundle, holdout = 0L, settings = LOGO_SETTINGS, stan_dir = "stan") {
  mod <- logo_model(model, stan_dir)
  seed <- logo_seed(model, holdout, settings)
  t0 <- proc.time()[["elapsed"]]
  # `sig_figs`: digits of the CSV output.  NULL keeps the CmdStan default.
  # A check that compares log_prob() at the draws with log_lik_g needs 18, so that the CSV values are not rounded.
  fit <- mod$sample(data = logo_standata(model, bundle, holdout),
                    chains = settings$chains,
                    parallel_chains = settings$chains,
                    iter_warmup = settings$iter_warmup,
                    iter_sampling = settings$iter_sampling,
                    seed = seed, init = settings$init,
                    sig_figs = settings$sig_figs, # specify sig_figs = 18
                    refresh = 0, show_messages = FALSE)
  wall <- proc.time()[["elapsed"]] - t0
  d <- fit$diagnostic_summary(quiet = TRUE)
  list(fit = fit, model = model, holdout = as.integer(holdout), seed = seed,
       wall = wall,
       diag = list(n_divergent = sum(d$num_divergent),
                   n_max_treedepth = sum(d$num_max_treedepth),
                   min_ebfmi = suppressWarnings(min(d$ebfmi))))
}

#' S x G matrix of log_lik_g, chains stacked (chain 1 first).
logo_log_lik_g <- function(b) {
  x <- unclass(posterior::as_draws_matrix(b$fit$draws("log_lik_g")))
  dimnames(x) <- NULL
  x
}

#' The adapter for the loobench bridge functions: the unconstrained
#' draws of a bundle, and the log density (with Jacobian) at
#' unconstrained rows.  The bodies are those of loobench R/adapter.R
#' (`adapter_base()`, members `unconstrain` and `log_density`).
logo_adapter <- function() {
  list(
    unconstrain = function(bundle, draws = NULL) {
      ud <- if (is.null(draws)) bundle$fit$unconstrain_draws(format = "draws_matrix")
            else bundle$fit$unconstrain_draws(draws = draws, format = "draws_matrix")
      out <- unclass(posterior::as_draws_matrix(ud))
      if (!is.matrix(out)) out <- as.matrix(out)
      dimnames(out) <- list(NULL, colnames(out))
      out
    },
    log_density = function(bundle, upars) {
      upars <- as.matrix(upars)
      apply(upars, 1L, function(u)
        bundle$fit$log_prob(unconstrained_variables = as.numeric(u),
                            jacobian = TRUE))
    })
}

## ── Stage 1: full fits and PSIS-LOGO ────────────────────────────────

#' Fit one model to all data and compute PSIS-LOGO for every species.
#' Writes stage1/<model>.Rds and returns its content.
logo_stage1_model <- function(model, bundle, out_dir, settings = LOGO_SETTINGS, stan_dir = "stan") {
  path <- logo_stage1_path(out_dir, model)
  if (file.exists(path)) return(readRDS(path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  full <- logo_fit(model, bundle, 0L, settings, stan_dir)
  llg <- logo_log_lik_g(full)
  p <- psis_elpd_fold(llg, n_chains = settings$chains)
  out <- list(model = model, species = bundle$species,
              ng = bundle$standata$ng,
              elpd_i = p$elpd, mcse_i = p$mcse, k_psis = p$k,
              elpd = sum(p$elpd), se = sqrt(length(p$elpd)) * stats::sd(p$elpd),
              ## the maximum over the full-fit draws of the log predictive
              ## density of each species: the full-data bound of brute force
              llg_max = apply(llg, 2L, max),
              seed = full$seed, wall = full$wall, diag = full$diag,
              settings = settings)
  saveRDS(out, path)
  out
}

#' The escalation set: the species with Pareto k > threshold in at least
#' one model, and the species in `extra`.  The Pareto k follows the
#' loobench convention (r_eff = 1, right tail), as in psis_elpd_fold().
logo_escalation_set <- function(out_dir, models, k_threshold = 0.7, extra = NULL) {
  s1 <- lapply(models, function(m) {
    p <- logo_stage1_path(out_dir, m)
    if (!file.exists(p)) stop("stage 1 missing for ", m, ": ", p)
    readRDS(p)
  })
  flagged <- unique(unlist(lapply(s1, function(x)
    x$species[x$k_psis > k_threshold])))
  sp <- union(flagged, extra)
  unknown <- setdiff(sp, s1[[1]]$species)
  if (length(unknown)) stop("unknown species: ", paste(unknown, collapse = ", "))
  ## the species order of the data
  s1[[1]]$species[s1[[1]]$species %in% sp]
}

## ── Stage 2: one task = one full fit + a batch of held-out species ──

#' The task list: every model x every escalated species, in batches of
#' `batch` species per task (`Inf` = one task per model).  The tasks of
#' the pruning models are first, because they take the longest; this
#' shortens the wall time when the tasks run in parallel.
logo_plan_tasks <- function(models, species, batch = Inf) {
  slow <- vapply(models, function(m) grepl("pruning", ANOBAN_MODELS[[m]]$stan),
                 logical(1))
  models <- c(models[slow], models[!slow])
  b <- if (is.finite(batch)) max(1L, as.integer(batch)) else length(species)
  tasks <- list()
  for (m in models) {
    chunks <- split(species, ceiling(seq_along(species) / b))
    for (ch in chunks) tasks[[length(tasks) + 1L]] <- list(model = m, species = ch)
  }
  tasks
}

#' One row of stage-2 results: one model, one held-out species.
.logo_fold_row <- function(model, species, g, full, fold, llg_max, settings) {
  ad <- logo_adapter()
  n_ch <- settings$chains
  cfg <- .lb_bridge_config(list())

  ## brute force: the log predictive density of the held-out species at
  ## the fold-refit draws; the full-data bound is its maximum at the
  ## full-fit draws
  ll_fold <- logo_log_lik_g(fold)[, g]
  bf <- .lb_bf_vector_stats(ll_fold, n_ch, y_max = Inf,
                            y_max_full = exp(llg_max[g]))

  ## CV bridge, plain and mean-MM components, from 6 log-density passes
  t0 <- proc.time()[["elapsed"]]
  qs <- .lb_bridge_qs(ad, full, fold, "logo", model)
  plain <- .lb_bridge_arm_plain(qs, cfg)
  mm <- .lb_bridge_arm_mm(qs, cfg)
  suspect <- .lb_bridge_suspect(plain, cfg)
  ## loobench METHODS.md section 11.1: the mean-MM component on a suspect
  ## fold (a plain value that is not finite always takes it), else the
  ## plain component
  cmp <- if (isTRUE(suspect)) mm else plain
  wall_cvb <- proc.time()[["elapsed"]] - t0

  ## bridge sampling of the two log marginal likelihoods, each against
  ## its own fitted normal proposal; the proposal seeds differ per fold
  t0 <- proc.time()[["elapsed"]]
  s0 <- logo_seed(model, g, settings)
  lm_full <- .lb_normal_bridge_logml(full, ad, n_ch, s0 + 2121L)
  lm_fold <- .lb_normal_bridge_logml(fold, ad, n_ch, s0 + 4242L)
  elpd_ml <- lm_full$logml - lm_fold$logml
  m_ml <- c(lm_full$mcse_logml, lm_fold$mcse_logml)
  wall_ml <- proc.time()[["elapsed"]] - t0

  num <- function(x) suppressWarnings(as.numeric(x %||% NA_real_))[1L]
  data.frame(
    model = model, species = species, g = g,
    ## composite (enters the comparison)
    elpd_cvbridge_mm = num(cmp$elpd), mcse_cvbridge_mm = num(cmp$mcse),
    mm_cvbridge_mm_applied = isTRUE(suspect),
    ## components
    elpd_cvbridge = num(plain$elpd), mcse_cvbridge = num(plain$mcse),
    gap_cvbridge = num(plain$gap), niter_cvbridge = num(plain$n_iter),
    error_cvbridge = as.character(plain$error %||% NA_character_),
    elpd_cvbridge_mmall = num(mm$elpd), mcse_cvbridge_mmall = num(mm$mcse),
    gap_cvbridge_mmall = num(mm$gap), niter_cvbridge_mmall = num(mm$n_iter),
    mm_delta_norm = num(mm$delta_norm),
    ## cross-check: bridge sampling of the marginal likelihoods
    elpd_mlbridge = num(elpd_ml),
    mcse_mlbridge = if (all(is.finite(m_ml))) sqrt(sum(m_ml^2)) else NA_real_,
    logml_full = num(lm_full$logml), logml_fold = num(lm_fold$logml),
    ## cross-check: brute force
    elpd_brute = num(bf$elpd), mcse_brute = num(bf$mcse), k_brute = num(bf$k),
    mcse_truncgpd_brute = num(bf$mcse_truncgpd),
    bias_jensen_truncgpd_brute = num(bf$bias_truncgpd),
    mcse_truncgpd_brute_hybrid = num(bf$mcse_truncgpd_hybrid),
    bias_jensen_truncgpd_brute_hybrid = num(bf$bias_truncgpd_hybrid),
    k_truncgpd_brute = num(bf$k_truncgpd),
    ## fits
    seed_full = full$seed, seed_fold = fold$seed,
    n_divergent_full = full$diag$n_divergent,
    n_divergent_fold = fold$diag$n_divergent,
    n_max_treedepth_fold = fold$diag$n_max_treedepth,
    wall_fit_full = full$wall, wall_fit_fold = fold$wall,
    wall_cvbridge = wall_cvb, wall_mlbridge = wall_ml,
    host = Sys.info()[["nodename"]],
    stringsAsFactors = FALSE)
}

#' Run one task: the full fit of `model` once, then for each species in
#' `species` that has no result file yet, the fold refit and the
#' estimates, written to stage2/<model>/<species>.Rds.
logo_task <- function(task, bundle, out_dir, settings = LOGO_SETTINGS, stan_dir = "stan") {
  model <- task$model
  todo <- task$species[!file.exists(
    logo_stage2_path(out_dir, model, task$species))]
  if (!length(todo)) return(invisible(character(0)))
  s1 <- readRDS(logo_stage1_path(out_dir, model))
  full <- logo_fit(model, bundle, 0L, settings, stan_dir)
  for (sp in todo) {
    g <- match(sp, bundle$species)
    fold <- logo_fit(model, bundle, g, settings, stan_dir)
    row <- .logo_fold_row(model, sp, g, full, fold, s1$llg_max, settings)
    row$elpd_psis <- s1$elpd_i[g]
    row$k_psis <- s1$k_psis[g]
    row$ng <- bundle$standata$ng[g]
    path <- logo_stage2_path(out_dir, model, sp)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(row, path)
    cat(sprintf("%s %-26s cvbridge_mm %9.2f (mcse %.3f, %s) | psis %9.2f | brute %9.2f | mlbridge %9.2f\n",
                model, sp, row$elpd_cvbridge_mm, row$mcse_cvbridge_mm,
                if (row$mm_cvbridge_mm_applied) "MM" else "plain",
                row$elpd_psis, row$elpd_brute, row$elpd_mlbridge))
    rm(fold)
    gc(verbose = FALSE)
  }
  invisible(todo)
}

#' Run tasks in parallel on this machine with parallel::mclapply.  It
#' uses forked processes, so it does not run in parallel on Windows.
#' Each task runs `settings$chains` chains in parallel, so
#' workers = floor(cores / chains) uses the machine.  The models are
#' compiled before the fork, so the workers have the compiled model
#' methods.
logo_run_local <- function(tasks, bundle_for, out_dir, workers = 1L, settings = LOGO_SETTINGS, stan_dir = "stan") {
  for (m in unique(vapply(tasks, `[[`, "", "model"))) logo_model(m, stan_dir)
  res <- parallel::mclapply(tasks, function(t) {
    tryCatch({
      logo_task(t, bundle_for(t$model), out_dir, settings, stan_dir)
      NA_character_
    }, error = function(e) sprintf("%s [%s]: %s", t$model,
                                   paste(t$species, collapse = ","),
                                   conditionMessage(e)))
  }, mc.cores = max(1L, as.integer(workers)), mc.preschedule = FALSE)
  err <- unlist(res)
  err <- err[!is.na(err)]
  if (length(err)) warning(paste(c("failed tasks:", err), collapse = "\n  "),
                           call. = FALSE)
  invisible(err)
}

## ── Stage 3: assemble and compare ───────────────────────────────────

#' The LOGO elpd of each species and model: PSIS, with the escalated
#' species replaced by `cvbridge_mm`.  Stops if a stage-2 result is
#' missing.
#'
#' @return list(E = G x models matrix, stage2 = data frame of all rows)
logo_assemble <- function(out_dir, models, species) {
  s1 <- lapply(models, function(m) readRDS(logo_stage1_path(out_dir, m)))
  names(s1) <- models
  miss <- character(0)
  rows <- list()
  for (m in models) for (sp in species) {
    p <- logo_stage2_path(out_dir, m, sp)
    if (file.exists(p)) rows[[length(rows) + 1L]] <- readRDS(p)
    else miss <- c(miss, sprintf("%s/%s", m, sp))
  }
  if (length(miss)) stop("stage-2 results missing: ", paste(miss, collapse = ", "))
  st2 <- do.call(rbind, rows)
  E <- sapply(models, function(m) {
    e <- s1[[m]]$elpd_i
    r <- st2[st2$model == m, ]
    e[r$g] <- r$elpd_cvbridge_mm
    e
  })
  if (any(!is.finite(E))) stop("non-finite per-species elpd after assembly")
  rownames(E) <- s1[[1]]$species
  list(E = E, stage2 = st2)
}

#' Model comparison in the layout of loo_compare(), and pairwise
#' differences.  se_diff = sqrt(G) sd(e_a - e_b) over the G species.
logo_compare <- function(E, pairs = NULL) {
  G <- nrow(E)
  tot <- colSums(E)
  ord <- order(tot, decreasing = TRUE)
  best <- ord[1]
  tab <- data.frame(
    model = colnames(E)[ord], elpd_logo = tot[ord],
    se = sqrt(G) * apply(E, 2, stats::sd)[ord],
    elpd_diff = tot[ord] - tot[best],
    se_diff = vapply(ord, function(i)
      if (i == best) 0 else sqrt(G) * stats::sd(E[, i] - E[, best]), 0),
    row.names = NULL)
  out <- list(table = tab)
  if (!is.null(pairs)) {
    out$pairs <- do.call(rbind, lapply(pairs, function(p) {
      d <- E[, p[1]] - E[, p[2]]
      data.frame(pair = sprintf("%s-%s", p[1], p[2]),
                 diff = sum(d), se_diff = sqrt(G) * stats::sd(d),
                 ratio = sum(d) / (sqrt(G) * stats::sd(d)),
                 stringsAsFactors = FALSE)
    }))
  }
  out
}

#' Diagnostics of the escalated folds, per model.
#'   * CV bridge: the MCSE of the composite, the number of folds that
#'     use the mean-MM component, and the number of folds where the
#'     composite still fails the bridge trust rule (not finite, OR
#'     mcse >= 0.3, OR cross-density gap > 50).
#'   * brute force (cross-check): over the escalated folds, the total of
#'     the truncated-GPD MCSE and Jensen bias at the hybrid bound, and
#'     rmse = sqrt(bias^2 + mcse^2).  Brute force is trusted when the
#'     rmse is below 1 nat.
#'   * the largest |cvbridge_mm - mlbridge|.
logo_diagnostics <- function(st2) {
  do.call(rbind, lapply(split(st2, st2$model), function(r) {
    gap_used <- ifelse(r$mm_cvbridge_mm_applied, r$gap_cvbridge_mmall,
                       r$gap_cvbridge)
    still <- !is.finite(r$elpd_cvbridge_mm) | r$mcse_cvbridge_mm >= 0.3 |
      gap_used > 50
    mt <- sqrt(sum(r$mcse_truncgpd_brute_hybrid^2, na.rm = TRUE))
    bi <- sum(r$bias_jensen_truncgpd_brute_hybrid, na.rm = TRUE)
    data.frame(model = r$model[1], n_folds = nrow(r),
               n_mm = sum(r$mm_cvbridge_mm_applied),
               n_suspect_after = sum(still),
               max_mcse_cvbridge_mm = max(r$mcse_cvbridge_mm),
               brute_rmse = sqrt(bi^2 + mt^2),
               max_abs_cvb_minus_ml = max(abs(r$elpd_cvbridge_mm - r$elpd_mlbridge)),
               stringsAsFactors = FALSE)
  }))
}
