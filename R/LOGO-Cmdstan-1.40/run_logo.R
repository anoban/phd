## ── Leave-one-species-out CV (LOGO) for the eight models ────────────
##
## Stages (see R/logo.R):
##   1     full fit of each model; PSIS-LOGO of each species
##   plan  the escalation set (Pareto k > 0.7 in at least one model) and
##         the task list
##   2     a fold refit for each model and each escalated species; the CV
##         bridge (with the mean-MM component), brute force and mlbridge
##   3     assemble and compare
##   all   1, plan, 2 and 3 on this machine
##
## On one machine:
##   Rscript run_logo.R --stage all --workers 8
##
## On a SLURM cluster (see slurm/run_logo_array.sbatch):
##   sbatch --array=1-8 --export=ALL,STAGE=1 slurm/run_logo_array.sbatch
##   Rscript run_logo.R --stage plan --batch 1
##   sbatch --array=1-<n> --export=ALL,STAGE=2 slurm/run_logo_array.sbatch
##   Rscript run_logo.R --stage 3
##
## Options:
##   --stage 1|plan|2|3|all
##   --models M1,M2,...     default: all eight
##   --workers n            parallel tasks on this machine (default 1);
##                          each task runs 4 chains in parallel
##   --task i               run only task i (stage 1: model i; stage 2:
##                          task i of tasks.Rds); for SLURM arrays
##   --batch k              species per stage-2 task (default: all species
##                          of a model in one task, which then needs only
##                          one full fit)
##   --extra sp1,sp2        add species to the escalation set
##   --species sp1,sp2      use exactly these species instead of the
##                          escalation set (a rerun, or a smoke test)
##   --out dir              default output/logo (output/logo_quick with --quick)
##   --quick                short chains, for a smoke test only
##   --dry-run              print the plan and stop
##
## A rerun skips every model and species that already has a result file.

library("cmdstanr")
library("posterior")
library("loo")

source("./R/build_data.R")
source("./R/fit_models.R")
source("./R/cv_core.R")
source("./R/logo.R")

options(mc.cores = 8)
options(loo.cores = 8)

split_csv <- function(x) if (is.null(x) || !nzchar(x)) NULL else strsplit(x, ",")[[1]]

STAGE <- "all"
if (!STAGE %in% c("1", "plan", "2", "3", "all")) stop("stage must be one of 1, plan, 2, 3, all")

MODELS <- c("M1", "M2", "M3", "M4", "N1", "N2", "N3", "N4")
bad <- setdiff(MODELS, names(ANOBAN_MODELS))
if (length(bad)) stop("unknown model(s): ", paste(bad, collapse = ", "))

WORKERS <- 32 / 4
TASK <- NA
BATCH <- Inf
EXTRA <- NULL
SPECIES <- NULL
OUT <- "./output/logo"
SETTINGS <- LOGO_SETTINGS

## Data: one bundle per right-hand side (M* use `state`, N* use `1`).
bundle_for <- local({
  cache <- list()
  function(model) {
    rhs <- ANOBAN_MODELS[[model]]$predictor
    if (is.null(cache[[rhs]]))
      cache[[rhs]] <<- build_anoban_data("data/stan.csv", "data/FRED4_1301.tre", log_response = TRUE, predictor = rhs)
    cache[[rhs]]
  }
})

dir.create(OUT, recursive = TRUE, showWarnings = TRUE)
cat(sprintf("LOGO | stage %s | out %s | %d chains x (%d + %d)\n", STAGE, OUT, SETTINGS$chains, SETTINGS$iter_warmup, SETTINGS$iter_sampling))

## ── stage 1 ─────────────────────────────────────────────────────────
run_stage1 <- function() {
  todo <- if (is.na(TASK)) MODELS else MODELS[TASK]
  if (anyNA(todo)) stop("--task ", TASK, " is beyond the ", length(MODELS), " models")
  cat(sprintf("stage 1: %s\n", paste(todo, collapse = ", ")))

  for (m in todo) logo_model(m)                  # compile before forking
  res <- parallel::mclapply(todo, function(m) tryCatch({
    s <- logo_stage1_model(m, bundle_for(m), OUT, SETTINGS) # logo_stage1_model() will call saveRDS() to serialize the fitted model objects
    sprintf("%s: elpd_logo %.1f (se %.1f) | k > 0.7: %d | max k %.2f | %.0f s | divergent %d",
            m, s$elpd, s$se, sum(s$k_psis > SETTINGS$k_threshold), max(s$k_psis), s$wall, s$diag$n_divergent)
  }, error = function(e) sprintf("%s: FAILED: %s", m, conditionMessage(e))),
  mc.cores = max(1L, WORKERS), mc.preschedule = FALSE)

  cat(paste0("  ", unlist(res), collapse = "\n"), "\n")
}

## ── plan ────────────────────────────────────────────────────────────
run_plan <- function() {
  have <- file.exists(vapply(MODELS, function(m) logo_stage1_path(OUT, m), ""))

  sp <- if (!is.null(SPECIES)) SPECIES
        else logo_escalation_set(OUT, MODELS, SETTINGS$k_threshold, EXTRA)
  unknown <- setdiff(sp, bundle_for(MODELS[1])$species)
  if (length(unknown)) stop("unknown species: ", paste(unknown, collapse = ", "))
  tasks <- logo_plan_tasks(MODELS, sp, BATCH)
  cat(sprintf("escalation set (%d species): %s\n", length(sp), paste(sp, collapse = ", ")))
  cat(sprintf("%d stage-2 tasks:\n", length(tasks)))
  for (i in seq_along(tasks))
    cat(sprintf("  %3d  %-3s %s\n", i, tasks[[i]]$model,
                paste(tasks[[i]]$species, collapse = ", ")))
  saveRDS(list(tasks = tasks, species = sp, models = MODELS, settings = SETTINGS), logo_tasks_path(OUT))
}

## ── stage 2 ─────────────────────────────────────────────────────────
run_stage2 <- function() {
  p <- logo_tasks_path(OUT)
  if (!file.exists(p)) stop("no task list; run --stage plan first")
  tl <- readRDS(p)
  if (!identical(tl$settings, SETTINGS))
    stop("the task list was planned with other sampler settings (--quick?); plan again or use the same settings")
  tasks <- if (is.na(TASK)) tl$tasks else tl$tasks[TASK]
  if (!length(tasks) || is.null(tasks[[1]])) stop("no task ", TASK)
  cat(sprintf("stage 2: %d task(s), %d worker(s)\n", length(tasks), WORKERS))

  if (length(tasks) == 1L) {
    logo_task(tasks[[1]], bundle_for(tasks[[1]]$model), OUT, SETTINGS)
  } else {
    err <- logo_run_local(tasks, bundle_for, OUT, WORKERS, SETTINGS)
    if (length(err)) quit(status = 1)
  }
}

## ── stage 3 ─────────────────────────────────────────────────────────
PAIRS <- list(c("M3", "M4"), c("M1", "M2"), c("N3", "N4"), c("N1", "N2"), # both phylogenetic and non-phylogenetic random effects vs only non-phylogenetic random effects
              c("M3", "M1"), c("M4", "M2"), c("N3", "N1"), c("N4", "N2"), # cross trait interactions vs no cross trait interactions
              c("M3", "N3"), c("M1", "N1"), c("M4", "N4"), c("M2", "N2")) # fixed effect vs no fixed effect

run_stage3 <- function() {
  tl <- readRDS(logo_tasks_path(OUT))
  a <- logo_assemble(OUT, tl$models, tl$species)
  pairs <- Filter(function(p) all(p %in% tl$models), PAIRS)
  cmp <- logo_compare(a$E, pairs)
  dg <- logo_diagnostics(a$stage2)
  cat(sprintf("\nLOGO model comparison, %d species (%d escalated: CV bridge + MM)\n\n", nrow(a$E), length(tl$species)))
  print(format(cmp$table, digits = 1, nsmall = 1), row.names = FALSE)
  if (!is.null(cmp$pairs)) {
    cat("\npairwise (first minus second):\n")
    print(format(cmp$pairs, digits = 1, nsmall = 1), row.names = FALSE)
  }
  cat("\ndiagnostics of the escalated folds:\n")
  print(format(dg, digits = 2), row.names = FALSE)
  saveRDS(list(E = a$E, stage2 = a$stage2, compare = cmp, diagnostics = dg, settings = tl$settings), file.path(OUT, "logo_results.Rds"))
  utils::write.csv(cmp$table, file.path(OUT, "logo_compare.csv"), row.names = FALSE)
  if (!is.null(cmp$pairs))
    utils::write.csv(cmp$pairs, file.path(OUT, "logo_pairs.csv"), row.names = FALSE)
  utils::write.csv(a$stage2, file.path(OUT, "logo_escalated_folds.csv"), row.names = FALSE)
  cat(sprintf("\nwritten to %s/\n", OUT))
}

switch(STAGE,
  "1" = run_stage1(),
  "plan" = run_plan(),
  "2" = run_stage2(),
  "3" = run_stage3(),
  "all" = { run_stage1(); run_plan(); run_stage2(); run_stage3() }
)
