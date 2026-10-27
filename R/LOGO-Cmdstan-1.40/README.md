# anoban: trivariate species models with the species effects integrated out

Stan and R code for the models M1–M4 and N1–N4.  The code integrates
the species-level effects out analytically instead of sampling them.
The code combines brms-generated code, code from a related research
project on LOO computation, new code by Claude, and edits by Aki.
Claude generated this README; Aki checked and edited it and all code
comments.  The LOGO code (`run_logo.R`, `R/logo.R`, `R/cv_core.R`,
`slurm/`) and its README section were added in October 2026.

The integration has two benefits.

- Speed.  brms samples the species-level effects: 3903 of them in
  M2/M4 (3 responses × 1301 species) and 7806 in M1/M3.  Here they are
  integrated out, and 27–36 free parameters remain.  Measured by
  `run_all.R` on one cluster node, 4 chains × (1000 warmup + 1000
  sampling), log responses.  Wall times vary by a factor of about 2
  across machines; the ESS values do not.

  | | M1 | M2 | M3 | M4 |
  |---|---:|---:|---:|---:|
  | wall (s) | 659 | 32 | 1016 | 60 |
  | min ess_bulk | 2714 | 3207 | 3656 | 4072 |

  None of the fits has divergences or treedepth saturation.  For
  comparison, the brms M2, measured elsewhere, took 691 s to sample and
  1091 s more for `loo()` (raw response scale, so the timings are not
  directly comparable).  For M1/M3 the difference is larger than speed:
  a version that sampled the 7806 phylogenetic effects hit the maximum
  treedepth on every transition and gave `ess_bulk = 2`.
- Reliable LOO.  In the brms parameterisation, 732 of the 1301 species
  have exactly one observation, so leaving out an observation changes
  the posterior of its own species effect a lot.  Measured on M2: 787
  of 5218 folds have Pareto k > 0.7, and `p_loo = 2942`.  With the
  species effects integrated out, the same model gives 0–2 folds with
  Pareto k > 0.7 and `p_loo` of 53–62.

The models are unchanged: the same data model, the same priors and the
same posterior for the population-level parameters.  Only the
computation is different.

---

## Quick start

### Requirements

Install these before you run the code.  The second column gives the
versions that the code was tested with.

| software | tested with | needed for |
|---|---|---|
| R | 4.6.1 | all scripts |
| CmdStan | 2.39.0 | all scripts |
| `cmdstanr` | 0.9.0 | all scripts |
| `posterior` | 1.7.1 | all scripts |
| `loo` | 2.10.1 | all scripts |
| `ape` | 5.8.1 | all scripts (reads the tree) |
| `matrixStats` | 1.5.0 | `run_logo.R` |
| `bayesplot` | 1.15.0 | `run_loo_pit.R` |
| `ggplot2` | 4.0.3 | `run_loo_pit.R` |
| `patchwork` | 1.3.2 | `run_loo_pit.R` |

brms is not needed.  To install the R packages and CmdStan:

```r
install.packages(c("posterior", "loo", "ape", "matrixStats",
                   "bayesplot", "ggplot2", "patchwork"))
install.packages("cmdstanr",
                 repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()        # installs the newest CmdStan release
cmdstanr::cmdstan_version()        # check the CmdStan version
```

`cmdstanr` is not on CRAN; the command above installs it from the Stan
R-universe.  `install_cmdstan()` needs a C++ toolchain; see
`cmdstanr::check_cmdstan_toolchain()`.  The LOGO code also supports
cmdstanr 1.0 and later, but it was tested only with 0.9.  `run_logo.R`
runs folds in parallel with `parallel::mclapply()`, which uses forked
processes.  On Windows, use `--workers 1`.

### Run LOO and LOO-PIT

You can also run the scripts in an interactive R session instead of
with `Rscript`.

```r
Rscript run_all.R                 # fit M1-M4 and N1-N4 on log responses, compute LOO
Rscript run_all.R --raw           # original (untransformed) response scale
Rscript run_all.R --models M3     # one model only
Rscript run_all.R --quick         # short chains, for a smoke test
```

The results are in `output/`: one `_fit.rds` and one `_loo.rds` per
model, and `loo_compare.txt`.  `print(loo)` prints the elpd_loo,
`p_loo` and the Pareto k table of each model when the model finishes.

LOO-PIT calibration plots, after `run_all.R`:

```r
Rscript run_loo_pit.R             # all fitted models
Rscript run_loo_pit.R --models M3
```

Leave-one-species-out CV: see [LOGO](#leave-one-species-out-cv-logo).

Run the checks of the derivations (see [below](#what-validater-checks)):

```r
Rscript R/validate.R              # all checks, compiles Stan (about 5 min)
Rscript R/validate.R --no-stan    # R-level checks only, seconds
```

---

## Response scale

`run_all.R` uses log responses by default.  Use `--raw` for the
original scale.

`F00727` has skewness 4.3 and a maximum of 840, against a median of 20.
A normal model for it on the raw scale is badly misspecified.  Three
diagnostics show this:

| | raw scale | log responses |
|---|---:|---:|
| elpd_loo (M3, both on the original `y` scale) | −20150 | −11759 |
| `p_loo` (M3, against 36 free parameters) | 114 | 62 |
| folds with Pareto k > 0.7 (M3, of 5218) | 4 | 2 |
| largest LOO-PIT departure from uniform, F00727 | 0.169 | 0.053 |

The fold counts are small integers and change by one or two between
runs; the other rows are stable.

The log models are about 8400 nats better, which is 18 paired standard
errors.  `loo_compare_original_scale()` applies the change of variables
to each observation,

```
elpd_i^Y = elpd_i^Z − Σ_d log y_id
```

so both scales give densities for the same random variable `y`, and
`loo_compare()` still reports a paired SE.

All other code in this package works the same on either scale.

---

## The models

| | this package | brms formula |
|---|---|---|
| M1 | `anoban_pruning.stan`, `cor_t=0, cor_p=0` | `~ state + (1 \| gr(binominal, cov = corrmat)) + (1 \| taxa)` |
| M2 | `anoban_marginal.stan`, `cor_t=0` | `~ state + (1 \| taxa)` |
| M3 | `anoban_pruning.stan`, `cor_t=1, cor_p=1` | `~ state + (1\|p\| gr(binominal, cov = corrmat)) + (1\|q\| taxa)` |
| M4 | `anoban_marginal.stan`, `cor_t=1` | `~ state + (1\|q\| taxa)` |

N1–N4 are M1–M4 with `~ 1` instead of `~ state`.  M1, M3, N1 and N3
are the pruning models (`anoban_pruning.stan`, with the phylogenetic
effect); M2, M4, N2 and N4 are the non-phylogenetic models
(`anoban_marginal.stan`).  All models use
`mvbind(F00727, F00679, F00709)` and `set_rescor(TRUE)`.
`cor_t` / `cor_p` = 1 means that the effect is correlated across
responses, i.e. brms `(1|q|...)` instead of `(1|...)`.

The priors are the brms 2.23.1 defaults, computed on the response scale
of the fit: `student_t(3, median(y), max(2.5, mad(y)))` for the
intercepts, half-`student_t(3, 0, max(2.5, mad(y)))` for `sigma` and
the group SDs, `lkj(1)` for all correlation matrices, and flat for
`b`.  On the raw scale this gives the brms values exactly (locations
`20.5, 0.5, 0.3`, scales `17.1, 2.5, 2.5`).

The code does not build `corrmat`.  The pruning model uses the tree
directly, so there is no 1301 × 1301 matrix.

---

## How it works

Both species effects are integrated out in closed form, without
approximation.

*The iid species effect and the residuals* (all models).  `state` is
constant within species, and there are no observation-level
predictors.  So for species `g` with `n` observations the likelihood
splits exactly:

```
log p(y_g) = N_D(ȳ_g | μ_g, Σ_r/n + Σ_t)                    between
           − ½[(n−1) log|Σ_r| + tr(Σ_r⁻¹ S_g) + D(n−1) log 2π]   within
           − (D/2) log n
```

Summed over species, the within part depends on the data only through
the pooled scatter `W = Σ_g S_g` and `N − G`, which R computes once.
The likelihood therefore has 1301 three-dimensional terms instead of
5218 observation terms.

*The phylogenetic effect* (M1/M3, N1/N3) by Felsenstein pruning.
`A = vcv.phylo(tree, corr = TRUE)` is the covariance of Brownian motion
that starts at 0 at the root, with branch lengths divided by the tree
height.  So `A` is the marginal covariance of a Markov process on the
tree.  Pruning uses this: messages go up the tree, gain `t·Σ_p` on each
edge and are multiplied at each node.  The integration costs
O(#nodes · D³), which is 3 × 3 algebra at 2377 nodes, instead of a
dense 3903 × 3903 Cholesky factorization per gradient (about 2 × 10¹⁰
flops).  Polytomies need no special handling.

*The LOO factors.*  `log_lik[i] = log p(y_i | y_{g∖i}, θ)` is the
density of observation i given the other observations of its species,
with the species effects integrated out.  This part has not been used
elsewhere; it will be part of a paper to be put on arXiv in the fall.
A downward (root-to-tips) pruning pass gives `p(u_g | y_{−g}, θ)` for
all 1301 species in one pass, so the cost is O(1) per observation.

---

## LOO-PIT: is the model calibrated?

`elpd_loo` compares models with each other.  It does not tell whether
a model is good.  LOO-PIT does: each observation is mapped through its
own leave-one-out predictive CDF.  For a calibrated model the values
are uniform on (0, 1), and the shape of the departure shows how the
model is wrong.

`run_loo_pit.R` makes one ECDF panel per response with
`bayesplot::ppc_loo_pit_ecdf(method = "correlated", test = "POT",
plot_diff = TRUE)`.  The `"correlated"` method takes into account that
the LOO-PIT values are not independent (Tesso and Vehtari, 2026).  All
panels have the same y-axis, so you can compare the magnitudes across
responses.

The PIT values are computed analytically.  The LOO predictive
distribution here is a finite mixture of normal distributions over the
draws,

```
PIT_id = Σ_s w_s · Φ((y_id − m_sid) / σ_sid)
```

with `w_s` the PSIS weights.  Stan gives the per-draw `Φ(·)` term and R
gives the weights.  So no `y_rep` is simulated, and the diagnostic has
no posterior predictive Monte Carlo error.

To read the plots: the line is the ECDF minus its uniform expectation,
so a calibrated model gives a flat line near zero, inside the band.  An
S-shape (first up, then down) means that the PIT values are too
concentrated near 0.5: the predictive distribution is too wide.  On the
log scale all three responses still show a mild S-shape, so the model
is much better calibrated than on the raw scale, but not well
calibrated.  Even on the log scale there are some outliers (the red
dots near 0 and 1), so the tails of the predictive distribution are
too thin.

---

## Leave-one-species-out CV (LOGO)

`run_logo.R` estimates the elpd for predicting each species when the
model is fitted without it.  With the species effects integrated out,
each posterior draw gives `log p(y_g | y_{−g}, θ)` of every species in
closed form.  The generated quantities give it as `log_lik_g[G]`; in the
pruning model it uses the downward pass.

The workflow has three stages:

1. Fit each model to all data.  PSIS-LOGO (PSIS on the columns of
   `log_lik_g`) gives the elpd, its MCSE and the Pareto k of each
   species.
2. The escalation set is the set of species with Pareto k > 0.7 in at
   least one model.  For these species PSIS is not reliable: a species
   with up to 242 observations is a large change to the posterior.  In
   our runs 5–9 species per model had Pareto k > 0.7, with a maximum of
   2.7–3.9, and the escalation set had 12 species.  For every model and
   every species in
   the escalation set, the code fits the model without the species
   (the fold refit) and computes the CV bridge between the full fit and
   the fold refit.  Where the bridge trust rule flags the plain CV
   bridge (the value is not finite, OR its MCSE ≥ 0.3, OR the
   cross-density gap > 50 nats), it uses the mean-MM CV bridge, which
   moves the full-fit draws to the mean of the fold-refit draws.  This
   composite is `cvbridge_mm`.  Two cross-checks use the same fits:
   brute force, and bridge sampling of the marginal likelihoods
   (`mlbridge`).  Only `cvbridge_mm` enters the model comparison.
3. For each model, the PSIS-LOGO values with the escalated species
   replaced by `cvbridge_mm`; then the model comparison.

The fold refit holds the species out by zero weight (the Stan data
field `holdout`): the data and the parameters stay the same, as the CV
bridge requires.

The estimators and their diagnostics are functions of the R package
loobench (to be published in a few weeks), copied unchanged into
`R/cv_core.R`.  Their definitions are in loobench's `METHODS.md`.

On one machine (Linux or macOS; the parallel tasks use forked
processes):

```r
Rscript run_logo.R --stage all --workers 8
```

Each task runs 4 chains in parallel, so `--workers` = (number of
cores) / 4 uses the machine.  On a SLURM cluster, one array element per
task:

```sh
cp slurm/site_setup.sh.example slurm/site_setup.sh   # edit for your cluster
sbatch --array=1-8 --export=ALL,STAGE=1 slurm/run_logo_array.sbatch
Rscript run_logo.R --stage plan --batch 1     # prints the number of tasks n
sbatch --array=1-<n> --export=ALL,STAGE=2 slurm/run_logo_array.sbatch
Rscript run_logo.R --stage 3
```

A rerun skips every model and species that already has a result file,
so you can resubmit a failed task.  The header of `run_logo.R` lists
the options.

The results are in `output/logo/`: `stage1/<model>.rds`,
`stage2/<model>/<species>.rds`, `logo_compare.csv` (the comparison
table), `logo_pairs.csv` (pairwise differences),
`logo_escalated_folds.csv` (all estimates and diagnostics of the
escalated folds) and `logo_results.rds` (all of these).

Cost, with 4 chains × (1000 + 1000) iterations: a fit of a pruning
model takes about 15 minutes, and a fit of a non-phylogenetic model
about 30 seconds (Aalto Triton, one task per job; CSC Roihu, 8 fits
in parallel).  Each escalated fold of a pruning model needs one fold
refit and about 5–6 minutes for the bridges.  When many tasks share
one node, each fit can be much slower: with 24 tasks (96 chains) in
one 96-core job on Roihu, a fit took 2.5–3.5 times as long, and the
whole run (8 models, 12 escalated species, `--workers 24 --batch 1`)
took 4 h 16 min.  So measure with a few tasks before you choose the
number of workers.

---

## Assumptions of this code

`build_anoban_data()` checks the first two and stops with an error if
they do not hold.  The third is a modelling decision.

1. No observation-level predictors.  The sufficient-statistic form needs
   a design matrix that is constant within species, and `state` is.  If
   you add a covariate that varies within a species, this code does not
   apply and the derivation must be done again.  You then get an error,
   not a wrong answer.
2. `binominal` and `taxa` are the same grouping factor.  They are
   identical in `stanforum.csv`, so M1/M3 divide one species effect into
   a phylogenetically correlated part and an iid part.  If `taxa` were a
   coarser level (genus, family), the block structure would change.
3. The tree is ultrametric, which `FRED4_1301.tre` is (root-to-tip
   height 400.7877 for every tip).  The fixed-root scaling needs this.

`state` is species-constant and very unbalanced: 1098 of 1301 species
are `AM`, and `ErM` has 9, `NM` 16 and `EcMAM` 20.  The coefficients of
the rare states rest on little data.

---

## Files

```
run_all.R                          driver: build -> fit -> LOO -> compare
run_loo_pit.R                      LOO-PIT calibration plots
run_logo.R                         driver: leave-one-species-out CV (LOGO)
R/build_data.R                     csv + tree -> Stan data (+ the checks)
R/fit_models.R                     model table, compile and sample
R/loo_compare_scales.R             comparison of log- and raw-response fits
R/loo_pit.R                        analytic LOO-PIT and the ECDF plots
R/logo.R                           LOGO: fits, estimators, comparison
R/cv_core.R                        CV estimators, copied from loobench
R/validate.R                       checks of the derivations
stan/anoban_functions.stanfunctions  block_marginal, loo_factor_vec, loo_pit_vec
stan/anoban_marginal.stan          M2 / M4 / N2 / N4
stan/anoban_pruning.stan           M1 / M3 / N1 / N3
slurm/run_logo_array.sbatch        SLURM array script for run_logo.R
slurm/site_setup.sh.example        example cluster setup (Aalto Triton)
data/stanforum.csv, data/FRED4_1301.tre
```

`run_loo_pit.R` reads the saved fits and runs only the generated
quantities again (`generate_quantities()`, no new sampling), so the main
fit does not have to store the per-draw PIT terms (`N × D` numbers per
draw).

## What `validate.R` checks

The likelihoods are derived, so `validate.R` checks them against
direct computations:

| check | against |
|---|---|
| compound-symmetry formula | a dense 3n_g multivariate normal, every species |
| pooled-W model-block form | the sum of per-species terms |
| O(1) rank-one LOO factor | O(n) recomputation, group sizes 1 / 2 / 5 / 242 |
| Felsenstein pruning | the dense `A ⊗ Σ_p` construction on the full tree |
| LOO predictive (for LOO-PIT) | the ratio-of-marginals `log_lik` |
| Stan `log_lik`, `loo_pit_draw`, `lp_between` | the R reference, at real draws |

All agree to machine precision.  The LOO predictive of the fifth check
comes from conditioning on the species effect given the other
observations of the species, while `log_lik` is a ratio of species
marginals.  These are two independent derivations of the same density,
and they agree to 5e-15, so an error in either derivation would show
in this check.

`validate.R` does not yet cover the LOGO code (`log_lik_g`, `holdout`).
