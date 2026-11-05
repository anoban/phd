# anoban — trivariate species models with the species effects integrated out

Stan + R code for the four models, rewritten so that the species-level
effects are integrated out analytically rather than sampled. The
resulting code is mix from brms generated code, code from a related
LOO computational methods research project, new code by Claude and
edit by Aki. This README was Claude generated. Aki has checked and
edited this README and all code comments.

Two things this buys:

- **Speed.** brms samples the species-level effects: 3903 of them in
  M2/M4 (3 responses × 1301 species) and 7806 in M1/M3. Here they are
  integrated out, leaving 27–36 free parameters. Measured by
  `run_all.R` on one cluster node, 4 chains × 1000 warmup + 1000
  sampling, log responses (wall times vary by a factor of ~2 across
  machines, the ESS values do not):

  | | M1 | M2 | M3 | M4 |
  |---|---:|---:|---:|---:|
  | wall (s) | 659 | 32 | 1016 | 60 |
  | min ess_bulk | 2714 | 3207 | 3656 | 4072 |

  with no divergences and no treedepth saturation in any of them. For
  reference, the brms M2 measured elsewhere took 691 s to sample and a
  further 1091 s for `loo()` (raw response scale, so not a like-for-like
  timing). M1/M3 are the ones that change character: sampling 7806
  phylogenetic effects instead hit the maximum treedepth on every
  transition and returned `ess_bulk = 2`.
- **Trustworthy LOO.** In the brms parameterisation, 732 of the 1301
  species are observed exactly once, so leaving out an observation
  destabilises its own species effect. Measured on M2: **787 of 5218
  folds with Pareto k > 0.7**, and `p_loo = 2942`. With the species
  effects integrated out the same model gives **0–2 folds with
  k > 0.7** and `p_loo` of 53–62.

The models themselves are unchanged: same data model, same priors, same
posterior for the parameters you care about. Only the computation
differs.

---

## Quick start

Requirements: `cmdstanr` (with CmdStan installed), `loo`, `posterior`,
`ape`; plus `bayesplot`, `ggplot2` and `patchwork` for the LOO-PIT
plots. No brms needed. Instead of `RScript` the code can be run
also interactively.

```r
Rscript run_all.R                 # fit M1-M4 on log responses, compute LOO
Rscript run_all.R --raw           # original (untransformed) response scale
Rscript run_all.R --models M3     # one model only
Rscript run_all.R --quick         # short chains, for a smoke test
```

Results land in `output/`: one `_fit.rds` and `_loo.rds` per model,
plus `loo_compare.txt`. Per-model elpd_loo, `p_loo` and the Pareto-k
table are printed by `print(loo)` as each model finishes.

LOO-PIT calibration plots, after `run_all.R`:

```r
Rscript run_loo_pit.R             # all fitted models
Rscript run_loo_pit.R --models M3
```

Check the code does what it claims before trusting it:

```r
Rscript R/validate.R              # all checks, compiles Stan (~5 min)
Rscript R/validate.R --no-stan    # R-level checks only, seconds
```

---

## Response scale

**`run_all.R` uses log responses by default.** Pass `--raw` for the
original scale.

`F00727` has skewness 4.3 and a maximum of 840 against a
median of 20. Fitting a Gaussian to that on the raw scale is badly
misspecified, and it shows up in three independent ways:

| | raw scale | log responses |
|---|---:|---:|
| elpd_loo (M3, both on the original `y` scale) | −20150 | **−11759** |
| `p_loo` (M3, against 36 free parameters) | 114 | 62 |
| folds with Pareto k > 0.7 (M3, of 5218) | 4 | **2** |
| worst LOO-PIT departure from uniform, F00727 | 0.169 | **0.053** |

(the fold counts are small integers and move by one or two between
runs; the other rows are stable)

The log models are ~8400 nats better, which is 18 paired standard
errors. The comparison is like-for-like: `loo_compare_original_scale()`
applies the change-of-variables pointwise,

```
elpd_i^Y = elpd_i^Z − Σ_d log y_id
```

so both scales are scored as densities for the same random variable `y`
and `loo_compare()` still reports a paired SE.

Everything else in this package works identically on either scale.

---

## The four models

| | this package | your brms formula |
|---|---|---|
| M1 | `anoban_pruning.stan`, `cor_t=0, cor_p=0` | `~ state + (1 \| gr(binominal, cov = corrmat)) + (1 \| taxa)` |
| M2 | `anoban_marginal.stan`, `cor_t=0` | `~ state + (1 \| taxa)` |
| M3 | `anoban_pruning.stan`, `cor_t=1, cor_p=1` | `~ state + (1\|p\| gr(binominal, cov = corrmat)) + (1\|q\| taxa)` |
| M4 | `anoban_marginal.stan`, `cor_t=1` | `~ state + (1\|q\| taxa)` |

All with `mvbind(F00727, F00679, F00709)` and `set_rescor(TRUE)`.
`cor_t` / `cor_p` = 1 means the effect is correlated **across
responses**, i.e. brms `(1|q|...)` rather than `(1|...)`.

Priors are brms 2.23.1's defaults, recomputed on whichever response
scale you fit: `student_t(3, median(y), max(2.5, mad(y)))` for the
intercepts, half-`student_t(3, 0, max(2.5, mad(y)))` for `sigma` and the
group SDs, `lkj(1)` for all correlation matrices, flat for `b`. On the
raw scale this reproduces brms exactly (locations `20.5, 0.5, 0.3`,
scales `17.1, 2.5, 2.5`).

`corrmat` is never built. The pruning model works from the tree
directly, so there is no 1301×1301 matrix anywhere.

---

## How it works

**Both species effects are integrated out in closed form.** Nothing is
approximated.

*The iid species effect and the residuals* (all four models). `state` is
constant within species and there are no observation-level predictors,
so for species `g` with `n` observations the likelihood splits exactly:

```
log p(y_g) = N_D(ȳ_g | μ_g, Σ_r/n + Σ_t)                    between
           − ½[(n−1) log|Σ_r| + tr(Σ_r⁻¹ S_g) + D(n−1) log 2π]   within
           − (D/2) log n
```

Summed over species, the within part depends on the data only through
the pooled scatter `W = Σ_g S_g` and `N − G`, both computed once in R.
The fitting problem shrinks from 5218 observations to 1301
three-dimensional units.

*The phylogenetic effect* (M1/M3) by **Felsenstein pruning**.
`A = vcv.phylo(tree, corr = TRUE)` is the covariance of Brownian motion
started at 0 at the root with branch lengths divided by the tree height,
so the dense `A` is the marginal covariance of a Markov process on the
tree. Pruning exploits that: messages travel up the tree, gaining
`t·Σ_p` per edge and multiplying at each node, and the whole
marginalisation costs O(#nodes · D³) — 2377 nodes of 3×3 algebra —
instead of a dense 3903×3903 Cholesky per gradient (~2×10¹⁰ flops).
Polytomies need no special handling.

*The LOO factors.* `log_lik[i] = p(y_i | y_{g∖i}, θ)` — the density of
observation i given the **other** observations of its species, species
effects integrated out (this part has not been used elsewhere, but
will be part of a paper to be arxived in fall). A downward
(root-to-tips) pruning pass supplies `p(u_g | y_{−g}, θ)` for all 1301
species in one sweep, so this costs O(1) per observation.

---

## LOO-PIT — is the model calibrated?

`elpd_loo` ranks models against each other; it does not say whether any
of them is any good. LOO-PIT does. Each observation is mapped through
its own leave-one-out predictive CDF; under a calibrated model the
values are Uniform(0, 1), and the shape of the departure says how the
model is wrong.

`run_loo_pit.R` produces one ECDF panel per response using
`bayesplot::ppc_loo_pit_ecdf(method = "correlated", test = "POT",
plot_diff = TRUE)`. The `"correlated"` method accounts for the LOO-PIT
values not being independent (see, Tesso and Vehtari, 2026).
All panels share a y-axis, so the magnitudes are comparable across
responses.

**The PIT values are computed analytically.** Because the LOO predictive
here is a finite mixture of Gaussians over the draws,

```
PIT_id = Σ_s w_s · Φ((y_id − m_sid) / σ_sid)
```

with `w_s` the PSIS weights. Stan emits the per-draw `Φ(·)` term and R
supplies the weights, so no `y_rep` is simulated and the diagnostic
carries no posterior-predictive Monte Carlo noise.

Reading the plots: the line is the ECDF **minus** its uniform
expectation, so calibrated means flat near zero and inside the
band. An S-shape (rising then falling) means PIT values pile up near
0.5 — the predictive is too **wide**.  On the log scale all three
still show a mild S-shape, so the model is much better calibrated, not
well calibrated. Even with the log-scale there seem to be some
outliers (shown as red dots near 0 and 1).

## Assumptions this code makes

Two are asserted in `build_anoban_data()` and will stop with an error;
the third is your modelling decision.

1. **No observation-level predictors.** The sufficient-statistic
   collapse needs the design matrix to be constant within species.
   `state` is. If you add a covariate that varies *within* a species,
   this code no longer applies and the derivation must be redone — you
   would get an error, not a wrong answer.
2. **`binominal` and `taxa` are the same grouping factor.** They are
   identical in `stan.csv`, so M1/M3 decompose one species effect
   into a phylogenetically correlated part and an iid part. If `taxa`
   were ever a coarser level (genus, family) the block structure would
   change.
3. **The tree is ultrametric**, which `FRED4_1301.tre` is (root-to-tip
   height 400.7877 for every tip). This is what makes the fixed-root
   scaling valid.

Also worth knowing: `state` is species-constant and very unbalanced —
1098 of 1301 species are `AM`, and `ErM` has 9, `NM` 16, `EcMAM` 20. The
coefficients for the rare states rest on very little data.

---

## Files

```
run_all.R                          driver: build -> fit -> LOO -> compare
run_loo_pit.R                      LOO-PIT calibration plots
R/build_data.R                     csv + tree -> Stan data (+ the assertions)
R/fit_models.R                     model table, compile and sample
R/loo_compare_scales.R             comparing log- and raw-response fits
R/loo_pit.R                        analytic LOO-PIT + the ECDF plots
R/validate.R                       exactness checks
stan/anoban_functions.stanfunctions  block_marginal, loo_factor_vec, loo_pit_vec
stan/anoban_marginal.stan          M2 / M4
stan/anoban_pruning.stan           M1 / M3
data/stan.csv, data/FRED4_1301.tre
```

`run_loo_pit.R` reads the saved fits and re-runs only the
generated-quantities block (`generate_quantities`, no re-sampling), so
the per-draw PIT terms — `N × D` numbers per draw — never have to be
carried by the main fit.

## What `validate.R` checks

The likelihoods here are derived, so they are checked against
brute-force references rather than trusted:

| check | against |
|---|---|
| compound-symmetry formula | a dense 3n_g multivariate normal, every species |
| pooled-W model-block form | the sum of per-species terms |
| O(1) rank-one LOO factor | O(n) recomputation, group sizes 1 / 2 / 5 / 242 |
| Felsenstein pruning | the dense `A ⊗ Σ_p` construction on the full tree |
| LOO predictive (for LOO-PIT) | the ratio-of-marginals `log_lik` |
| Stan `log_lik`, `loo_pit_draw`, `lp_between` | the R reference, at real draws |

All pass at machine precision. The fifth is worth singling out: the
LOO predictive used for LOO-PIT is derived by conditioning on the
species effect given the retained observations, while `log_lik` is a
ratio of species marginals. They are independent routes to the same
density, and they agree to 5e-15 — so a mistake in either derivation
would show up rather than propagate quietly into the diagnostics.

## Getting LOGO as well

Leave-one-species-out is one line away: the per-species marginal
`log p(y_g | y_{−g}, θ)` is already computed inside `loo_factor_vec`
(it is the `A_n` term plus the within-species part). Emitting it as
`log_lik_g[G]` in generated quantities and passing that to `loo::loo()`
gives PSIS-LOGO. Note that LOGO is harder than LOO here — leaving out a
species with up to 242 observations is a large perturbation, and the
Pareto k values are correspondingly worse (3.3–3.7 in the experiments).
