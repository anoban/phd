// ── M2 / M4: trivariate species model, species effect integrated out ──
//
//   cor_t = 0  ->  M2   mvbind(...) ~ state + (1   | taxa)
//   cor_t = 1  ->  M4   mvbind(...) ~ state + (1|q| taxa)
//
// Model (species g, observation i, D = 3 responses):
//
//   y_gi = mu_g + v_g + e_gi
//   e_gi ~ N(0, Sigma_r),  Sigma_r = diag(sigma) Lrescor Lrescor' diag(sigma)
//   v_g  ~ N(0, Sigma_t)   iid species effect; diagonal if cor_t = 0
//
// `state` is constant within species and there are no
// observation-level predictors, so mu_g is species-constant and the
// likelihood splits exactly (see block_marginal in
// anoban_functions.stanfunctions).  Summed over species, the
// within-species part depends on the data only through the pooled
// scatter W = sum_g S_g and N - G, both computed once in R.
//
// The species effects are therefore NOT sampled: this model has ~33
// parameters where the equivalent brms model has 3933.  Only 38
// distinct group sizes occur, so the between-species Cholesky factors
// are built 38 times per gradient evaluation, not 1301 times.
//
// Priors are the brms 2.23.1 defaults for this data.

functions {
#include anoban_functions.stanfunctions
}

data {
  int<lower=1> N;                        // observations
  int<lower=1> G;                        // species
  int<lower=1> D;                        // responses
  int<lower=0> Kc;                       // centered population-level predictors
                                         // (0 for the intercept-only N models)

  matrix[G, Kc] Xc;                      // species-level centered design
  vector[Kc] means_X;                    // column means used for centering
  matrix[N, D] Y;                        // responses
  matrix[G, D] ybar;                     // species means
  array[N] int<lower=1, upper=G> gidx;   // species of each observation
  array[G] int<lower=1> ng;              // observations per species

  int<lower=1> U;                        // number of distinct group sizes
  array[U] int<lower=1> n_uniq;
  array[G] int<lower=1, upper=U> uidx;

  int<lower=1> max_ng;
  array[G, max_ng] int obs_idx;          // ragged species -> obs map, padded

  matrix[D, D] W;                        // pooled within-species scatter
  int<lower=0> N_minus_G;
  real sum_half_D_log_ng;                // sum_g 0.5 D log n_g

  vector[D] prior_loc;                   // brms default intercept locations
  vector[D] prior_scale;                 // brms default scales (per response)

  int<lower=0, upper=1> cor_t;           // 0 = M2, 1 = M4
  int<lower=0, upper=1> prior_only;
  int<lower=0, upper=1> compute_loo;     // 0 skips log_lik (faster)
  // Per-draw LOO predictive CDF values for analytic LOO-PIT.  N x D
  // per draw, so leave it off during fitting and switch it on in a
  // generate_quantities() pass (see R/loo_pit.R).
  int<lower=0, upper=1> compute_pit;
}

transformed data {
  matrix[D, N] Yt = Y';
  array[G] vector[D] ybar_arr;
  matrix[D, D] Lw_data = cholesky_decompose(W);
  real within_const = -0.5 * D * N_minus_G * log(2 * pi()) - sum_half_D_log_ng;

  // Species grouped by size, so the between-species density can be
  // evaluated with one vectorised call per distinct size.
  array[U] int u_count = rep_array(0, U);
  int max_u_count;
  for (g in 1 : G) u_count[uidx[g]] += 1;
  max_u_count = max(u_count);
  array[U, max_u_count] int u_members = rep_array(0, U, max_u_count);
  {
    array[U] int ptr = rep_array(0, U);
    for (g in 1 : G) {
      ptr[uidx[g]] += 1;
      u_members[uidx[g], ptr[uidx[g]]] = g;
    }
  }
  for (g in 1 : G) ybar_arr[g] = ybar[g]';
}

parameters {
  matrix[Kc, D] b;                                  // population-level effects
  vector[D] Intercept;                              // on centered predictors
  vector<lower=0>[D] sigma;                         // residual scales
  cholesky_factor_corr[D] Lrescor;                  // residual correlations
  vector<lower=0>[D] sd_t;                          // species effect scales
  array[cor_t] cholesky_factor_corr[D] L_t;         // present only if cor_t = 1
}

transformed parameters {
  matrix[D, D] LSigma_r = diag_pre_multiply(sigma, Lrescor);
  matrix[D, D] LSigma_t;
  real lprior = 0;

  if (cor_t == 1) {
    LSigma_t = diag_pre_multiply(sd_t, L_t[1]);
  } else {
    LSigma_t = diag_matrix(sd_t);
  }

  for (d in 1 : D) {
    lprior += student_t_lpdf(Intercept[d] | 3, prior_loc[d], prior_scale[d]);
    lprior += student_t_lpdf(sigma[d] | 3, 0, prior_scale[d])
              - student_t_lccdf(0 | 3, 0, prior_scale[d]);
    lprior += student_t_lpdf(sd_t[d] | 3, 0, prior_scale[d])
              - student_t_lccdf(0 | 3, 0, prior_scale[d]);
  }
  lprior += lkj_corr_cholesky_lpdf(Lrescor | 1);
  if (cor_t == 1) lprior += lkj_corr_cholesky_lpdf(L_t[1] | 1);
}

model {
  if (!prior_only) {
    matrix[D, D] Sigma_r = multiply_lower_tri_self_transpose(LSigma_r);
    matrix[D, D] Sigma_t = multiply_lower_tri_self_transpose(LSigma_t);
    matrix[G, D] MU = rep_matrix(Intercept', G) + Xc * b;
    array[G] vector[D] mu;
    for (g in 1 : G) mu[g] = MU[g]';

    // between-species part, one vectorised call per distinct group size
    for (u in 1 : U) {
      array[u_count[u]] int idx = u_members[u, 1 : u_count[u]];
      matrix[D, D] Lb = cholesky_decompose(Sigma_r / n_uniq[u] + Sigma_t);
      target += multi_normal_cholesky_lpdf(ybar_arr[idx] | mu[idx], Lb);
    }

    // within-species part: depends on the data only through W and N-G
    target += -0.5 * (2 * N_minus_G * sum(log(diagonal(LSigma_r)))
                      + sum(square(mdivide_left_tri_low(LSigma_r, Lw_data))))
              + within_const;
  }

  target += lprior;
}

generated quantities {
  vector[D] b_Intercept = Intercept - (means_X' * b)';
  matrix[D, D] Rescor = multiply_lower_tri_self_transpose(Lrescor);
  matrix[D, D] Sigma_t_out = multiply_lower_tri_self_transpose(LSigma_t);

  // Per-observation LOO factors log p(y_i | y_{g\i}, theta): what
  // loo::loo() needs.  The species effect is integrated out, which is
  // what keeps the Pareto-k values low.
  vector[compute_loo ? N : 0] log_lik;
  matrix[compute_pit ? N : 0, D] loo_pit_draw;

  if (compute_loo || compute_pit) {
    matrix[D, D] Sigma_r = multiply_lower_tri_self_transpose(LSigma_r);
    matrix[D, D] Sigma_t = multiply_lower_tri_self_transpose(LSigma_t);
    matrix[G, D] MU = rep_matrix(Intercept', G) + Xc * b;
    array[G] matrix[D, D] Sigma_b = rep_array(Sigma_t, G);
    array[G] vector[D] mu_arr;
    for (g in 1 : G) mu_arr[g] = MU[g]';
    if (compute_loo) {
      log_lik = loo_factor_vec(N, Yt, obs_idx, ng, ybar_arr, mu_arr,
                               LSigma_r, Sigma_r, Sigma_b);
    }
    if (compute_pit) {
      loo_pit_draw = loo_pit_vec(N, Yt, obs_idx, ng, ybar_arr, mu_arr,
                                 Sigma_r, Sigma_b);
    }
  }
}
