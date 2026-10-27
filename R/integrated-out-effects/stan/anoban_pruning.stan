// ── M1 / M3: both species-level effects integrated out ──
//
//   cor_t = 0, cor_p = 0  ->  M1   ~ state + (1   | gr(binominal, cov)) + (1   | taxa)
//   cor_t = 1, cor_p = 1  ->  M3   ~ state + (1|p| gr(binominal, cov)) + (1|q| taxa)
//
// Model (species g, observation i, D = 3 responses):
//
//   y_gi = mu_g + u_g + v_g + e_gi
//   e_gi ~ N(0, Sigma_r)          residual
//   v_g  ~ N(0, Sigma_t)          iid species effect
//   u    ~ N(0, A (x) Sigma_p)    phylogenetic effect
//
// v_g and the residuals are integrated out by the compound-symmetry
// decomposition (see anoban_functions.stanfunctions).  The
// phylogenetic effect is integrated out by Felsenstein pruning.
//
// Why pruning: A = vcv.phylo(tree, corr = TRUE) is the covariance of
// Brownian motion started at 0 at the root with branch lengths divided
// by the (constant, because the tree is ultrametric) root-to-tip
// height.  Marginalising u directly would need a dense 3903 x 3903
// Cholesky per gradient evaluation (~2e10 flops); pruning does it
// exactly in O(#nodes * D^3) — 2377 nodes of 3x3 algebra.  The dense
// matrix A is never formed: only the tree arrays are needed.
//
// Messages are carried in canonical (information) form:
//
//   leaf g           (m, V) = (ybar_g - mu_g, Sigma_r/n_g + Sigma_t)
//   up an edge of t  V <- V + t Sigma_p                (convolution)
//   into a parent    P += V^-1, h += V^-1 m, g += const  (product)
//   at a node        V = P^-1, m = V h, + product constant
//   at the root      u_root = 0, so the canonical form evaluates to 1
//                    and the log density is the accumulated constant
//
// Generated quantities run the downward (root-to-tips) pass, which
// yields p(u_g | y_-g, theta) = N(m_g, V_g) for every tip in a single
// sweep — because "all of the parent's other children" is a
// subtraction of the stored upward contributions, not a
// recomputation.  The LOO factors then follow with Sigma_b[g] =
// Sigma_t + V_g and the species mean shifted by m_g.
//
// This model has ~42 parameters where the equivalent brms model has
// 7847.  Sampling the 3903 phylogenetic effects instead (non-centred)
// was tried and fails: it hit the maximum treedepth on 50 of 50
// transitions and returned ess_bulk = 2.
//
// Priors are the brms 2.23.1 defaults for this data.

functions {
#include anoban_functions.stanfunctions
}

data {
  int<lower=1> N;
  int<lower=1> G;
  int<lower=1> D;
  int<lower=0> Kc;

  matrix[G, Kc] Xc;
  vector[Kc] means_X;
  matrix[N, D] Y;
  matrix[G, D] ybar;
  array[N] int<lower=1, upper=G> gidx;
  array[G] int<lower=1> ng;

  int<lower=1> U;
  array[U] int<lower=1> n_uniq;
  array[G] int<lower=1, upper=U> uidx;

  int<lower=1> max_ng;
  array[G, max_ng] int obs_idx;

  matrix[D, D] W;
  int<lower=0> N_minus_G;
  real sum_half_D_log_ng;

  // ── tree ──
  int<lower=1> n_node;                   // tips + internal nodes
  int<lower=1> root;
  array[n_node] int<lower=0> parent;     // 0 for the root
  vector<lower=0>[n_node] edge_len;      // scaled by tree height; 0 for root
  array[n_node] int<lower=0, upper=1> is_leaf;
  array[n_node] int<lower=0> species_of; // 0 for internal nodes
  array[n_node] int<lower=1> post_order; // children before parents
  array[n_node] int<lower=1> pre_order;  // parents before children

  vector[D] prior_loc;
  vector[D] prior_scale;

  int<lower=0, upper=1> cor_t;
  int<lower=0, upper=1> cor_p;
  int<lower=0, upper=1> prior_only;
  int<lower=0, upper=1> compute_loo;
  // Per-draw LOO predictive CDF values for analytic LOO-PIT; see
  // anoban_marginal.stan and R/loo_pit.R.
  int<lower=0, upper=1> compute_pit;
}

transformed data {
  matrix[D, N] Yt = Y';
  array[G] vector[D] ybar_arr;
  matrix[D, D] Lw_data = cholesky_decompose(W);
  real within_const = -0.5 * D * N_minus_G * log(2 * pi()) - sum_half_D_log_ng;
  real ld2pi = D * log(2 * pi());
  for (g in 1 : G) ybar_arr[g] = ybar[g]';
}

parameters {
  matrix[Kc, D] b;
  vector[D] Intercept;
  vector<lower=0>[D] sigma;
  cholesky_factor_corr[D] Lrescor;

  vector<lower=0>[D] sd_t;                      // iid species effect
  array[cor_t] cholesky_factor_corr[D] L_t;

  vector<lower=0>[D] sd_p;                      // phylogenetic effect
  array[cor_p] cholesky_factor_corr[D] L_p;
}

transformed parameters {
  matrix[D, D] LSigma_r = diag_pre_multiply(sigma, Lrescor);
  matrix[D, D] LSigma_t;
  matrix[D, D] LSigma_p;
  real lprior = 0;

  if (cor_t == 1) {
    LSigma_t = diag_pre_multiply(sd_t, L_t[1]);
  } else {
    LSigma_t = diag_matrix(sd_t);
  }
  if (cor_p == 1) {
    LSigma_p = diag_pre_multiply(sd_p, L_p[1]);
  } else {
    LSigma_p = diag_matrix(sd_p);
  }

  for (d in 1 : D) {
    lprior += student_t_lpdf(Intercept[d] | 3, prior_loc[d], prior_scale[d]);
    lprior += student_t_lpdf(sigma[d] | 3, 0, prior_scale[d])
              - student_t_lccdf(0 | 3, 0, prior_scale[d]);
    lprior += student_t_lpdf(sd_t[d] | 3, 0, prior_scale[d])
              - student_t_lccdf(0 | 3, 0, prior_scale[d]);
    lprior += student_t_lpdf(sd_p[d] | 3, 0, prior_scale[d])
              - student_t_lccdf(0 | 3, 0, prior_scale[d]);
  }
  lprior += lkj_corr_cholesky_lpdf(Lrescor | 1);
  if (cor_t == 1) lprior += lkj_corr_cholesky_lpdf(L_t[1] | 1);
  if (cor_p == 1) lprior += lkj_corr_cholesky_lpdf(L_p[1] | 1);
}

model {
  if (!prior_only) {
    matrix[D, D] Sigma_r = multiply_lower_tri_self_transpose(LSigma_r);
    matrix[D, D] Sigma_t = multiply_lower_tri_self_transpose(LSigma_t);
    matrix[D, D] Sigma_p = multiply_lower_tri_self_transpose(LSigma_p);
    matrix[G, D] MU = rep_matrix(Intercept', G) + Xc * b;

    array[n_node] matrix[D, D] Pacc = rep_array(rep_matrix(0.0, D, D), n_node);
    array[n_node] vector[D] hacc = rep_array(rep_vector(0.0, D), n_node);
    vector[n_node] gacc = rep_vector(0.0, n_node);

    // ── upward pass ──
    for (ii in 1 : n_node) {
      int v = post_order[ii];
      if (v == root) continue;
      matrix[D, D] V;
      vector[D] m;
      real K;
      if (is_leaf[v] == 1) {
        int g = species_of[v];
        V = Sigma_r / ng[g] + Sigma_t;
        m = ybar_arr[g] - MU[g]';
        K = 0;
      } else {
        matrix[D, D] Lp = cholesky_decompose(Pacc[v]);
        V = chol2inv(Lp);
        m = V * hacc[v];
        K = gacc[v] + 0.5 * dot_product(m, hacc[v])
            + 0.5 * (ld2pi - 2 * sum(log(diagonal(Lp))));
      }
      matrix[D, D] Lup = cholesky_decompose(V + edge_len[v] * Sigma_p);
      matrix[D, D] Pw = chol2inv(Lup);
      vector[D] hw = Pw * m;
      int p = parent[v];
      Pacc[p] += Pw;
      hacc[p] += hw;
      gacc[p] += K - 0.5 * (ld2pi + 2 * sum(log(diagonal(Lup)))
                            + dot_product(m, hw));
    }

    // u_root = 0, so the canonical form contributes nothing further.
    target += gacc[root];

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
  matrix[D, D] Sigma_p_out = multiply_lower_tri_self_transpose(LSigma_p);
  real lp_between;                       // upward-pass log density (validation)
  vector[compute_loo ? N : 0] log_lik;   // exact integrated LOO factors
  matrix[compute_pit ? N : 0, D] loo_pit_draw;

  {
    matrix[D, D] Sigma_r = multiply_lower_tri_self_transpose(LSigma_r);
    matrix[D, D] Sigma_t = multiply_lower_tri_self_transpose(LSigma_t);
    matrix[D, D] Sigma_p = multiply_lower_tri_self_transpose(LSigma_p);
    matrix[G, D] MU = rep_matrix(Intercept', G) + Xc * b;

    array[n_node] matrix[D, D] Pacc = rep_array(rep_matrix(0.0, D, D), n_node);
    array[n_node] vector[D] hacc = rep_array(rep_vector(0.0, D), n_node);
    vector[n_node] gacc = rep_vector(0.0, n_node);
    array[n_node] matrix[D, D] Pc = rep_array(rep_matrix(0.0, D, D), n_node);
    array[n_node] vector[D] hc = rep_array(rep_vector(0.0, D), n_node);

    // ── upward pass (same recursion as the model block; repeated
    //    because Stan functions cannot return several typed
    //    accumulators, and the downward pass needs all of them) ──
    for (ii in 1 : n_node) {
      int v = post_order[ii];
      if (v == root) continue;
      matrix[D, D] V;
      vector[D] m;
      real K;
      if (is_leaf[v] == 1) {
        int g = species_of[v];
        V = Sigma_r / ng[g] + Sigma_t;
        m = ybar_arr[g] - MU[g]';
        K = 0;
      } else {
        matrix[D, D] Lp = cholesky_decompose(Pacc[v]);
        V = chol2inv(Lp);
        m = V * hacc[v];
        K = gacc[v] + 0.5 * dot_product(m, hacc[v])
            + 0.5 * (ld2pi - 2 * sum(log(diagonal(Lp))));
      }
      matrix[D, D] Lup = cholesky_decompose(V + edge_len[v] * Sigma_p);
      Pc[v] = chol2inv(Lup);
      hc[v] = Pc[v] * m;
      int p = parent[v];
      Pacc[p] += Pc[v];
      hacc[p] += hc[v];
      gacc[p] += K - 0.5 * (ld2pi + 2 * sum(log(diagonal(Lup)))
                            + dot_product(m, hc[v]));
    }
    lp_between = gacc[root];

    if (compute_loo || compute_pit) {
      // ── downward pass: p(u_v | everything outside v's subtree) ──
      array[n_node] matrix[D, D] Pd = rep_array(rep_matrix(0.0, D, D), n_node);
      array[n_node] vector[D] hd = rep_array(rep_vector(0.0, D), n_node);
      array[G] matrix[D, D] Sigma_b;
      array[G] vector[D] mu_arr;

      for (ii in 1 : n_node) {
        int v = pre_order[ii];
        if (v == root) continue;
        int p = parent[v];
        matrix[D, D] Vb;
        vector[D] mb;
        if (p == root) {
          // u_root is fixed at 0, so the siblings carry no information.
          Vb = rep_matrix(0.0, D, D);
          mb = rep_vector(0.0, D);
        } else {
          // "all of p's other children" = stored total minus this one.
          matrix[D, D] Lb = cholesky_decompose(Pacc[p] - Pc[v] + Pd[p]);
          Vb = chol2inv(Lb);
          mb = Vb * (hacc[p] - hc[v] + hd[p]);
        }
        matrix[D, D] Vv = Vb + edge_len[v] * Sigma_p;
        Pd[v] = inverse_spd(Vv);
        hd[v] = Pd[v] * mb;
        if (is_leaf[v] == 1) {
          int g = species_of[v];
          Sigma_b[g] = Sigma_t + Vv;       // p(u_g | y_-g) folded in
          mu_arr[g] = MU[g]' + mb;
        }
      }

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
}
