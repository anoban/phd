## ── Build the Stan data from the csv + the tree ─────────────────────
##
## Produces everything both Stan models need.  The marginal model
## (M2/M4) ignores the tree entries; extra entries in the data list are
## harmless.
##
## The species-level sufficient statistics are what make these models
## cheap.  `state` is constant within species and there are no
## observation-level predictors, so the within-species part of the
## likelihood depends on the data only through the pooled scatter
## W = sum_g S_g and N - G, and the between-species part only through
## the species means.  build_anoban_data() asserts this; see the README.

suppressPackageStartupMessages(library(ape))

#' Build the Stan data list.
#'
#' @param csv_path  path to stan.csv
#' @param tree_path path to FRED4_1301.tre
#' @param log_response model log(y) instead of y (default TRUE; see README)
#' @param responses response column names
#' @param group_col species column
#' @param predictor right-hand side, must be species-constant
#' @return list(standata, species, log_response, predictor_names)
build_anoban_data <- function(csv_path, tree_path, log_response = TRUE, responses = c("F00727", "F00679", "F00709"), group_col = "binominal", predictor = "state") {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  D <- length(responses)

  ## ── structural checks ────────────────────────────────────────────
  ## `predictor = "1"` is the intercept-only right-hand side and names
  ## no column.
  pred_cols <- if (identical(trimws(predictor), "1")) character(0) else predictor
  need <- c(responses, group_col, pred_cols)
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("missing columns: ", paste(miss, collapse = ", "))
  if (anyNA(df[, need])) stop("missing values are not supported")
  ## M1/M3 put a phylogenetic and an iid effect on the same grouping
  ## factor.  In this data set `binominal` and `taxa` are identical, so
  ## that is the intended decomposition.  If they ever differ, the
  ## block structure changes and these models do not apply.
  if ("taxa" %in% names(df) &&
      !identical(as.character(df[[group_col]]), as.character(df$taxa))) {
    stop("`", group_col, "` and `taxa` differ; see the README — these ",
         "models assume they are the same grouping factor")
  }

  species <- sort(unique(as.character(df[[group_col]])))
  G <- length(species)
  N <- nrow(df)
  gidx <- match(as.character(df[[group_col]]), species)

  ## ── responses ────────────────────────────────────────────────────
  Y <- as.matrix(df[, responses])
  if (log_response) {
    if (any(Y <= 0)) stop("log_response = TRUE needs strictly positive responses")
    Y <- log(Y)
  }

  ## ── design matrix, centered as brms does ─────────────────────────
  X <- stats::model.matrix(stats::as.formula(paste("~", predictor)), data = df)
  means_X <- colMeans(X[, -1, drop = FALSE])
  Xc <- sweep(X[, -1, drop = FALSE], 2, means_X)
  Kc <- ncol(Xc)

  ## The sufficient-statistic collapse requires a species-constant
  ## design.  Fail loudly rather than return a silently wrong model.
  ## Vacuous when there are no predictors at all (`predictor = "1"`),
  ## where unique() on a zero-column matrix returns zero rows and would
  ## otherwise trip the check.
  const_ok <- Kc == 0L || all(vapply(split(seq_len(N), gidx), function(ii) {
    nrow(unique(Xc[ii, , drop = FALSE])) == 1L
  }, logical(1)))
  if (!const_ok) {
    stop("the design matrix varies within a species. These models assume\n",
         "  species-constant predictors (see the README). An\n",
         "  observation-level covariate needs a different derivation.")
  }
  first_row <- vapply(split(seq_len(N), gidx), `[`, integer(1), 1L)
  Xc_sp <- Xc[first_row[as.character(seq_len(G))], , drop = FALSE]
  dimnames(Xc_sp) <- NULL

  ## ── sufficient statistics ────────────────────────────────────────
  ng <- as.integer(table(factor(gidx, levels = seq_len(G))))
  ybar <- t(vapply(seq_len(G), function(g) colMeans(Y[gidx == g, , drop = FALSE]),
                   numeric(D)))
  W <- crossprod(Y - ybar[gidx, , drop = FALSE])   # sum_g S_g

  n_uniq <- sort(unique(ng))
  uidx <- match(ng, n_uniq)

  max_ng <- max(ng)
  obs_idx <- matrix(0L, nrow = G, ncol = max_ng)
  ptr <- integer(G)
  for (n in seq_len(N)) {
    g <- gidx[n]
    ptr[g] <- ptr[g] + 1L
    obs_idx[g, ptr[g]] <- n
  }

  ## ── brms default priors, recomputed on the modelled scale ────────
  ## student_t(3, median(y), max(2.5, mad(y))), rounded to 1 decimal.
  prior_loc <- round(apply(Y, 2, stats::median), 1)
  prior_scale <- round(pmax(2.5, apply(Y, 2, stats::mad)), 1)

  ## ── tree arrays for the pruning model ────────────────────────────
  ta <- build_tree_arrays(read.tree(tree_path), species)

  standata <- list(
    N = N, G = G, D = D, Kc = Kc,
    Xc = unname(as.matrix(Xc_sp)), means_X = as.numeric(means_X),
    Y = unname(as.matrix(Y)), ybar = unname(as.matrix(ybar)),
    gidx = gidx, ng = ng,
    U = length(n_uniq), n_uniq = n_uniq, uidx = uidx,
    max_ng = max_ng, obs_idx = obs_idx,
    W = W, N_minus_G = N - G,
    sum_half_D_log_ng = 0.5 * D * sum(log(ng)),
    n_node = ta$n_node, root = ta$root, parent = ta$parent,
    edge_len = ta$edge_len, is_leaf = ta$is_leaf,
    species_of = ta$species_of,
    post_order = ta$post_order, pre_order = ta$pre_order,
    prior_loc = as.numeric(prior_loc),
    prior_scale = as.numeric(prior_scale),
    ## Both Stan models declare these; set defaults here so every caller
    ## gets a complete data list.  run_loo_pit.R switches compute_pit on.
    prior_only = 0L, compute_loo = 1L, compute_pit = 0L
  )

  list(standata = standata, species = species,
       log_response = log_response,
       predictor_names = colnames(Xc),
       tree_height = ta$height)
}

#' Node-level traversal arrays for Felsenstein pruning.
#'
#' `u ~ N(0, A (x) Sigma_p)` with A = vcv.phylo(corr = TRUE) is Brownian
#' motion started at 0 at the root, with branch lengths divided by the
#' (constant, because the tree is ultrametric) root-to-tip height.  Only
#' these arrays are needed — the dense A is never formed.
build_tree_arrays <- function(tree, species) {
  stopifnot(is.rooted(tree), setequal(tree$tip.label, species))
  if (!is.ultrametric(tree)) {
    stop("the tree is not ultrametric; the fixed-root-at-zero scaling ",
         "used here assumes it is")
  }

  n_tip <- length(tree$tip.label)
  n_node <- n_tip + tree$Nnode
  root <- n_tip + 1L
  ## root-to-tip distance, constant because the tree is ultrametric.
  ## Not via vcv.phylo(): no need to build a 1301 x 1301 matrix for one
  ## number.
  height <- mean(node.depth.edgelength(tree)[seq_len(n_tip)])

  parent <- integer(n_node)
  edge_len <- numeric(n_node)
  parent[tree$edge[, 2]] <- tree$edge[, 1]
  edge_len[tree$edge[, 2]] <- tree$edge.length / height
  if (parent[root] != 0L) stop("root index is not the parentless node")

  species_of <- integer(n_node)
  species_of[seq_len(n_tip)] <- match(tree$tip.label, species)
  if (anyNA(species_of[seq_len(n_tip)])) stop("tip labels not in species")

  post_edges <- reorder(tree, "postorder")$edge
  post_order <- c(post_edges[, 2], root)
  if (length(post_order) != n_node || anyDuplicated(post_order)) {
    stop("post-order traversal is not a permutation of the nodes")
  }

  list(n_tip = n_tip, n_node = n_node, root = root, height = height,
       parent = parent, edge_len = edge_len,
       is_leaf = as.integer(seq_len(n_node) <= n_tip),
       species_of = species_of,
       post_order = post_order, pre_order = rev(post_order))
}
