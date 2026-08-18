# ============================================================
# NUMERICAL UTILITIES
# ============================================================

trapz_weights <- function(grid) {
  grid <- as.numeric(grid)
  
  if (length(grid) < 2L) {
    stop("grid must contain at least two points.")
  }
  
  if (any(!is.finite(grid))) {
    stop("grid contains non-finite values.")
  }
  
  if (any(diff(grid) <= 0)) {
    stop("grid must be strictly increasing.")
  }
  
  m <- length(grid)
  w <- numeric(m)
  
  w[1] <- (grid[2] - grid[1]) / 2
  w[m] <- (grid[m] - grid[m - 1]) / 2
  
  if (m > 2L) {
    w[2:(m - 1)] <-
      (grid[3:m] - grid[1:(m - 2)]) / 2
  }
  
  w
}


integrate_grid <- function(values, grid) {
  sum(
    trapz_weights(grid) * values,
    na.rm = TRUE
  )
}


ise_function <- function(estimate, truth, grid) {
  integrate_grid(
    (estimate - truth)^2,
    grid
  )
}


ise_surface <- function(estimate,
                        truth,
                        grid_s,
                        grid_t) {
  ws <- trapz_weights(grid_s)
  wt <- trapz_weights(grid_t)
  
  sum(
    outer(ws, wt) * (estimate - truth)^2,
    na.rm = TRUE
  )
}


log_sum_exp <- function(x) {
  xmax <- max(x)
  xmax + log(sum(exp(x - xmax)))
}


softmax <- function(x) {
  xmax <- max(x)
  z <- exp(x - xmax)
  z / sum(z)
}


standardize_matrix <- function(X,
                               tolerance = 1e-10) {
  X <- as.matrix(X)
  
  if (nrow(X) < 2L) {
    stop("At least two rows are needed.")
  }
  
  center <- colMeans(X)
  scale_value <- apply(X, 2, stats::sd)
  
  keep <-
    is.finite(scale_value) &
    scale_value > tolerance
  
  if (!any(keep)) {
    stop("No nonconstant columns remain.")
  }
  
  if (!all(keep)) {
    warning(
      sum(!keep),
      " constant or nearly constant columns were removed."
    )
  }
  
  X_keep <- X[, keep, drop = FALSE]
  center_keep <- center[keep]
  scale_keep <- scale_value[keep]
  
  X_standardized <- sweep(
    X_keep,
    2,
    center_keep,
    "-"
  )
  
  X_standardized <- sweep(
    X_standardized,
    2,
    scale_keep,
    "/"
  )
  
  list(
    data = X_standardized,
    center = center_keep,
    scale = scale_keep,
    kept_columns = keep
  )
}


# ============================================================
# DENSE FPCA
# ============================================================

dense_fpca <- function(curves,
                       grid,
                       pve_keep = 0.999,
                       max_components = NULL,
                       tolerance = 1e-10) {
  curves <- as.matrix(curves)
  grid <- as.numeric(grid)
  
  n <- nrow(curves)
  m <- ncol(curves)
  
  if (length(grid) != m) {
    stop("length(grid) must equal ncol(curves).")
  }
  
  if (n < 2L) {
    stop("At least two subjects are required.")
  }
  
  if (any(!is.finite(curves))) {
    stop("curves contains missing or non-finite values.")
  }
  
  if (pve_keep <= 0 || pve_keep > 1) {
    stop("pve_keep must be in (0,1].")
  }
  
  integration_weights <- trapz_weights(grid)
  sqrt_weights <- sqrt(integration_weights)
  
  mean_function <- colMeans(curves)
  
  centered_curves <- sweep(
    curves,
    2,
    mean_function,
    "-"
  )
  
  weighted_curves <- sweep(
    centered_curves,
    2,
    sqrt_weights,
    "*"
  )
  
  max_rank <- min(n - 1L, m)
  
  if (!is.null(max_components)) {
    max_rank <- min(
      max_rank,
      as.integer(max_components)
    )
  }
  
  svd_fit <- svd(
    weighted_curves,
    nu = max_rank,
    nv = max_rank
  )
  
  eigenvalues_all <- svd_fit$d^2 / (n - 1)
  total_variance <- sum(eigenvalues_all)
  
  positive <- eigenvalues_all > tolerance
  
  if (!any(positive)) {
    stop("The curves have essentially zero variation.")
  }
  
  eigenvalues_positive <- eigenvalues_all[positive]
  
  V <- svd_fit$v[
    ,
    positive,
    drop = FALSE
  ]
  
  cumulative_pve <-
    cumsum(eigenvalues_positive) /
    total_variance
  
  npc <- which(cumulative_pve >= pve_keep)[1]
  
  if (is.na(npc)) {
    npc <- length(eigenvalues_positive)
  }
  
  eigenvalues <- eigenvalues_positive[
    seq_len(npc)
  ]
  
  V <- V[
    ,
    seq_len(npc),
    drop = FALSE
  ]
  
  eigenfunctions <- sweep(
    V,
    1,
    sqrt_weights,
    "/"
  )
  
  scores <- weighted_curves %*% V
  
  colnames(scores) <- paste0(
    "FPC",
    seq_len(npc)
  )
  
  colnames(eigenfunctions) <- paste0(
    "phi",
    seq_len(npc)
  )
  
  result <- list(
    mean = mean_function,
    eigenfunctions = eigenfunctions,
    scores = scores,
    eigenvalues = eigenvalues,
    total_variance = total_variance,
    pve = cumsum(eigenvalues) / total_variance,
    grid = grid,
    integration_weights = integration_weights,
    centered_curves = centered_curves
  )
  
  class(result) <- "dense_fpca"
  result
}


select_ncomp <- function(fpca_object,
                         pve = 0.95) {
  if (pve <= 0 || pve > 1) {
    stop("pve must be in (0,1].")
  }
  
  selected <- which(
    fpca_object$pve >= pve
  )[1]
  
  if (is.na(selected)) {
    selected <- ncol(
      fpca_object$scores
    )
  }
  
  selected
}


reconstruct_fpca <- function(fpca_object,
                             npc = NULL) {
  if (is.null(npc)) {
    npc <- ncol(
      fpca_object$scores
    )
  }
  
  reconstruction <-
    fpca_object$scores[
      ,
      seq_len(npc),
      drop = FALSE
    ] %*%
    t(
      fpca_object$eigenfunctions[
        ,
        seq_len(npc),
        drop = FALSE
      ]
    )
  
  sweep(
    reconstruction,
    2,
    fpca_object$mean,
    "+"
  )
}


# ============================================================
# CONFOUNDER PREPARATION
# ============================================================

scalar_covariate_matrix <- function(C) {
  if (is.null(C)) {
    return(NULL)
  }
  
  if (is.data.frame(C)) {
    X <- stats::model.matrix(
      ~ .,
      data = C,
      na.action = stats::na.fail
    )
    
    if ("(Intercept)" %in% colnames(X)) {
      X <- X[
        ,
        colnames(X) != "(Intercept)",
        drop = FALSE
      ]
    }
    
    return(X)
  }
  
  as.matrix(C)
}


prepare_confounders <- function(
    scalar_covariates = NULL,
    functional_covariates = NULL) {
  C_parts <- list()
  functional_fpca <- list()
  
  C_scalar <- scalar_covariate_matrix(
    scalar_covariates
  )
  
  if (!is.null(C_scalar)) {
    C_parts[["scalar"]] <- C_scalar
  }
  
  if (!is.null(functional_covariates)) {
    for (j in seq_along(functional_covariates)) {
      current <- functional_covariates[[j]]
      
      pve_j <- if (is.null(current$pve)) {
        0.95
      } else {
        current$pve
      }
      
      fpca_j <- dense_fpca(
        curves = current$curves,
        grid = current$grid,
        pve_keep = max(pve_j, 0.999)
      )
      
      K_j <- select_ncomp(
        fpca_j,
        pve_j
      )
      
      score_j <- fpca_j$scores[
        ,
        seq_len(K_j),
        drop = FALSE
      ]
      
      colnames(score_j) <- paste0(
        "D",
        j,
        "_FPC",
        seq_len(K_j)
      )
      
      C_parts[[paste0(
        "functional_",
        j
      )]] <- score_j
      
      functional_fpca[[j]] <- fpca_j
    }
  }
  
  if (length(C_parts) == 0L) {
    stop(
      "At least one scalar or functional ",
      "confounder is required."
    )
  }
  
  C_combined <- do.call(
    cbind,
    C_parts
  )
  
  list(
    C = C_combined,
    functional_fpca = functional_fpca
  )
}


# ============================================================
# FUNCTIONAL PROPENSITY SCORE WEIGHTS
# ============================================================

construct_fps_moments <- function(A, C) {
  A <- as.matrix(A)
  C <- as.matrix(C)
  
  if (nrow(A) != nrow(C)) {
    stop("A and C must have the same number of rows.")
  }
  
  if (any(!is.finite(A)) || any(!is.finite(C))) {
    stop("A and C must not contain missing values.")
  }
  
  A_standardized <- standardize_matrix(A)
  C_standardized <- standardize_matrix(C)
  
  Az <- A_standardized$data
  Cz <- C_standardized$data
  
  colnames(Az) <- paste0(
    "A_",
    seq_len(ncol(Az))
  )
  
  colnames(Cz) <- paste0(
    "C_",
    seq_len(ncol(Cz))
  )
  
  interaction_list <- vector(
    "list",
    ncol(Cz)
  )
  
  for (j in seq_len(ncol(Cz))) {
    interaction_list[[j]] <- sweep(
      Az,
      1,
      Cz[, j],
      "*"
    )
    
    colnames(interaction_list[[j]]) <-
      paste0(
        colnames(Az),
        ":",
        colnames(Cz)[j]
      )
  }
  
  AC <- do.call(
    cbind,
    interaction_list
  )
  
  G <- cbind(Az, Cz, AC)
  
  moment_scale <- apply(
    G,
    2,
    stats::sd
  )
  
  keep <-
    is.finite(moment_scale) &
    moment_scale > 1e-10
  
  if (!any(keep)) {
    stop("No balancing moments remain.")
  }
  
  if (!all(keep)) {
    warning(
      sum(!keep),
      " constant balancing moments were removed."
    )
  }
  
  G_reduced <- G[
    ,
    keep,
    drop = FALSE
  ]
  
  moment_scale <- moment_scale[keep]
  
  G_optimization <- sweep(
    G_reduced,
    2,
    moment_scale,
    "/"
  )
  
  list(
    G = G_reduced,
    G_optimization = G_optimization,
    A_standardized = Az,
    C_standardized = Cz,
    A_scaling = A_standardized,
    C_scaling = C_standardized,
    retained_moments = keep,
    moment_scale = moment_scale
  )
}


fps_weights <- function(A,
                        C,
                        ridge = 0,
                        maxit = 5000,
                        reltol = 1e-10,
                        warn = TRUE) {
  moments <- construct_fps_moments(
    A,
    C
  )
  
  Gopt <- moments$G_optimization
  Graw <- moments$G
  
  d <- ncol(Gopt)
  
  objective <- function(theta) {
    eta <- -drop(Gopt %*% theta)
    
    log_sum_exp(eta) +
      0.5 * ridge * sum(theta^2)
  }
  
  gradient <- function(theta) {
    eta <- -drop(Gopt %*% theta)
    weights <- softmax(eta)
    
    -drop(crossprod(Gopt, weights)) +
      ridge * theta
  }
  
  optimization <- stats::optim(
    par = rep(0, d),
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = list(
      maxit = maxit,
      reltol = reltol
    )
  )
  
  theta_hat <- optimization$par
  eta_hat <- -drop(Gopt %*% theta_hat)
  weights <- softmax(eta_hat)
  
  weighted_moments <- drop(
    crossprod(Graw, weights)
  )
  
  maximum_imbalance <- max(
    abs(weighted_moments)
  )
  
  effective_sample_size <-
    1 / sum(weights^2)
  
  if (warn && optimization$convergence != 0) {
    warning(
      "FPS optimization did not converge successfully: ",
      optimization$message
    )
  }
  
  if (
    warn &&
    effective_sample_size < 0.1 * length(weights)
  ) {
    warning(
      "Effective sample size is below 10% of n."
    )
  }
  
  result <- list(
    weights = weights,
    theta = theta_hat,
    convergence = optimization$convergence,
    message = optimization$message,
    objective = optimization$value,
    weighted_moments = weighted_moments,
    maximum_imbalance = maximum_imbalance,
    effective_sample_size = effective_sample_size,
    moments = moments,
    optimization = optimization
  )
  
  class(result) <- "fps_weights"
  result
}


# ============================================================
# BALANCE AND WEIGHT DIAGNOSTICS
# ============================================================

weighted_mean <- function(x, w) {
  sum(w * x) / sum(w)
}


weighted_covariance <- function(x, y, w) {
  w <- w / sum(w)
  
  mx <- sum(w * x)
  my <- sum(w * y)
  
  sum(
    w * (x - mx) * (y - my)
  )
}


weighted_correlation <- function(x, y, w) {
  covariance <- weighted_covariance(
    x,
    y,
    w
  )
  
  vx <- weighted_covariance(
    x,
    x,
    w
  )
  
  vy <- weighted_covariance(
    y,
    y,
    w
  )
  
  if (
    !is.finite(vx) ||
    !is.finite(vy) ||
    vx <= 0 ||
    vy <= 0
  ) {
    return(NA_real_)
  }
  
  covariance / sqrt(vx * vy)
}


fps_balance_table <- function(A, C, weights) {
  A <- as.matrix(A)
  C <- as.matrix(C)
  
  if (nrow(A) != nrow(C)) {
    stop("A and C must have the same number of rows.")
  }
  
  if (length(weights) != nrow(A)) {
    stop("weights must contain one value per subject.")
  }
  
  unweighted <- matrix(
    NA_real_,
    nrow = ncol(A),
    ncol = ncol(C)
  )
  
  weighted <- unweighted
  
  for (k in seq_len(ncol(A))) {
    for (j in seq_len(ncol(C))) {
      unweighted[k, j] <- stats::cor(
        A[, k],
        C[, j]
      )
      
      weighted[k, j] <- weighted_correlation(
        A[, k],
        C[, j],
        weights
      )
    }
  }
  
  rownames(unweighted) <- colnames(A)
  colnames(unweighted) <- colnames(C)
  
  rownames(weighted) <- colnames(A)
  colnames(weighted) <- colnames(C)
  
  list(
    unweighted = unweighted,
    weighted = weighted,
    maximum_before = max(
      abs(unweighted),
      na.rm = TRUE
    ),
    maximum_after = max(
      abs(weighted),
      na.rm = TRUE
    ),
    proportion_below_0.1_before = mean(
      abs(unweighted) < 0.1,
      na.rm = TRUE
    ),
    proportion_below_0.1_after = mean(
      abs(weighted) < 0.1,
      na.rm = TRUE
    )
  )
}


plot_balance <- function(balance_object) {
  before <- as.vector(
    abs(balance_object$unweighted)
  )
  
  after <- as.vector(
    abs(balance_object$weighted)
  )
  
  graphics::plot(
    before,
    after,
    pch = 19,
    col = grDevices::rgb(
      0.1,
      0.3,
      0.8,
      0.55
    ),
    xlab = "Absolute correlation before weighting",
    ylab = "Absolute correlation after weighting",
    main = "Treatment-score/confounder balance"
  )
  
  graphics::abline(
    0,
    1,
    lty = 2
  )
  
  graphics::abline(
    h = 0.1,
    col = "red",
    lty = 3
  )
}


weight_diagnostics <- function(fps_object) {
  w <- fps_object$weights
  n <- length(w)
  normalized_w <- n * w
  
  data.frame(
    n = n,
    ESS = 1 / sum(w^2),
    ESS_fraction =
      (1 / sum(w^2)) / n,
    minimum = min(normalized_w),
    q01 = unname(
      stats::quantile(
        normalized_w,
        0.01
      )
    ),
    q05 = unname(
      stats::quantile(
        normalized_w,
        0.05
      )
    ),
    median = stats::median(normalized_w),
    mean = mean(normalized_w),
    q95 = unname(
      stats::quantile(
        normalized_w,
        0.95
      )
    ),
    q99 = unname(
      stats::quantile(
        normalized_w,
        0.99
      )
    ),
    maximum = max(normalized_w),
    maximum_moment_imbalance =
      fps_object$maximum_imbalance,
    optimizer_convergence =
      fps_object$convergence
  )
}


# ============================================================
# MARGINAL STRUCTURAL OUTCOME MODELS
# ============================================================

fit_scalar_functional_msm <- function(
    Y,
    treatment_fpca,
    weights,
    pve_outcome_model = 0.95) {
  Lstar <- select_ncomp(
    treatment_fpca,
    pve_outcome_model
  )
  
  A <- treatment_fpca$scores[
    ,
    seq_len(Lstar),
    drop = FALSE
  ]
  
  design <- cbind(
    "(Intercept)" = 1,
    A
  )
  
  fit <- stats::lm.wfit(
    x = design,
    y = Y,
    w = length(weights) * weights
  )
  
  coefficients <- fit$coefficients
  beta_scores <- coefficients[-1]
  
  phi <- treatment_fpca$eigenfunctions[
    ,
    seq_len(Lstar),
    drop = FALSE
  ]
  
  effect_function <- drop(
    phi %*% beta_scores
  )
  
  list(
    coefficients = coefficients,
    score_coefficients = beta_scores,
    effect_function = effect_function,
    intercept = coefficients[1],
    grid = treatment_fpca$grid,
    Lstar = Lstar,
    fit = fit
  )
}


fit_binary_functional_msm <- function(
    Y,
    treatment_fpca,
    weights,
    pve_outcome_model = 0.95) {
  if (!all(Y %in% c(0, 1))) {
    stop("Binary outcome Y must contain only 0 and 1.")
  }
  
  Lstar <- select_ncomp(
    treatment_fpca,
    pve_outcome_model
  )
  
  A <- treatment_fpca$scores[
    ,
    seq_len(Lstar),
    drop = FALSE
  ]
  
  design <- cbind(
    "(Intercept)" = 1,
    A
  )
  
  fit <- suppressWarnings(
    stats::glm.fit(
      x = design,
      y = Y,
      weights = length(weights) * weights,
      family = stats::binomial()
    )
  )
  
  coefficients <- fit$coefficients
  
  if (any(!is.finite(coefficients))) {
    stop(
      "Binary MSM produced non-finite coefficients, ",
      "possibly because of separation."
    )
  }
  
  beta_scores <- coefficients[-1]
  
  phi <- treatment_fpca$eigenfunctions[
    ,
    seq_len(Lstar),
    drop = FALSE
  ]
  
  effect_function <- drop(
    phi %*% beta_scores
  )
  
  list(
    coefficients = coefficients,
    score_coefficients = beta_scores,
    effect_function = effect_function,
    intercept = coefficients[1],
    grid = treatment_fpca$grid,
    Lstar = Lstar,
    fit = fit
  )
}


fit_function_on_function_msm <- function(
    outcome_curves,
    outcome_grid,
    treatment_fpca,
    weights,
    pve_treatment_model = 0.95,
    pve_outcome_model = 0.95) {
  outcome_curves <- as.matrix(
    outcome_curves
  )
  
  if (
    nrow(outcome_curves) !=
    nrow(treatment_fpca$scores)
  ) {
    stop(
      "Treatment and outcome data must contain ",
      "the same subjects in the same order."
    )
  }
  
  outcome_fpca <- dense_fpca(
    curves = outcome_curves,
    grid = outcome_grid,
    pve_keep = max(
      pve_outcome_model,
      0.999
    )
  )
  
  Ls <- select_ncomp(
    treatment_fpca,
    pve_treatment_model
  )
  
  Lt <- select_ncomp(
    outcome_fpca,
    pve_outcome_model
  )
  
  A <- treatment_fpca$scores[
    ,
    seq_len(Ls),
    drop = FALSE
  ]
  
  CY <- outcome_fpca$scores[
    ,
    seq_len(Lt),
    drop = FALSE
  ]
  
  design <- cbind(
    "(Intercept)" = 1,
    A
  )
  
  weighted_design <- sweep(
    design,
    1,
    weights,
    "*"
  )
  
  weighted_CY <- sweep(
    CY,
    1,
    weights,
    "*"
  )
  
  XtWX <- crossprod(
    design,
    weighted_design
  )
  
  XtWY <- crossprod(
    design,
    weighted_CY
  )
  
  coefficient_matrix <- qr.solve(
    XtWX,
    XtWY
  )
  
  intercept_scores <-
    coefficient_matrix[1, ]
  
  beta_matrix <-
    coefficient_matrix[
      -1,
      ,
      drop = FALSE
    ]
  
  phi_x <- treatment_fpca$eigenfunctions[
    ,
    seq_len(Ls),
    drop = FALSE
  ]
  
  psi_y <- outcome_fpca$eigenfunctions[
    ,
    seq_len(Lt),
    drop = FALSE
  ]
  
  effect_surface <-
    phi_x %*%
    beta_matrix %*%
    t(psi_y)
  
  intercept_function <-
    outcome_fpca$mean +
    drop(
      psi_y %*% intercept_scores
    )
  
  list(
    coefficient_matrix = coefficient_matrix,
    beta_matrix = beta_matrix,
    effect_surface = effect_surface,
    intercept_function = intercept_function,
    treatment_grid = treatment_fpca$grid,
    outcome_grid = outcome_grid,
    treatment_npc = Ls,
    outcome_npc = Lt,
    outcome_fpca = outcome_fpca
  )
}


# ============================================================
# SIMULTANEOUS CONFIDENCE BANDS
# ============================================================

simultaneous_bootstrap_band <- function(
    estimate,
    bootstrap_estimates,
    alpha = 0.05,
    studentized = TRUE,
    minimum_sd = NULL,
    center_bootstrap = FALSE) {
  estimate <- as.numeric(estimate)
  bootstrap_estimates <- as.matrix(
    bootstrap_estimates
  )
  
  if (
    ncol(bootstrap_estimates) !=
    length(estimate)
  ) {
    stop(
      "Bootstrap estimates are not evaluated ",
      "on the same grid as estimate."
    )
  }
  
  valid_rows <- apply(
    bootstrap_estimates,
    1,
    function(x) all(is.finite(x))
  )
  
  bootstrap_estimates <-
    bootstrap_estimates[
      valid_rows,
      ,
      drop = FALSE
    ]
  
  B_valid <- nrow(
    bootstrap_estimates
  )
  
  if (B_valid < 20L) {
    stop(
      "Fewer than 20 valid bootstrap replicates."
    )
  }
  
  if (B_valid < 500L) {
    warning(
      "Only ",
      B_valid,
      " valid bootstrap replicates. ",
      "Use at least 500, preferably 1000, ",
      "for final simultaneous inference."
    )
  }
  
  bootstrap_center <- if (center_bootstrap) {
    colMeans(bootstrap_estimates)
  } else {
    estimate
  }
  
  deviations <- sweep(
    bootstrap_estimates,
    2,
    bootstrap_center,
    "-"
  )
  
  pointwise_se <- apply(
    bootstrap_estimates,
    2,
    stats::sd
  )
  
  if (studentized) {
    if (is.null(minimum_sd)) {
      positive_sd <- pointwise_se[
        is.finite(pointwise_se) &
          pointwise_se > 0
      ]
      
      if (length(positive_sd) == 0L) {
        stop(
          "All bootstrap standard errors are zero."
        )
      }
      
      minimum_sd <- max(
        1e-10,
        stats::quantile(
          positive_sd,
          0.01,
          names = FALSE
        ) * 1e-3
      )
    }
    
    pointwise_se <- pmax(
      pointwise_se,
      minimum_sd
    )
    
    standardized_deviations <- sweep(
      abs(deviations),
      2,
      pointwise_se,
      "/"
    )
    
    sup_statistics <- apply(
      standardized_deviations,
      1,
      max
    )
    
    critical_value <- unname(
      stats::quantile(
        sup_statistics,
        probs = 1 - alpha,
        type = 8
      )
    )
    
    lower <-
      estimate -
      critical_value * pointwise_se
    
    upper <-
      estimate +
      critical_value * pointwise_se
  } else {
    sup_statistics <- apply(
      abs(deviations),
      1,
      max
    )
    
    critical_value <- unname(
      stats::quantile(
        sup_statistics,
        probs = 1 - alpha,
        type = 8
      )
    )
    
    lower <- estimate - critical_value
    upper <- estimate + critical_value
  }
  
  list(
    estimate = estimate,
    lower = lower,
    upper = upper,
    pointwise_se = pointwise_se,
    critical_value = critical_value,
    sup_statistics = sup_statistics,
    alpha = alpha,
    simultaneous_level = 1 - alpha,
    valid_bootstrap_replicates = B_valid,
    failed_bootstrap_replicates =
      sum(!valid_rows),
    studentized = studentized
  )
}


simultaneous_surface_band <- function(
    estimate_surface,
    bootstrap_surfaces,
    alpha = 0.05,
    studentized = TRUE) {
  estimate_surface <- as.matrix(
    estimate_surface
  )
  
  if (
    length(dim(bootstrap_surfaces)) != 3L
  ) {
    stop(
      "bootstrap_surfaces must have dimensions ",
      "B x length(grid_s) x length(grid_t)."
    )
  }
  
  B <- dim(bootstrap_surfaces)[1]
  ns <- nrow(estimate_surface)
  nt <- ncol(estimate_surface)
  
  if (
    any(
      dim(bootstrap_surfaces)[2:3] !=
      c(ns, nt)
    )
  ) {
    stop(
      "Bootstrap surfaces have incompatible dimensions."
    )
  }
  
  bootstrap_matrix <- matrix(
    NA_real_,
    nrow = B,
    ncol = ns * nt
  )
  
  for (b in seq_len(B)) {
    bootstrap_matrix[b, ] <-
      as.vector(
        bootstrap_surfaces[b, , ]
      )
  }
  
  band <- simultaneous_bootstrap_band(
    estimate = as.vector(estimate_surface),
    bootstrap_estimates = bootstrap_matrix,
    alpha = alpha,
    studentized = studentized
  )
  
  list(
    estimate = estimate_surface,
    lower = matrix(
      band$lower,
      nrow = ns,
      ncol = nt
    ),
    upper = matrix(
      band$upper,
      nrow = ns,
      ncol = nt
    ),
    standard_error = matrix(
      band$pointwise_se,
      nrow = ns,
      ncol = nt
    ),
    critical_value = band$critical_value,
    significant = matrix(
      band$lower > 0 |
        band$upper < 0,
      nrow = ns,
      ncol = nt
    ),
    details = band
  )
}
