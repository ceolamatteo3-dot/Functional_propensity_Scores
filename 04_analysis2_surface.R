# ============================================================
# STUDY 2:
# FUNCTIONAL BMI EXPOSURE -> LATER FUNCTIONAL OUTCOME
# ============================================================

analysis2_directory <- file.path(
  CFG$out_root,
  "analysis2"
)

set.seed(CFG$seed + 1L)

confounder_names <- confounder_formula(
  CFG$conf_set
)


# ------------------------------------------------------------
# Treatment FPCA and FPS weights
# ------------------------------------------------------------

fpca2 <- dense_fpca(
  curves = study2$X,
  grid = study2$grid_s,
  pve_keep = 0.999
)

L_balance2 <- select_ncomp(
  fpca2,
  CFG$pve_balance
)

A2 <- fpca2$scores[
  ,
  seq_len(L_balance2),
  drop = FALSE
]

confounders2 <- prepare_confounders(
  scalar_covariates = as.data.frame(
    study2$dat[
      ,
      ..confounder_names
    ]
  )
)

fps2 <- fps_weights(
  A = A2,
  C = confounders2$C,
  ridge = CFG$fps_ridge
)

diagnostics2 <- weight_diagnostics(
  fps2
)

utils::write.csv(
  diagnostics2,
  file.path(
    analysis2_directory,
    "weight_diagnostics.csv"
  ),
  row.names = FALSE
)

log_line(
  "S2 ESS = ",
  round(diagnostics2$ESS, 1),
  " (",
  round(
    100 * diagnostics2$ESS_fraction,
    1
  ),
  "% of n); max imbalance = ",
  signif(
    diagnostics2$maximum_moment_imbalance,
    3
  )
)

if (
  diagnostics2$ESS_fraction <
  CFG$min_ess_frac
) {
  log_line(
    "WARNING: S2 ESS fraction is below ",
    CFG$min_ess_frac
  )
}


# ------------------------------------------------------------
# Balance
# ------------------------------------------------------------

balance2 <- fps_balance_table(
  A = A2,
  C = confounders2$C,
  weights = fps2$weights
)

utils::write.csv(
  round(
    balance2$unweighted,
    4
  ),
  file.path(
    analysis2_directory,
    "balance_before.csv"
  )
)

utils::write.csv(
  round(
    balance2$weighted,
    4
  ),
  file.path(
    analysis2_directory,
    "balance_after.csv"
  )
)

grDevices::png(
  file.path(
    analysis2_directory,
    "01_balance.png"
  ),
  width = 1300,
  height = 1000,
  res = 160
)

plot_balance(balance2)

grDevices::dev.off()


# ------------------------------------------------------------
# Weighted and unweighted function-on-function models
# ------------------------------------------------------------

uniform_weights2 <- rep(
  1 / nrow(study2$dat),
  nrow(study2$dat)
)

fit2_weighted <- fit_function_on_function_msm(
  outcome_curves = study2$Y,
  outcome_grid = study2$grid_t,
  treatment_fpca = fpca2,
  weights = fps2$weights,
  pve_treatment_model =
    CFG$pve_balance,
  pve_outcome_model =
    CFG$pve_outcome
)

fit2_unweighted <- fit_function_on_function_msm(
  outcome_curves = study2$Y,
  outcome_grid = study2$grid_t,
  treatment_fpca = fpca2,
  weights = uniform_weights2,
  pve_treatment_model =
    CFG$pve_balance,
  pve_outcome_model =
    CFG$pve_outcome
)

surface_range <- range(
  fit2_weighted$effect_surface,
  fit2_unweighted$effect_surface
)


# ------------------------------------------------------------
# Surface plots
# ------------------------------------------------------------

grDevices::png(
  file.path(
    analysis2_directory,
    "02_surface_weighted.png"
  ),
  width = 1400,
  height = 1100,
  res = 160
)

graphics::filled.contour(
  x = study2$grid_s,
  y = study2$grid_t,
  z = fit2_weighted$effect_surface,
  zlim = surface_range,
  color.palette = function(n) {
    grDevices::hcl.colors(
      n,
      "RdBu",
      rev = TRUE
    )
  },
  xlab = "BMI exposure age s",
  ylab = "Outcome age t",
  main = paste0(
    "FPS-weighted BMI effect surface: ",
    study2$outcome_name
  )
)

grDevices::dev.off()

grDevices::png(
  file.path(
    analysis2_directory,
    "03_surface_unweighted.png"
  ),
  width = 1400,
  height = 1100,
  res = 160
)

graphics::filled.contour(
  x = study2$grid_s,
  y = study2$grid_t,
  z = fit2_unweighted$effect_surface,
  zlim = surface_range,
  color.palette = function(n) {
    grDevices::hcl.colors(
      n,
      "RdBu",
      rev = TRUE
    )
  },
  xlab = "BMI exposure age s",
  ylab = "Outcome age t",
  main = paste0(
    "Unweighted BMI effect surface: ",
    study2$outcome_name
  )
)

grDevices::dev.off()


# ------------------------------------------------------------
# Household-cluster full-pipeline bootstrap
# ------------------------------------------------------------

bootstrap_study2 <- function(B) {
  clusters <- unique(
    study2$dat$hhid
  )
  
  bootstrap_surfaces <- array(
    NA_real_,
    dim = c(
      B,
      length(study2$grid_s),
      length(study2$grid_t)
    )
  )
  
  for (b in seq_len(B)) {
    sampled_clusters <- sample(
      clusters,
      size = length(clusters),
      replace = TRUE
    )
    
    index <- unlist(
      lapply(
        sampled_clusters,
        function(cluster) {
          which(
            study2$dat$hhid == cluster
          )
        }
      ),
      use.names = FALSE
    )
    
    result_b <- tryCatch(
      {
        X_b <- study2$X[
          index,
          ,
          drop = FALSE
        ]
        
        X_b <- sweep(
          X_b,
          2,
          colMeans(X_b),
          "-"
        )
        
        fpca_b <- dense_fpca(
          curves = X_b,
          grid = study2$grid_s,
          pve_keep = 0.999
        )
        
        L_b <- select_ncomp(
          fpca_b,
          CFG$pve_balance
        )
        
        confounders_b <- prepare_confounders(
          scalar_covariates = as.data.frame(
            study2$dat[
              index,
              ..confounder_names
            ]
          )
        )
        
        fps_b <- fps_weights(
          A = fpca_b$scores[
            ,
            seq_len(L_b),
            drop = FALSE
          ],
          C = confounders_b$C,
          ridge = CFG$fps_ridge,
          warn = FALSE
        )
        
        fit_b <- fit_function_on_function_msm(
          outcome_curves =
            study2$Y[
              index,
              ,
              drop = FALSE
            ],
          outcome_grid =
            study2$grid_t,
          treatment_fpca =
            fpca_b,
          weights =
            fps_b$weights,
          pve_treatment_model =
            CFG$pve_balance,
          pve_outcome_model =
            CFG$pve_outcome
        )
        
        fit_b$effect_surface
      },
      error = function(error) {
        NULL
      }
    )
    
    if (!is.null(result_b)) {
      bootstrap_surfaces[b, , ] <-
        result_b
    }
    
    if (
      b %% 25 == 0 ||
      b == B
    ) {
      log_line(
        "S2 bootstrap ",
        b,
        "/",
        B
      )
    }
  }
  
  bootstrap_surfaces
}


bootstrap2 <- bootstrap_study2(
  CFG$n_boot
)

saveRDS(
  bootstrap2,
  file.path(
    analysis2_directory,
    "bootstrap_surfaces.rds"
  )
)


# ------------------------------------------------------------
# Simultaneous surface band
# ------------------------------------------------------------

surface_band2 <- simultaneous_surface_band(
  estimate_surface =
    fit2_weighted$effect_surface,
  bootstrap_surfaces =
    bootstrap2,
  alpha = 0.05,
  studentized = TRUE
)

grDevices::png(
  file.path(
    analysis2_directory,
    "04_significance_map.png"
  ),
  width = 1400,
  height = 1100,
  res = 160
)

graphics::image(
  x = study2$grid_s,
  y = study2$grid_t,
  z = surface_band2$significant + 0,
  col = c(
    "grey92",
    "firebrick"
  ),
  xlab = "BMI exposure age s",
  ylab = "Outcome age t",
  main = paste(
    "Regions where the simultaneous 95% band",
    "excludes zero"
  )
)

graphics::contour(
  x = study2$grid_s,
  y = study2$grid_t,
  z = fit2_weighted$effect_surface,
  add = TRUE
)

graphics::box()

grDevices::dev.off()


# ------------------------------------------------------------
# Effect-surface slices
# ------------------------------------------------------------

extract_slice_band <- function(
    bootstrap_array,
    estimate,
    fixed_dimension,
    index) {
  draws <- if (
    fixed_dimension == "t"
  ) {
    bootstrap_array[, , index]
  } else {
    bootstrap_array[, index, ]
  }
  
  if (is.null(dim(draws))) {
    draws <- matrix(
      draws,
      nrow = 1
    )
  }
  
  valid <- apply(
    draws,
    1,
    function(row) {
      all(is.finite(row))
    }
  )
  
  draws <- draws[
    valid,
    ,
    drop = FALSE
  ]
  
  if (nrow(draws) < 20L) {
    stop(
      "Fewer than 20 valid bootstrap draws for a slice."
    )
  }
  
  lower_percentile <- apply(
    draws,
    2,
    stats::quantile,
    0.025,
    na.rm = TRUE
  )
  
  upper_percentile <- apply(
    draws,
    2,
    stats::quantile,
    0.975,
    na.rm = TRUE
  )
  
  list(
    lower =
      2 * estimate -
      upper_percentile,
    upper =
      2 * estimate -
      lower_percentile
  )
}


fixed_outcome_ages <- c(
  CFG$s2_out_lo,
  mean(
    c(
      CFG$s2_out_lo,
      CFG$s2_out_hi
    )
  ),
  CFG$s2_out_hi
)

fixed_exposure_ages <- c(
  CFG$s2_expo_lo,
  mean(
    c(
      CFG$s2_expo_lo,
      CFG$s2_expo_hi
    )
  ),
  CFG$s2_expo_hi
)

grDevices::png(
  file.path(
    analysis2_directory,
    "05_surface_slices.png"
  ),
  width = 2000,
  height = 1300,
  res = 160
)

graphics::par(
  mfrow = c(2, 3)
)

for (outcome_age in fixed_outcome_ages) {
  j <- which.min(
    abs(
      study2$grid_t -
        outcome_age
    )
  )
  
  estimate <- fit2_weighted$
    effect_surface[, j]
  
  interval <- extract_slice_band(
    bootstrap_array = bootstrap2,
    estimate = estimate,
    fixed_dimension = "t",
    index = j
  )
  
  graphics::plot(
    study2$grid_s,
    estimate,
    type = "n",
    ylim = range(
      interval$lower,
      interval$upper,
      0,
      finite = TRUE
    ),
    xlab = "BMI exposure age s",
    ylab = expression(
      hat(mu)(s, t[0])
    ),
    main = paste0(
      "Outcome age t = ",
      round(
        study2$grid_t[j],
        1
      )
    )
  )
  
  graphics::polygon(
    c(
      study2$grid_s,
      rev(study2$grid_s)
    ),
    c(
      interval$lower,
      rev(interval$upper)
    ),
    border = NA,
    col = grDevices::rgb(
      0.2,
      0.4,
      0.9,
      0.2
    )
  )
  
  graphics::lines(
    study2$grid_s,
    estimate,
    col = "navy",
    lwd = 3
  )
  
  graphics::lines(
    study2$grid_s,
    fit2_unweighted$
      effect_surface[, j],
    col = "firebrick",
    lty = 2,
    lwd = 2
  )
  
  graphics::abline(
    h = 0,
    lty = 3
  )
}

for (exposure_age in fixed_exposure_ages) {
  i <- which.min(
    abs(
      study2$grid_s -
        exposure_age
    )
  )
  
  estimate <- fit2_weighted$
    effect_surface[i, ]
  
  interval <- extract_slice_band(
    bootstrap_array = bootstrap2,
    estimate = estimate,
    fixed_dimension = "s",
    index = i
  )
  
  graphics::plot(
    study2$grid_t,
    estimate,
    type = "n",
    ylim = range(
      interval$lower,
      interval$upper,
      0,
      finite = TRUE
    ),
    xlab = "Outcome age t",
    ylab = expression(
      hat(mu)(s[0], t)
    ),
    main = paste0(
      "Exposure age s = ",
      round(
        study2$grid_s[i],
        1
      )
    )
  )
  
  graphics::polygon(
    c(
      study2$grid_t,
      rev(study2$grid_t)
    ),
    c(
      interval$lower,
      rev(interval$upper)
    ),
    border = NA,
    col = grDevices::rgb(
      0.2,
      0.4,
      0.9,
      0.2
    )
  )
  
  graphics::lines(
    study2$grid_t,
    estimate,
    col = "navy",
    lwd = 3
  )
  
  graphics::lines(
    study2$grid_t,
    fit2_unweighted$
      effect_surface[i, ],
    col = "firebrick",
    lty = 2,
    lwd = 2
  )
  
  graphics::abline(
    h = 0,
    lty = 3
  )
}

grDevices::dev.off()


# ------------------------------------------------------------
# Outcome intercept function
# ------------------------------------------------------------

grDevices::png(
  file.path(
    analysis2_directory,
    "06_intercept_function.png"
  ),
  width = 1300,
  height = 900,
  res = 160
)

graphics::plot(
  study2$grid_t,
  fit2_weighted$intercept_function,
  type = "l",
  lwd = 3,
  col = "navy",
  xlab = "Outcome age t",
  ylab = expression(
    hat(mu)[0](t)
  ),
  main = paste(
    "Estimated mean",
    study2$outcome_name,
    "trajectory"
  )
)

graphics::lines(
  study2$grid_t,
  colMeans(study2$Y),
  col = "grey50",
  lty = 2,
  lwd = 2
)

graphics::legend(
  "topright",
  c(
    "Model intercept",
    "Empirical mean"
  ),
  col = c(
    "navy",
    "grey50"
  ),
  lty = c(
    1,
    2
  ),
  lwd = c(
    3,
    2
  ),
  bty = "n"
)

grDevices::dev.off()


# ------------------------------------------------------------
# Save numerical surface
# ------------------------------------------------------------

surface_results <- expand.grid(
  s = study2$grid_s,
  t = study2$grid_t
)

surface_results$mu_weighted <-
  as.vector(
    fit2_weighted$effect_surface
  )

surface_results$mu_unweighted <-
  as.vector(
    fit2_unweighted$effect_surface
  )

surface_results$simultaneous_lower <-
  as.vector(
    surface_band2$lower
  )

surface_results$simultaneous_upper <-
  as.vector(
    surface_band2$upper
  )

surface_results$simultaneous_excludes_zero <-
  as.integer(
    as.vector(
      surface_band2$significant
    )
  )

utils::write.csv(
  surface_results,
  file.path(
    analysis2_directory,
    "effect_surface.csv"
  ),
  row.names = FALSE
)

writeLines(
  c(
    paste(
      "n:",
      nrow(study2$dat)
    ),
    paste(
      "functional outcome:",
      study2$outcome_name
    ),
    paste(
      "exposure age window:",
      CFG$s2_expo_lo,
      "-",
      CFG$s2_expo_hi
    ),
    paste(
      "outcome age window:",
      CFG$s2_out_lo,
      "-",
      CFG$s2_out_hi
    ),
    paste(
      "treatment FPCs:",
      fit2_weighted$treatment_npc
    ),
    paste(
      "outcome FPCs:",
      fit2_weighted$outcome_npc
    ),
    paste(
      "ESS:",
      round(
        diagnostics2$ESS,
        1
      )
    ),
    paste(
      "ESS fraction:",
      round(
        diagnostics2$ESS_fraction,
        4
      )
    ),
    paste(
      "maximum absolute correlation before:",
      round(
        balance2$maximum_before,
        4
      )
    ),
    paste(
      "maximum absolute correlation after:",
      round(
        balance2$maximum_after,
        4
      )
    ),
    paste(
      "simultaneous surface critical value:",
      round(
        surface_band2$critical_value,
        3
      )
    ),
    paste(
      "proportion of surface excluding zero:",
      round(
        mean(surface_band2$significant),
        4
      )
    ),
    "",
    paste(
      "The exposure and outcome age domains are disjoint.",
      "All exposure ages precede all outcome ages."
    ),
    paste(
      "Therefore, no additional historical or",
      "non-anticipativity restriction is needed."
    ),
    paste(
      "The analysis is conditional on surviving and being",
      "observed through the outcome window."
    )
  ),
  file.path(
    analysis2_directory,
    "summary.txt"
  )
)

log_line(
  "Study 2 analysis complete."
)
