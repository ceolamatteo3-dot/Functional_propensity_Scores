# ============================================================
# STUDY 1:
# FUNCTIONAL BMI EXPOSURE -> INCIDENT DIABETES
# ============================================================

analysis1_directory <- file.path(
  CFG$out_root,
  "analysis1"
)

set.seed(CFG$seed)

confounder_names <- confounder_formula(
  CFG$conf_set
)


# ------------------------------------------------------------
# Treatment FPCA and FPS weights
# ------------------------------------------------------------

fpca1 <- dense_fpca(
  curves = study1$X,
  grid = study1$grid,
  pve_keep = 0.999
)

L_balance1 <- select_ncomp(
  fpca1,
  CFG$pve_balance
)

A1 <- fpca1$scores[
  ,
  seq_len(L_balance1),
  drop = FALSE
]

confounders1 <- prepare_confounders(
  scalar_covariates = as.data.frame(
    study1$dat[
      ,
      ..confounder_names
    ]
  )
)

fps1 <- fps_weights(
  A = A1,
  C = confounders1$C,
  ridge = CFG$fps_ridge
)

diagnostics1 <- weight_diagnostics(
  fps1
)

utils::write.csv(
  diagnostics1,
  file.path(
    analysis1_directory,
    "weight_diagnostics.csv"
  ),
  row.names = FALSE
)

log_line(
  "S1 ESS = ",
  round(diagnostics1$ESS, 1),
  " (",
  round(
    100 * diagnostics1$ESS_fraction,
    1
  ),
  "% of n); max imbalance = ",
  signif(
    diagnostics1$maximum_moment_imbalance,
    3
  )
)

if (
  diagnostics1$ESS_fraction <
  CFG$min_ess_frac
) {
  log_line(
    "WARNING: S1 ESS fraction is below ",
    CFG$min_ess_frac,
    ". Consider fewer confounders, fewer FPCs, ",
    "or positive fps_ridge."
  )
}


# ------------------------------------------------------------
# Balance
# ------------------------------------------------------------

balance1 <- fps_balance_table(
  A = A1,
  C = confounders1$C,
  weights = fps1$weights
)

utils::write.csv(
  round(
    balance1$unweighted,
    4
  ),
  file.path(
    analysis1_directory,
    "balance_before.csv"
  )
)

utils::write.csv(
  round(
    balance1$weighted,
    4
  ),
  file.path(
    analysis1_directory,
    "balance_after.csv"
  )
)

grDevices::png(
  file.path(
    analysis1_directory,
    "01_balance.png"
  ),
  width = 1300,
  height = 1000,
  res = 160
)

plot_balance(balance1)

grDevices::dev.off()


# ------------------------------------------------------------
# Weight distribution
# ------------------------------------------------------------

grDevices::png(
  file.path(
    analysis1_directory,
    "02_weights.png"
  ),
  width = 1400,
  height = 850,
  res = 160
)

graphics::hist(
  length(fps1$weights) * fps1$weights,
  breaks = 40,
  col = "grey80",
  xlab = "n times FPS weight",
  main = "Normalized functional propensity-score weights"
)

graphics::abline(
  v = 1,
  col = "red",
  lty = 2
)

grDevices::dev.off()


# ------------------------------------------------------------
# Weighted and unweighted outcome models
# ------------------------------------------------------------

Y1 <- study1$dat$Y

uniform_weights1 <- rep(
  1 / length(Y1),
  length(Y1)
)

fit1_weighted <- fit_binary_functional_msm(
  Y = Y1,
  treatment_fpca = fpca1,
  weights = fps1$weights,
  pve_outcome_model = CFG$pve_outcome
)

fit1_unweighted <- fit_binary_functional_msm(
  Y = Y1,
  treatment_fpca = fpca1,
  weights = uniform_weights1,
  pve_outcome_model = CFG$pve_outcome
)


# ------------------------------------------------------------
# Household-cluster bootstrap
# ------------------------------------------------------------

bootstrap_study1 <- function(B) {
  clusters <- unique(
    study1$dat$hhid
  )
  
  effect_bootstrap <- matrix(
    NA_real_,
    nrow = B,
    ncol = length(study1$grid)
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
            study1$dat$hhid == cluster
          )
        }
      ),
      use.names = FALSE
    )
    
    effect_bootstrap[b, ] <- tryCatch(
      {
        X_b <- study1$X[
          index,
          ,
          drop = FALSE
        ]
        
        # Re-centre within the bootstrap sample.
        X_b <- sweep(
          X_b,
          2,
          colMeans(X_b),
          "-"
        )
        
        fpca_b <- dense_fpca(
          curves = X_b,
          grid = study1$grid,
          pve_keep = 0.999
        )
        
        L_b <- select_ncomp(
          fpca_b,
          CFG$pve_balance
        )
        
        confounders_b <- prepare_confounders(
          scalar_covariates = as.data.frame(
            study1$dat[
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
        
        fit_b <- fit_binary_functional_msm(
          Y = Y1[index],
          treatment_fpca = fpca_b,
          weights = fps_b$weights,
          pve_outcome_model =
            CFG$pve_outcome
        )
        
        fit_b$effect_function
      },
      error = function(error) {
        rep(
          NA_real_,
          length(study1$grid)
        )
      }
    )
    
    if (
      b %% 25 == 0 ||
      b == B
    ) {
      log_line(
        "S1 bootstrap ",
        b,
        "/",
        B
      )
    }
  }
  
  effect_bootstrap
}


bootstrap1 <- bootstrap_study1(
  CFG$n_boot
)

saveRDS(
  bootstrap1,
  file.path(
    analysis1_directory,
    "bootstrap_draws.rds"
  )
)


# ------------------------------------------------------------
# Pointwise and simultaneous intervals
# ------------------------------------------------------------

valid_counts1 <- colSums(
  is.finite(bootstrap1)
)

if (any(valid_counts1 < 20L)) {
  stop(
    "Fewer than 20 valid Study 1 bootstrap draws ",
    "at one or more grid points."
  )
}

lower_percentile1 <- apply(
  bootstrap1,
  2,
  stats::quantile,
  probs = 0.025,
  na.rm = TRUE
)

upper_percentile1 <- apply(
  bootstrap1,
  2,
  stats::quantile,
  probs = 0.975,
  na.rm = TRUE
)

# Reverse-percentile/basic bootstrap intervals.
lower_basic1 <-
  2 * fit1_weighted$effect_function -
  upper_percentile1

upper_basic1 <-
  2 * fit1_weighted$effect_function -
  lower_percentile1

simultaneous_band1 <-
  simultaneous_bootstrap_band(
    estimate =
      fit1_weighted$effect_function,
    bootstrap_estimates =
      bootstrap1,
    alpha = 0.05,
    studentized = TRUE
  )


# ------------------------------------------------------------
# Effect-function plot
# ------------------------------------------------------------

grDevices::png(
  file.path(
    analysis1_directory,
    "03_effect_function.png"
  ),
  width = 1600,
  height = 1050,
  res = 160
)

graphics::plot(
  study1$grid,
  fit1_weighted$effect_function,
  type = "n",
  ylim = range(
    lower_basic1,
    upper_basic1,
    simultaneous_band1$lower,
    simultaneous_band1$upper,
    fit1_unweighted$effect_function,
    0,
    finite = TRUE
  ),
  xlab = "BMI exposure age",
  ylab = expression(
    hat(mu)(t)~
      "log-odds per one-unit centred BMI"
  ),
  main = paste(
    "Functional BMI effect on incident diabetes",
    "with pointwise and simultaneous inference"
  )
)

# Simultaneous band.
graphics::polygon(
  c(
    study1$grid,
    rev(study1$grid)
  ),
  c(
    simultaneous_band1$lower,
    rev(simultaneous_band1$upper)
  ),
  border = NA,
  col = grDevices::rgb(
    0.2,
    0.4,
    0.9,
    0.13
  )
)

# Pointwise band.
graphics::polygon(
  c(
    study1$grid,
    rev(study1$grid)
  ),
  c(
    lower_basic1,
    rev(upper_basic1)
  ),
  border = NA,
  col = grDevices::rgb(
    0.2,
    0.4,
    0.9,
    0.25
  )
)

graphics::lines(
  study1$grid,
  fit1_weighted$effect_function,
  col = "navy",
  lwd = 3
)

graphics::lines(
  study1$grid,
  fit1_unweighted$effect_function,
  col = "firebrick",
  lwd = 2,
  lty = 2
)

graphics::abline(
  h = 0,
  lty = 3
)

graphics::legend(
  "topright",
  legend = c(
    "FPS weighted",
    "Unweighted",
    "95% pointwise band",
    "95% simultaneous band"
  ),
  col = c(
    "navy",
    "firebrick",
    grDevices::rgb(
      0.2,
      0.4,
      0.9,
      0.50
    ),
    grDevices::rgb(
      0.2,
      0.4,
      0.9,
      0.25
    )
  ),
  lty = c(
    1,
    2,
    1,
    1
  ),
  lwd = c(
    3,
    2,
    8,
    8
  ),
  bty = "n"
)

grDevices::dev.off()


# ------------------------------------------------------------
# Save numerical results
# ------------------------------------------------------------

results1 <- data.frame(
  age = study1$grid,
  mu_weighted =
    fit1_weighted$effect_function,
  mu_unweighted =
    fit1_unweighted$effect_function,
  pointwise_lower =
    lower_basic1,
  pointwise_upper =
    upper_basic1,
  simultaneous_lower =
    simultaneous_band1$lower,
  simultaneous_upper =
    simultaneous_band1$upper,
  pointwise_excludes_zero =
    as.integer(
      lower_basic1 > 0 |
        upper_basic1 < 0
    ),
  simultaneous_excludes_zero =
    as.integer(
      simultaneous_band1$lower > 0 |
        simultaneous_band1$upper < 0
    )
)

utils::write.csv(
  results1,
  file.path(
    analysis1_directory,
    "effect_function.csv"
  ),
  row.names = FALSE
)

coefficient_names <- names(
  fit1_weighted$coefficients
)

utils::write.csv(
  data.frame(
    term = coefficient_names,
    weighted = as.numeric(
      fit1_weighted$coefficients
    ),
    unweighted = as.numeric(
      fit1_unweighted$coefficients
    )
  ),
  file.path(
    analysis1_directory,
    "score_coefficients.csv"
  ),
  row.names = FALSE
)

writeLines(
  c(
    paste(
      "n:",
      length(Y1)
    ),
    paste(
      "incident diabetes events:",
      sum(Y1)
    ),
    paste(
      "treatment FPCs in balance model:",
      L_balance1
    ),
    paste(
      "treatment FPCs in outcome model:",
      fit1_weighted$Lstar
    ),
    paste(
      "confounder design columns:",
      ncol(confounders1$C)
    ),
    paste(
      "nominal balancing constraints:",
      L_balance1 +
        ncol(confounders1$C) +
        L_balance1 *
        ncol(confounders1$C)
    ),
    paste(
      "ESS:",
      round(
        diagnostics1$ESS,
        1
      )
    ),
    paste(
      "ESS fraction:",
      round(
        diagnostics1$ESS_fraction,
        4
      )
    ),
    paste(
      "maximum absolute correlation before:",
      round(
        balance1$maximum_before,
        4
      )
    ),
    paste(
      "maximum absolute correlation after:",
      round(
        balance1$maximum_after,
        4
      )
    ),
    paste(
      "simultaneous critical value:",
      round(
        simultaneous_band1$critical_value,
        3
      )
    ),
    paste(
      "valid simultaneous bootstrap replicates:",
      simultaneous_band1$
        valid_bootstrap_replicates
    ),
    paste(
      "failed simultaneous bootstrap replicates:",
      simultaneous_band1$
        failed_bootstrap_replicates
    ),
    "",
    paste(
      "The effect function is on the marginal",
      "log-odds scale."
    ),
    paste(
      "It represents the effect per one kg/m^2",
      "deviation from the cohort mean BMI curve",
      "at exposure age t."
    )
  ),
  file.path(
    analysis1_directory,
    "summary.txt"
  )
)

log_line(
  "Study 1 analysis complete."
)
