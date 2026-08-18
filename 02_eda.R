# ============================================================
# EXPLORATORY ANALYSIS
# ============================================================

eda_directory <- file.path(
  CFG$out_root,
  "eda"
)

open_png <- function(filename,
                     width = 1600,
                     height = 1100) {
  grDevices::png(
    file.path(
      eda_directory,
      filename
    ),
    width = width,
    height = height,
    res = 160
  )
}

confounder_names <- confounder_formula(
  CFG$conf_set
)


# ------------------------------------------------------------
# Exposure observation counts
# ------------------------------------------------------------

observation_counts <- long[
  hhidpn %in% study1$dat$hhidpn &
    age >= CFG$s1_expo_lo &
    age <= CFG$s1_expo_hi &
    is.finite(bmi),
  .(
    n_bmi_observations = uniqueN(age),
    minimum_age = min(age),
    maximum_age = max(age),
    observed_span = max(age) - min(age)
  ),
  by = hhidpn
]

utils::write.csv(
  observation_counts,
  file.path(
    CFG$out_root,
    "tables",
    "study1_exposure_observation_counts.csv"
  ),
  row.names = FALSE
)

open_png("01_observations_per_subject.png")

graphics::hist(
  observation_counts$n_bmi_observations,
  breaks = seq(
    0.5,
    max(
      observation_counts$n_bmi_observations
    ) + 0.5,
    by = 1
  ),
  col = "grey80",
  xlab = "BMI observations in exposure window",
  ylab = "Subjects",
  main = "Sampling density of functional BMI exposure"
)

grDevices::dev.off()


# ------------------------------------------------------------
# Reconstructed BMI curves
# ------------------------------------------------------------

open_png("02_bmi_trajectories.png")

graphics::par(
  mfrow = c(1, 2)
)

graphics::plot(
  NA,
  xlim = range(study1$grid),
  ylim = range(study1$X_raw),
  xlab = "Age",
  ylab = expression(BMI~(kg/m^2)),
  main = "Reconstructed BMI trajectories"
)

set.seed(CFG$seed)

sampled_rows <- sample(
  seq_len(nrow(study1$X_raw)),
  min(
    250,
    nrow(study1$X_raw)
  )
)

for (i in sampled_rows) {
  graphics::lines(
    study1$grid,
    study1$X_raw[i, ],
    col = grDevices::rgb(
      0.4,
      0.4,
      0.4,
      0.15
    )
  )
}

graphics::lines(
  study1$grid,
  study1$mean_curve,
  col = "firebrick",
  lwd = 3
)

graphics::plot(
  NA,
  xlim = range(study1$grid),
  ylim = range(study1$X),
  xlab = "Age",
  ylab = "Centred BMI",
  main = "Pointwise-centred model input"
)

for (i in sampled_rows) {
  graphics::lines(
    study1$grid,
    study1$X[i, ],
    col = grDevices::rgb(
      0.2,
      0.3,
      0.7,
      0.15
    )
  )
}

graphics::abline(
  h = 0,
  lty = 2
)

grDevices::dev.off()


# ------------------------------------------------------------
# FPCA
# ------------------------------------------------------------

fpca1_eda <- dense_fpca(
  curves = study1$X,
  grid = study1$grid,
  pve_keep = 0.999
)

number_to_plot <- min(
  3,
  ncol(fpca1_eda$eigenfunctions)
)

open_png("03_exposure_fpca.png")

graphics::par(
  mfrow = c(1, 3)
)

graphics::plot(
  seq_along(fpca1_eda$eigenvalues),
  fpca1_eda$pve,
  type = "b",
  pch = 19,
  xlab = "Component",
  ylab = "Cumulative PVE",
  main = "Exposure variance explained",
  ylim = c(0, 1)
)

graphics::abline(
  h = c(0.95, 0.99),
  lty = c(2, 3),
  col = "red"
)

graphics::matplot(
  fpca1_eda$grid,
  fpca1_eda$eigenfunctions[
    ,
    seq_len(number_to_plot),
    drop = FALSE
  ],
  type = "l",
  lty = 1,
  lwd = 2,
  xlab = "Age",
  ylab = expression(phi[k](t)),
  main = "Leading eigenfunctions"
)

graphics::abline(
  h = 0,
  lty = 3
)

graphics::legend(
  "topright",
  paste0(
    "PC",
    seq_len(number_to_plot)
  ),
  col = seq_len(number_to_plot),
  lty = 1,
  lwd = 2,
  bty = "n"
)

second_score <- if (
  ncol(fpca1_eda$scores) >= 2
) {
  fpca1_eda$scores[, 2]
} else {
  rep(
    0,
    nrow(fpca1_eda$scores)
  )
}

graphics::plot(
  fpca1_eda$scores[, 1],
  second_score,
  pch = 19,
  col = grDevices::rgb(
    0.1,
    0.3,
    0.8,
    0.4
  ),
  xlab = "FPC 1 score",
  ylab = "FPC 2 score",
  main = "Exposure-score scatter"
)

grDevices::dev.off()


# ------------------------------------------------------------
# Perturbation interpretation
# ------------------------------------------------------------

open_png("04_fpc_perturbations.png")

graphics::par(
  mfrow = c(
    1,
    number_to_plot
  )
)

for (k in seq_len(number_to_plot)) {
  score_sd <- sqrt(
    fpca1_eda$eigenvalues[k]
  )
  
  plus_curve <-
    study1$mean_curve +
    2 * score_sd *
    fpca1_eda$eigenfunctions[, k]
  
  minus_curve <-
    study1$mean_curve -
    2 * score_sd *
    fpca1_eda$eigenfunctions[, k]
  
  graphics::plot(
    study1$grid,
    study1$mean_curve,
    type = "l",
    lwd = 2,
    ylim = range(
      study1$mean_curve,
      plus_curve,
      minus_curve
    ),
    xlab = "Age",
    ylab = "BMI",
    main = paste0(
      "PC",
      k,
      ": ",
      round(
        100 *
          fpca1_eda$eigenvalues[k] /
          fpca1_eda$total_variance,
        1
      ),
      "%"
    )
  )
  
  graphics::lines(
    study1$grid,
    plus_curve,
    col = "firebrick",
    lty = 2,
    lwd = 2
  )
  
  graphics::lines(
    study1$grid,
    minus_curve,
    col = "steelblue",
    lty = 2,
    lwd = 2
  )
}

grDevices::dev.off()


# ------------------------------------------------------------
# Confounding correlation heatmap
# ------------------------------------------------------------

L_balance_eda <- select_ncomp(
  fpca1_eda,
  CFG$pve_balance
)

A_eda <- fpca1_eda$scores[
  ,
  seq_len(L_balance_eda),
  drop = FALSE
]

C_eda <- prepare_confounders(
  scalar_covariates = as.data.frame(
    study1$dat[
      ,
      ..confounder_names
    ]
  )
)$C

pre_correlations <- stats::cor(
  A_eda,
  C_eda
)

utils::write.csv(
  round(
    pre_correlations,
    4
  ),
  file.path(
    CFG$out_root,
    "tables",
    "eda_pre_weight_correlations.csv"
  )
)

open_png(
  "05_confounding_heatmap.png",
  width = 1800,
  height = 900
)

graphics::image(
  seq_len(ncol(C_eda)),
  seq_len(nrow(pre_correlations)),
  t(abs(pre_correlations)),
  col = grDevices::hcl.colors(
    24,
    "YlOrRd",
    rev = TRUE
  ),
  xlab = "",
  ylab = "Treatment FPC",
  axes = FALSE,
  main = paste(
    "Absolute treatment-score/confounder correlation",
    "before weighting"
  )
)

graphics::axis(
  1,
  at = seq_len(ncol(C_eda)),
  labels = colnames(C_eda),
  las = 2,
  cex.axis = 0.6
)

graphics::axis(
  2,
  at = seq_len(nrow(pre_correlations)),
  labels = rownames(pre_correlations),
  las = 1
)

graphics::box()

grDevices::dev.off()


# ------------------------------------------------------------
# Outcome descriptions
# ------------------------------------------------------------

open_png(
  "06_outcome_descriptives.png",
  width = 1800,
  height = 900
)

graphics::par(
  mfrow = c(1, 3)
)

graphics::barplot(
  table(study1$dat$Y),
  names.arg = c(
    "No diabetes",
    "Incident diabetes"
  ),
  col = c(
    "grey80",
    "firebrick"
  ),
  ylab = "Subjects",
  main = "Study 1 outcome"
)

graphics::boxplot(
  fpca1_eda$scores[, 1] ~
    study1$dat$Y,
  names = c(
    "No",
    "Yes"
  ),
  col = c(
    "grey85",
    "lightpink"
  ),
  xlab = "Incident diabetes",
  ylab = "Exposure FPC 1",
  main = "Exposure level by outcome"
)

graphics::plot(
  NA,
  xlim = range(study2$grid_t),
  ylim = range(study2$Y),
  xlab = "Age",
  ylab = study2$outcome_name,
  main = paste(
    "Study 2 outcome:",
    study2$outcome_name
  )
)

sampled_outcomes <- sample(
  seq_len(nrow(study2$Y)),
  min(
    200,
    nrow(study2$Y)
  )
)

for (i in sampled_outcomes) {
  graphics::lines(
    study2$grid_t,
    study2$Y[i, ],
    col = grDevices::rgb(
      0.1,
      0.4,
      0.6,
      0.15
    )
  )
}

graphics::lines(
  study2$grid_t,
  colMeans(study2$Y),
  col = "navy",
  lwd = 3
)

grDevices::dev.off()


# ------------------------------------------------------------
# Baseline descriptive tables
# ------------------------------------------------------------

baseline_table <- function(dat, variables) {
  rows <- lapply(
    variables,
    function(variable) {
      x <- dat[[variable]]
      
      if (
        is.factor(x) ||
        length(unique(stats::na.omit(x))) <= 2
      ) {
        tabulation <- table(
          x,
          useNA = "ifany"
        )
        
        data.frame(
          variable = variable,
          level = names(tabulation),
          summary = paste0(
            as.integer(tabulation),
            " (",
            round(
              100 *
                as.integer(tabulation) /
                length(x),
              1
            ),
            "%)"
          )
        )
      } else {
        data.frame(
          variable = variable,
          level = "mean (sd)",
          summary = paste0(
            round(
              mean(x, na.rm = TRUE),
              2
            ),
            " (",
            round(
              stats::sd(
                x,
                na.rm = TRUE
              ),
              2
            ),
            ")"
          )
        )
      }
    }
  )
  
  do.call(
    rbind,
    rows
  )
}

utils::write.csv(
  baseline_table(
    study1$dat,
    c(
      "baseline_age",
      confounder_names,
      "hacohort"
    )
  ),
  file.path(
    CFG$out_root,
    "tables",
    "table1_study1.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  baseline_table(
    study2$dat,
    c(
      "baseline_age",
      confounder_names,
      "hacohort"
    )
  ),
  file.path(
    CFG$out_root,
    "tables",
    "table1_study2.csv"
  ),
  row.names = FALSE
)

log_line(
  "EDA complete. Exposure FPCs retained at PVE ",
  CFG$pve_balance,
  ": ",
  L_balance_eda
)
