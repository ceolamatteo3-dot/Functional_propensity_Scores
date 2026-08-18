# ============================================================
# PROJECT CONFIGURATION
# ============================================================

CFG <- list(
  
  # RAND HRS Stata dataset.
  data_file = file.path(
    "data",
    "randhrs1992_2022v1.dta"
  ),
  
  out_root = "outputs",
  
  # Eligible entry cohorts:
  # 3 = original HRS
  # 4 = War Babies
  # 5 = Early Baby Boomers
  # 6 = Mid Baby Boomers
  cohorts = c(3, 4, 5, 6),
  
  # ----------------------------------------------------------
  # Study 1
  # Functional BMI exposure followed by incident diabetes.
  # ----------------------------------------------------------
  
  s1_expo_lo = 52,
  s1_expo_hi = 64,
  s1_fup_hi  = 76,
  
  # ----------------------------------------------------------
  # Study 2
  # Functional BMI exposure and later functional outcome.
  # ----------------------------------------------------------
  
  s2_expo_lo = 52,
  s2_expo_hi = 62,
  s2_out_lo  = 64,
  s2_out_hi  = 74,
  
  # Baseline is the latest observed interview at or below
  # this age.
  baseline_age_max = 54,
  
  # Trajectory reconstruction requirements.
  min_expo_obs      = 4,
  min_out_obs       = 3,
  min_span_fraction = 0.70,
  n_grid            = 41,
  
  # Options: "spline3", "quadratic", "linear".
  smooth_method = "spline3",
  
  # Plausible BMI range.
  bmi_lo = 12,
  bmi_hi = 75,
  
  # Study 2 outcome.
  # Options: "cog27", "cesd", "adl5a".
  outcome2 = "cog27",
  
  # FPS settings.
  pve_balance = 0.95,
  pve_outcome = 0.95,
  
  # Options: "core", "extended".
  conf_set = "core",
  
  # Zero gives exact balance.
  # Positive values provide ridge-relaxed approximate balance.
  fps_ridge = 0,
  
  min_ess_frac = 0.10,
  
  # Use 25-50 for a test run, then 1000 for final inference.
  n_boot = 21,
  
  seed = 20260818
)

CFG$dirs <- file.path(
  CFG$out_root,
  c(
    "eda",
    "analysis1",
    "analysis2",
    "tables",
    "logs"
  )
)

for (directory in c(CFG$out_root, CFG$dirs)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

log_line <- function(...) {
  message_text <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "  ",
    paste0(..., collapse = "")
  )
  
  cat(message_text, "\n")
  
  cat(
    message_text,
    "\n",
    file = file.path(
      CFG$out_root,
      "logs",
      "run_log.txt"
    ),
    append = TRUE
  )
  
  invisible(NULL)
}
