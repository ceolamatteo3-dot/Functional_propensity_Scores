# ============================================================
# READ RAND HRS .DTA AND BUILD ANALYTIC COHORTS
# ============================================================

suppressPackageStartupMessages({
  library(haven)
  library(data.table)
  library(splines)
})

WAVES <- 1:16


# ============================================================
# ROBUST STATA IMPORT
# ============================================================

normalize_stata_column <- function(x) {
  # Convert Stata special/tagged missing values to ordinary NA.
  x <- haven::zap_missing(x)
  
  # Remove haven_labelled/value-label classes. This leaves the
  # underlying numeric codes used in the RAND documentation.
  if (inherits(x, "haven_labelled")) {
    x <- haven::zap_labels(x)
  }
  
  x
}


read_rand_hrs_dta <- function(path,
                              requested_names) {
  if (!file.exists(path)) {
    stop(
      "RAND HRS .dta file not found: ",
      normalizePath(
        path,
        winslash = "/",
        mustWork = FALSE
      )
    )
  }
  
  log_line(
    "Reading Stata metadata from: ",
    path
  )
  
  # Read metadata/header only.
  header <- haven::read_dta(
    path,
    n_max = 0,
    .name_repair = "minimal"
  )
  
  original_names <- names(header)
  lower_names <- tolower(original_names)
  
  if (anyDuplicated(lower_names)) {
    duplicated_names <- unique(
      lower_names[
        duplicated(lower_names)
      ]
    )
    
    stop(
      "Variable names are not unique after conversion ",
      "to lowercase: ",
      paste(
        duplicated_names,
        collapse = ", "
      )
    )
  }
  
  # Map lowercase names used by the R scripts to the original
  # names stored in the Stata file.
  name_map <- stats::setNames(
    original_names,
    lower_names
  )
  
  requested_lower <- unique(
    tolower(requested_names)
  )
  
  available_lower <- intersect(
    requested_lower,
    lower_names
  )
  
  missing_lower <- setdiff(
    requested_lower,
    lower_names
  )
  
  original_keep <- unname(
    name_map[available_lower]
  )
  
  writeLines(
    missing_lower,
    file.path(
      CFG$out_root,
      "logs",
      "variables_not_in_file.txt"
    )
  )
  
  data.frame(
    requested_name = requested_lower,
    available = requested_lower %in% lower_names,
    original_stata_name = unname(
      name_map[requested_lower]
    ),
    stringsAsFactors = FALSE
  ) |>
    utils::write.csv(
      file.path(
        CFG$out_root,
        "logs",
        "variable_import_map.csv"
      ),
      row.names = FALSE
    )
  
  log_line(
    "Variables in Stata file: ",
    length(original_names)
  )
  
  log_line(
    "Requested variables: ",
    length(requested_lower),
    "; available: ",
    length(available_lower),
    "; missing: ",
    length(missing_lower)
  )
  
  if (length(original_keep) == 0L) {
    stop(
      "None of the requested variables are in the .dta file."
    )
  }
  
  # Only the required columns are loaded from the very wide
  # Stata dataset.
  imported <- haven::read_dta(
    path,
    col_select = tidyselect::all_of(
      original_keep
    ),
    .name_repair = "minimal"
  )
  
  imported <- data.table::as.data.table(
    imported
  )
  
  # All downstream code uses lowercase RAND names.
  data.table::setnames(
    imported,
    tolower(names(imported))
  )
  
  for (variable in names(imported)) {
    data.table::set(
      imported,
      j = variable,
      value = normalize_stata_column(
        imported[[variable]]
      )
    )
  }
  
  imported
}


# ============================================================
# REQUESTED VARIABLES
# ============================================================

wave_variables <- c(
  "agey_e",
  "iwstat",
  "bmi",
  "wtresp",
  "mstat",
  "shlt",
  "diab",
  "diabe",
  "hibpe",
  "hearte",
  "stroke",
  "lunge",
  "cancre",
  "arthre",
  "psyche",
  "smokev",
  "smoken",
  "drinkd",
  "drinkr",
  "cesd",
  "cog27",
  "adl5a",
  "higov"
)

household_variables <- c(
  "atotb",
  "itot"
)

fixed_variables <- c(
  "hhidpn",
  "hhid",
  "pn",
  "hacohort",
  "racohbyr",
  "ragender",
  "raracem",
  "rahispan",
  "raeduc",
  "raedyrs",
  "rabyear",
  "radyear",
  "raestrat",
  "raehsamp"
)

respondent_wave_variables <- unlist(
  lapply(
    WAVES,
    function(wave) {
      paste0(
        "r",
        wave,
        wave_variables
      )
    }
  ),
  use.names = FALSE
)

household_wave_variables <- unlist(
  lapply(
    WAVES,
    function(wave) {
      paste0(
        "h",
        wave,
        household_variables
      )
    }
  ),
  use.names = FALSE
)

desired_variables <- unique(
  c(
    fixed_variables,
    respondent_wave_variables,
    household_wave_variables
  )
)


wide <- read_rand_hrs_dta(
  path = CFG$data_file,
  requested_names = desired_variables
)

log_line(
  "Rows imported from .dta file: ",
  nrow(wide)
)


# ============================================================
# REQUIRED VARIABLE VALIDATION
# ============================================================

required_fixed <- c(
  "hhidpn",
  "hhid",
  "hacohort",
  "ragender",
  "raracem",
  "rahispan",
  "raeduc",
  "rabyear"
)

missing_required_fixed <- setdiff(
  required_fixed,
  names(wide)
)

if (length(missing_required_fixed) > 0L) {
  stop(
    "Required fixed variables are missing: ",
    paste(
      missing_required_fixed,
      collapse = ", "
    )
  )
}

required_trajectory_patterns <- c(
  "^r[0-9]+agey_e$",
  "^r[0-9]+iwstat$",
  "^r[0-9]+bmi$"
)

for (pattern in required_trajectory_patterns) {
  if (!any(grepl(pattern, names(wide)))) {
    stop(
      "No variables match required pattern: ",
      pattern
    )
  }
}


# ============================================================
# SAFE VARIABLE ACCESSORS
# ============================================================

get_numeric_wave <- function(prefix,
                             wave,
                             suffix) {
  variable <- paste0(
    prefix,
    wave,
    suffix
  )
  
  if (!variable %in% names(wide)) {
    return(
      rep(
        NA_real_,
        nrow(wide)
      )
    )
  }
  
  value <- wide[[variable]]
  
  if (is.factor(value)) {
    value <- as.character(value)
  }
  
  suppressWarnings(
    as.numeric(value)
  )
}


get_fixed <- function(variable,
                      default = NA_real_) {
  if (variable %in% names(wide)) {
    return(wide[[variable]])
  }
  
  rep(
    default,
    nrow(wide)
  )
}


# ============================================================
# RESHAPE PERSON-WAVE VARIABLES TO LONG FORMAT
# ============================================================

long_list <- lapply(
  WAVES,
  function(wave) {
    data.table(
      hhidpn = get_fixed("hhidpn"),
      hhid = as.character(
        get_fixed("hhid", NA_character_)
      ),
      wave = wave,
      
      age = get_numeric_wave(
        "r",
        wave,
        "agey_e"
      ),
      
      iwstat = get_numeric_wave(
        "r",
        wave,
        "iwstat"
      ),
      
      bmi = get_numeric_wave(
        "r",
        wave,
        "bmi"
      ),
      
      svywt = get_numeric_wave(
        "r",
        wave,
        "wtresp"
      ),
      
      shlt = get_numeric_wave(
        "r",
        wave,
        "shlt"
      ),
      
      diab = get_numeric_wave(
        "r",
        wave,
        "diab"
      ),
      
      diabe = get_numeric_wave(
        "r",
        wave,
        "diabe"
      ),
      
      hibpe = get_numeric_wave(
        "r",
        wave,
        "hibpe"
      ),
      
      hearte = get_numeric_wave(
        "r",
        wave,
        "hearte"
      ),
      
      stroke = get_numeric_wave(
        "r",
        wave,
        "stroke"
      ),
      
      lunge = get_numeric_wave(
        "r",
        wave,
        "lunge"
      ),
      
      cancre = get_numeric_wave(
        "r",
        wave,
        "cancre"
      ),
      
      arthre = get_numeric_wave(
        "r",
        wave,
        "arthre"
      ),
      
      psyche = get_numeric_wave(
        "r",
        wave,
        "psyche"
      ),
      
      smokev = get_numeric_wave(
        "r",
        wave,
        "smokev"
      ),
      
      smoken = get_numeric_wave(
        "r",
        wave,
        "smoken"
      ),
      
      drinkd = get_numeric_wave(
        "r",
        wave,
        "drinkd"
      ),
      
      drinkr = get_numeric_wave(
        "r",
        wave,
        "drinkr"
      ),
      
      cesd = get_numeric_wave(
        "r",
        wave,
        "cesd"
      ),
      
      cog27 = get_numeric_wave(
        "r",
        wave,
        "cog27"
      ),
      
      adl5a = get_numeric_wave(
        "r",
        wave,
        "adl5a"
      ),
      
      higov = get_numeric_wave(
        "r",
        wave,
        "higov"
      ),
      
      atotb = get_numeric_wave(
        "h",
        wave,
        "atotb"
      ),
      
      itot = get_numeric_wave(
        "h",
        wave,
        "itot"
      )
    )
  }
)

long <- data.table::rbindlist(
  long_list,
  use.names = TRUE,
  fill = TRUE
)

# A live Core Interview is code 1.
long <- long[
  iwstat == 1 &
    is.finite(age)
]

long[
  !is.finite(bmi) |
    bmi < CFG$bmi_lo |
    bmi > CFG$bmi_hi,
  bmi := NA_real_
]

data.table::setorder(
  long,
  hhidpn,
  wave
)

log_line(
  "Live person-wave records: ",
  nrow(long)
)


# ============================================================
# PERSON-LEVEL FIXED CHARACTERISTICS
# ============================================================

person <- data.table(
  hhidpn = get_fixed("hhidpn"),
  hhid = as.character(
    get_fixed("hhid", NA_character_)
  ),
  hacohort = suppressWarnings(
    as.numeric(get_fixed("hacohort"))
  ),
  ragender = suppressWarnings(
    as.numeric(get_fixed("ragender"))
  ),
  raracem = suppressWarnings(
    as.numeric(get_fixed("raracem"))
  ),
  rahispan = suppressWarnings(
    as.numeric(get_fixed("rahispan"))
  ),
  raeduc = suppressWarnings(
    as.numeric(get_fixed("raeduc"))
  ),
  raedyrs = suppressWarnings(
    as.numeric(get_fixed("raedyrs"))
  ),
  rabyear = suppressWarnings(
    as.numeric(get_fixed("rabyear"))
  ),
  radyear = suppressWarnings(
    as.numeric(get_fixed("radyear"))
  ),
  raestrat = suppressWarnings(
    as.numeric(get_fixed("raestrat"))
  ),
  raehsamp = suppressWarnings(
    as.numeric(get_fixed("raehsamp"))
  )
)

# Use HHIDPN as a fallback cluster when HHID is unavailable.
person[
  is.na(hhid) |
    hhid == "",
  hhid := paste0(
    "person_",
    format(
      hhidpn,
      scientific = FALSE,
      trim = TRUE
    )
  )
]

person[, age_at_death := ifelse(
  is.finite(radyear) &
    is.finite(rabyear),
  radyear - rabyear,
  NA_real_
)]


# ============================================================
# BASELINE CONFOUNDER CONSTRUCTION
# ============================================================

build_baseline <- function(long_data,
                           person_data,
                           age_max) {
  base <- long_data[
    age <= age_max
  ]
  
  data.table::setorder(
    base,
    hhidpn,
    -age
  )
  
  # Latest live interview at or below baseline age.
  base <- base[
    ,
    .SD[1],
    by = hhidpn
  ]
  
  out <- merge(
    person_data,
    base,
    by = c("hhidpn", "hhid"),
    all = FALSE
  )
  
  out[, sex := factor(
    ifelse(
      ragender == 1,
      "male",
      "female"
    ),
    levels = c(
      "male",
      "female"
    )
  )]
  
  out[, race3 := factor(
    ifelse(
      raracem == 1,
      "white",
      ifelse(
        raracem == 2,
        "black",
        "other"
      )
    ),
    levels = c(
      "white",
      "black",
      "other"
    )
  )]
  
  out[, hisp := factor(
    ifelse(
      is.finite(rahispan) &
        rahispan > 0,
      "yes",
      "no"
    ),
    levels = c(
      "no",
      "yes"
    )
  )]
  
  out[, educ4 := factor(
    ifelse(
      raeduc %in% c(1, 2),
      "lt_hs",
      ifelse(
        raeduc == 3,
        "hs_ged",
        ifelse(
          raeduc == 4,
          "some_col",
          "college"
        )
      )
    ),
    levels = c(
      "lt_hs",
      "hs_ged",
      "some_col",
      "college"
    )
  )]
  
  out[, smoke3 := factor(
    ifelse(
      is.finite(smoken) &
        smoken == 1,
      "current",
      ifelse(
        is.finite(smokev) &
          smokev == 1,
        "former",
        "never"
      )
    ),
    levels = c(
      "never",
      "former",
      "current"
    )
  )]
  
  out[, drink_any := factor(
    ifelse(
      is.finite(drinkd),
      ifelse(
        drinkd > 0,
        "yes",
        "no"
      ),
      ifelse(
        is.finite(drinkr),
        ifelse(
          drinkr > 0,
          "yes",
          "no"
        ),
        NA_character_
      )
    ),
    levels = c(
      "no",
      "yes"
    )
  )]
  
  binary_conditions <- c(
    "hibpe",
    "hearte",
    "stroke",
    "lunge",
    "cancre",
    "arthre",
    "psyche"
  )
  
  for (variable in binary_conditions) {
    out[
      ,
      (variable) := as.numeric(
        get(variable) == 1
      )
    ]
  }
  
  out[, shlt_poor := as.numeric(
    shlt >= 4
  )]
  
  out[, cesd_base := cesd]
  
  out[, log_wealth := log1p(
    pmax(atotb, 0)
  )]
  
  out[, log_income := log1p(
    pmax(itot, 0)
  )]
  
  out[, baseline_age := age]
  out[, baseline_wave := wave]
  out[, svywt_base := svywt]
  
  out
}


baseline <- build_baseline(
  long_data = long,
  person_data = person,
  age_max = CFG$baseline_age_max
)

log_line(
  "Subjects with baseline interview at age <= ",
  CFG$baseline_age_max,
  ": ",
  nrow(baseline)
)


confounder_formula <- function(which_set) {
  core <- c(
    "sex",
    "race3",
    "hisp",
    "educ4",
    "smoke3",
    "hibpe",
    "hearte",
    "shlt_poor",
    "log_wealth"
  )
  
  extended <- c(
    core,
    "drink_any",
    "stroke",
    "lunge",
    "cancre",
    "arthre",
    "psyche",
    "cesd_base",
    "log_income"
  )
  
  if (which_set == "core") {
    core
  } else if (which_set == "extended") {
    extended
  } else {
    stop(
      "Unknown confounder set: ",
      which_set
    )
  }
}


# ============================================================
# TRAJECTORY RECONSTRUCTION
# ============================================================

fit_one_curve <- function(age,
                          outcome,
                          grid,
                          method) {
  ok <- is.finite(age) &
    is.finite(outcome)
  
  age <- age[ok]
  outcome <- outcome[ok]
  
  # Collapse multiple observations at the same reported age.
  collapsed <- stats::aggregate(
    outcome,
    by = list(age = age),
    FUN = mean
  )
  
  age <- collapsed$age
  outcome <- collapsed$x
  
  n_obs <- length(age)
  
  if (n_obs < 2L) {
    return(
      rep(
        NA_real_,
        length(grid)
      )
    )
  }
  
  method_effective <- method
  
  if (
    method_effective == "spline3" &&
    n_obs < 5L
  ) {
    method_effective <- "quadratic"
  }
  
  if (
    method_effective == "quadratic" &&
    n_obs < 3L
  ) {
    method_effective <- "linear"
  }
  
  model_data <- data.frame(
    age = age,
    outcome = outcome
  )
  
  fit <- switch(
    method_effective,
    
    spline3 = stats::lm(
      outcome ~ splines::ns(age, df = 3),
      data = model_data
    ),
    
    quadratic = stats::lm(
      outcome ~ stats::poly(
        age,
        2,
        raw = TRUE
      ),
      data = model_data
    ),
    
    linear = stats::lm(
      outcome ~ age,
      data = model_data
    ),
    
    stop(
      "Unknown smoothing method: ",
      method_effective
    )
  )
  
  # Avoid unsupported spline extrapolation. Predictions outside
  # a subject's observed range are held at the nearest boundary.
  prediction_age <- pmin(
    pmax(grid, min(age)),
    max(age)
  )
  
  as.numeric(
    stats::predict(
      fit,
      newdata = data.frame(
        age = prediction_age
      )
    )
  )
}


build_curves <- function(long_data,
                         variable,
                         lo,
                         hi,
                         min_obs,
                         min_span,
                         n_grid,
                         method) {
  if (!variable %in% names(long_data)) {
    stop(
      "Trajectory variable is unavailable: ",
      variable
    )
  }
  
  grid <- seq(
    lo,
    hi,
    length.out = n_grid
  )
  
  dt <- long_data[
    age >= lo &
      age <= hi &
      is.finite(get(variable))
  ]
  
  if (nrow(dt) == 0L) {
    stop(
      "No usable observations for ",
      variable,
      " in the requested age window."
    )
  }
  
  summary_data <- dt[
    ,
    .(
      n_obs = uniqueN(age),
      span = max(age) - min(age),
      minimum_age = min(age),
      maximum_age = max(age),
      lo_ok =
        min(age) <=
        lo + 0.35 * (hi - lo),
      hi_ok =
        max(age) >=
        hi - 0.35 * (hi - lo)
    ),
    by = hhidpn
  ]
  
  eligible_ids <- summary_data[
    n_obs >= min_obs &
      span >= min_span * (hi - lo) &
      lo_ok &
      hi_ok,
    hhidpn
  ]
  
  dt <- dt[
    hhidpn %in% eligible_ids
  ]
  
  ids <- sort(
    unique(dt$hhidpn)
  )
  
  X <- matrix(
    NA_real_,
    nrow = length(ids),
    ncol = n_grid,
    dimnames = list(
      as.character(ids),
      NULL
    )
  )
  
  for (i in seq_along(ids)) {
    subject <- dt[
      hhidpn == ids[i]
    ]
    
    X[i, ] <- fit_one_curve(
      age = subject$age,
      outcome = subject[[variable]],
      grid = grid,
      method = method
    )
  }
  
  complete <- apply(
    X,
    1,
    function(row) {
      all(is.finite(row))
    }
  )
  
  list(
    X = X[
      complete,
      ,
      drop = FALSE
    ],
    ids = ids[complete],
    grid = grid,
    diagnostics = summary_data
  )
}


# ============================================================
# STUDY 1 COHORT
# ============================================================

build_study1 <- function() {
  ids0 <- baseline[
    hacohort %in% CFG$cohorts,
    hhidpn
  ]
  
  log_line(
    "S1 step 1, eligible cohorts: ",
    length(ids0)
  )
  
  L <- long[
    hhidpn %in% ids0
  ]
  
  # Diabetes-free through the end of the exposure window.
  pre <- L[
    age <= CFG$s1_expo_hi,
    .(
      ever_dm_pre = any(
        diab == 1 |
          diabe == 1,
        na.rm = TRUE
      ),
      any_dm_information = any(
        is.finite(diab) |
          is.finite(diabe)
      )
    ),
    by = hhidpn
  ]
  
  ids1 <- pre[
    !ever_dm_pre &
      any_dm_information,
    hhidpn
  ]
  
  log_line(
    "S1 step 2, diabetes-free through age ",
    CFG$s1_expo_hi,
    ": ",
    length(ids1)
  )
  
  post <- L[
    hhidpn %in% ids1 &
      age > CFG$s1_expo_hi &
      age <= CFG$s1_fup_hi
  ]
  
  post_summary <- post[
    ,
    .(
      n_post = .N,
      any_dm_information = any(
        is.finite(diab) |
          is.finite(diabe)
      ),
      Y = as.numeric(
        any(
          diab == 1 |
            diabe == 1,
          na.rm = TRUE
        )
      ),
      age_last = max(age)
    ),
    by = hhidpn
  ]
  
  post_summary <- post_summary[
    n_post >= 1 &
      any_dm_information
  ]
  
  ids2 <- post_summary$hhidpn
  
  log_line(
    "S1 step 3, with informative follow-up: ",
    length(ids2)
  )
  
  curves <- build_curves(
    long_data = long[
      hhidpn %in% ids2
    ],
    variable = "bmi",
    lo = CFG$s1_expo_lo,
    hi = CFG$s1_expo_hi,
    min_obs = CFG$min_expo_obs,
    min_span = CFG$min_span_fraction,
    n_grid = CFG$n_grid,
    method = CFG$smooth_method
  )
  
  log_line(
    "S1 step 4, with reconstructable BMI curve: ",
    length(curves$ids)
  )
  
  dat <- merge(
    data.table(
      hhidpn = curves$ids
    ),
    post_summary,
    by = "hhidpn"
  )
  
  dat <- merge(
    dat,
    baseline,
    by = "hhidpn"
  )
  
  keep_id <- dat$hhidpn
  
  X_raw <- curves$X[
    as.character(keep_id),
    ,
    drop = FALSE
  ]
  
  confounders <- confounder_formula(
    CFG$conf_set
  )
  
  complete <- stats::complete.cases(
    dat[
      ,
      ..confounders
    ]
  )
  
  log_line(
    "S1 step 5, complete baseline confounders: ",
    sum(complete)
  )
  
  dat <- dat[complete]
  
  X_raw <- X_raw[
    complete,
    ,
    drop = FALSE
  ]
  
  if (nrow(dat) < 20L) {
    stop(
      "Study 1 has fewer than 20 complete subjects."
    )
  }
  
  if (sum(dat$Y) < 5L) {
    warning(
      "Study 1 has fewer than five incident diabetes events."
    )
  }
  
  mean_curve <- colMeans(X_raw)
  
  X_centered <- sweep(
    X_raw,
    2,
    mean_curve,
    "-"
  )
  
  list(
    dat = dat,
    X = X_centered,
    X_raw = X_raw,
    mean_curve = mean_curve,
    grid = curves$grid,
    trajectory_diagnostics =
      curves$diagnostics
  )
}


# ============================================================
# STUDY 2 COHORT
# ============================================================

build_study2 <- function() {
  ids0 <- baseline[
    hacohort %in% CFG$cohorts,
    hhidpn
  ]
  
  exposure <- build_curves(
    long_data = long[
      hhidpn %in% ids0
    ],
    variable = "bmi",
    lo = CFG$s2_expo_lo,
    hi = CFG$s2_expo_hi,
    min_obs = CFG$min_expo_obs,
    min_span = CFG$min_span_fraction,
    n_grid = CFG$n_grid,
    method = CFG$smooth_method
  )
  
  log_line(
    "S2 step 1, with BMI exposure curve: ",
    length(exposure$ids)
  )
  
  if (!CFG$outcome2 %in% names(long)) {
    stop(
      "Configured Study 2 outcome is unavailable: ",
      CFG$outcome2
    )
  }
  
  outcome <- build_curves(
    long_data = long[
      hhidpn %in% exposure$ids
    ],
    variable = CFG$outcome2,
    lo = CFG$s2_out_lo,
    hi = CFG$s2_out_hi,
    min_obs = CFG$min_out_obs,
    min_span = CFG$min_span_fraction,
    n_grid = CFG$n_grid,
    method = CFG$smooth_method
  )
  
  log_line(
    "S2 step 2, with ",
    CFG$outcome2,
    " outcome curve: ",
    length(outcome$ids)
  )
  
  ids <- intersect(
    exposure$ids,
    outcome$ids
  )
  
  alive_ids <- baseline[
    hhidpn %in% ids &
      (
        !is.finite(age_at_death) |
          age_at_death > CFG$s2_out_hi
      ),
    hhidpn
  ]
  
  ids <- intersect(
    ids,
    alive_ids
  )
  
  log_line(
    "S2 step 3, alive through age ",
    CFG$s2_out_hi,
    ": ",
    length(ids)
  )
  
  dat <- baseline[
    hhidpn %in% ids
  ]
  
  confounders <- confounder_formula(
    CFG$conf_set
  )
  
  complete <- stats::complete.cases(
    dat[
      ,
      ..confounders
    ]
  )
  
  dat <- dat[complete]
  
  log_line(
    "S2 step 4, complete baseline confounders: ",
    nrow(dat)
  )
  
  if (nrow(dat) < 20L) {
    stop(
      "Study 2 has fewer than 20 complete subjects."
    )
  }
  
  key <- as.character(
    dat$hhidpn
  )
  
  X_raw <- exposure$X[
    key,
    ,
    drop = FALSE
  ]
  
  Y_raw <- outcome$X[
    key,
    ,
    drop = FALSE
  ]
  
  mean_x <- colMeans(X_raw)
  
  X_centered <- sweep(
    X_raw,
    2,
    mean_x,
    "-"
  )
  
  list(
    dat = dat,
    X = X_centered,
    X_raw = X_raw,
    Y = Y_raw,
    mean_curve = mean_x,
    grid_s = exposure$grid,
    grid_t = outcome$grid,
    outcome_name = CFG$outcome2,
    exposure_diagnostics =
      exposure$diagnostics,
    outcome_diagnostics =
      outcome$diagnostics
  )
}


study1 <- build_study1()
study2 <- build_study2()

saveRDS(
  list(
    study1 = study1,
    study2 = study2,
    long = long,
    baseline = baseline,
    CFG = CFG
  ),
  file.path(
    CFG$out_root,
    "cohorts.rds"
  )
)

log_line(
  "Study 1 n = ",
  nrow(study1$dat),
  "; events = ",
  sum(study1$dat$Y),
  " (",
  round(
    100 * mean(study1$dat$Y),
    1
  ),
  "%)"
)

log_line(
  "Study 2 n = ",
  nrow(study2$dat)
)
