# Causal Inference with Functional Exposures: Evaluating Midlife BMI Trajectories on Metabolic and Cognitive Aging in RAND HRS

This repository implements a **Functional Causal Inference** pipeline using **Functional Propensity Score (FPS) weighting** and **Marginal Structural Models (MSMs)** to evaluate the causal impact of midlife Body Mass Index (BMI) trajectories on late-life health outcomes.

Using longitudinal data from the **RAND Health and Retirement Study (HRS 1992–2022)**, this project investigates:

- **Study 1 — Scalar MSM**:  
  The effect of midlife BMI trajectories ($t \in [52, 64]$) on the risk of late-onset incident Type 2 Diabetes during follow-up ($t \in (64, 76]$).

- **Study 2 — Function-on-Function MSM**:  
  The causal effect surface relating midlife BMI trajectories ($s \in [52, 62]$) to subsequent cognitive performance trajectories (`cog27`, $t \in [64, 74]$).

---

## Overview

Estimating causal effects of continuous longitudinal trajectories requires addressing infinite-dimensional exposures while adjusting for baseline confounding.

The **Functional Propensity Score (FPS)** solves the dual Lagrangian formulation of empirical likelihood balancing weights, eliminating linear dependence between treatment Functional Principal Component (FPC) scores and scalar baseline confounders, including demographics, socioeconomic status, and cardiovascular comorbidities.

A key design feature of **Study 1** is the explicit separation between the exposure window and the outcome window to avoid reverse causation. By requiring participants to be diabetes-free through age 64 and by modeling BMI exposure only before the follow-up period, the analysis avoids post-diagnosis weight contamination and reverse causality.

Household-cluster studentized bootstrap algorithms are used to construct both **pointwise** and **simultaneous uniform confidence bands** across continuous exposure or outcome domains.

---

## Analytical Pipeline

```text
┌────────────────────────────────────────────────────────┐
│          1. Longitudinal Smoothing                     │
│  - Reconstruct individual continuous trajectories      │
│    via natural splines / polynomials                   │
└──────────────────────────┬─────────────────────────────┘
                           │
                           ▼
┌────────────────────────────────────────────────────────┐
│          2. Dimension Reduction via FPCA               │
│  - Extract dominant eigenfunctions ϕ_k(t) and scores  │
│    A_k                                                 │
│  - Retain components up to fixed PVE threshold         │
└──────────────────────────┬─────────────────────────────┘
                           │
                           ▼
┌────────────────────────────────────────────────────────┐
│          3. Functional Propensity Scores (FPS)         │
│  - Dual unconstrained convex optimization              │
│  - Exact balance: Cor_w(A_k, C_j) ≈ 0                  │
└──────────────────────────┬─────────────────────────────┘
                           │
                           ▼
┌────────────────────────────────────────────────────────┐
│          4. Marginal Structural Outcome Models         │
│  - Study 1: Weighted Logistic MSM                      │
│  - Study 2: Weighted Function-on-Function WLS Surface  │
│  - Cluster-bootstrap pointwise & simultaneous bands    │
└────────────────────────────────────────────────────────┘
```

---

## Study 1: Incident Late-Onset Diabetes

### Outcome

Incident late-onset Type 2 Diabetes during the follow-up window:

$$
t \in (64, 76]
$$

### Exposure

Midlife BMI trajectories:

$$
t \in [52, 64]
$$

### Design

This study explicitly departs from concurrent exposure-outcome setups to avoid reverse causation.

Key design choices include:

- Strict separation of the exposure window $[52, 64]$ from the follow-up window $(64, 76]$.
- Participants must remain diabetes-free through age 64.
- Post-diagnosis weight changes are excluded from the exposure construction.
- Baseline confounders are balanced using FPS weighting.

### Main Results

The estimated effect function $\hat{\mu}(t)$ increases monotonically across the midlife exposure window.

Excess body mass immediately prior to the follow-up window, approximately ages 60–64, carries the highest marginal risk contribution:

$$
\hat{\mu}(64) \approx +0.023
$$

corresponding to an odds multiplier of approximately:

$$
\text{OR} \approx 1.023
$$

per unit increase in centered BMI.

This effect is statistically significant under both pointwise and simultaneous uniform confidence bands.

Early-midlife BMI, around ages 52–54, shows a near-zero or slightly negative point estimate:

$$
\hat{\mu}(52) \approx -0.011
$$

This pattern reflects conditioning on remaining disease-free through age 64. Individuals with elevated BMI at age 52 who reach age 64 without diabetes may represent a metabolically resilient subpopulation, whereas recent late-midlife weight gain appears to act as the primary acute trigger for late-onset diabetes incidence.

The integral over the full exposure window is strictly positive:

$$
\int_{52}^{64} \hat{\mu}(t)\,dt > 0
$$

indicating that sustained, uniform elevation in BMI across midlife increases overall lifetime risk of incident diabetes.

---

## Study 2: Cognitive Trajectories

### Outcome

Late-life cognitive performance trajectories:

$$
t \in [64, 74]
$$

Outcome measure: `cog27`.

### Exposure

Midlife BMI trajectories:

$$
s \in [52, 62]
$$

### Main Results

The estimated causal surface $\hat{\mu}(s, t)$ reveals an exposure-timing gradient.

Higher BMI in early midlife:

$$
s \approx 52\text{–}54
$$

is associated with lower cognitive scores across older ages.

In contrast, higher BMI in late midlife:

$$
s \approx 60\text{–}62
$$

exhibits a positive marginal association with cognitive performance in older age, peaking around outcome age:

$$
t = 74
$$

These findings are consistent with the literature on the “obesity paradox” in cognitive aging:

- Midlife adiposity may impose long-term vascular and cognitive burdens.
- Modest body mass in older age may provide metabolic or frailty reserve.
- Higher late-life BMI may partly reflect protection against preclinical neurodegenerative weight loss.

Comparisons between naive unweighted regressions and FPS-weighted models show substantial confounding bias in unadjusted estimates. FPS propensity score reweighting substantially mitigates this bias.

---

## Methodological Components

### 1. Longitudinal Smoothing

Individual BMI trajectories are reconstructed from discrete longitudinal observations using spline or polynomial smoothing.

### 2. Functional Principal Component Analysis

FPCA is used to reduce infinite-dimensional functional exposures to a finite set of dominant functional principal component scores:

$$
A_k
$$

corresponding to eigenfunctions:

$$
\phi_k(t)
$$

Components are retained according to a fixed proportion of variance explained threshold.

### 3. Functional Propensity Score Weighting

FPS weights are obtained through a dual unconstrained convex optimization problem.

The balancing target is:

$$
\text{Cor}_w(A_k, C_j) \approx 0
$$

where:

- $A_k$ are treatment FPC scores.
- $C_j$ are scalar baseline confounders.
- $w$ denotes the estimated propensity score weights.

### 4. Marginal Structural Models

Two MSM frameworks are used:

- **Study 1**: Weighted logistic MSM for binary incident diabetes.
- **Study 2**: Weighted function-on-function weighted least squares surface model for cognitive trajectories.

### 5. Bootstrap Inference

Household-cluster studentized bootstrap methods are used to construct:

- Pointwise confidence bands.
- Simultaneous uniform confidence bands.

---

## Repository Structure

```text
├── R/
│   ├── 00_config.R             # Global parameters, grids, and seeds
│   ├── 01_build_cohort.R       # RAND HRS extraction, inclusion filters, curve fitting
│   ├── 02_eda.R                # FPCA diagnostics, perturbation plots, confounding heatmaps
│   ├── 03_analysis1_scalar.R   # Study 1: Binary MSM, FPS weighting, bootstrap bands
│   ├── 04_analysis2_surface.R  # Study 2: Function-on-function MSM, 2D effect surface & slices
│   ├── fps_functions.R         # Core numerical utilities, dense FPCA, and WLS solvers
│   └── fps_extra.R             # Survey-weighted base-measure extensions
│
├── data/                       # Place raw RAND HRS .dta file here
├── outputs/                    # Generated figures, tables, and serialized objects
├── run_all.R                   # Master execution script
└── README.md
```

---

## Requirements

The pipeline requires:

- **R** ≥ 4.1.0

Required R packages:

```r
install.packages(c(
  "haven",
  "data.table",
  "tidyselect",
  "splines"
))
```

Note: `splines` is included with base R, but is listed here because it is used in the pipeline.

---

## Data

This analysis uses the **RAND HRS Longitudinal File (1992–2022, v1)**.

Required file:

```text
randhrs1992_2022v1.dta
```

Download the file from the Health and Retirement Study data portal, then place it in the `data/` folder:

```text
data/randhrs1992_2022v1.dta
```

---

## Usage

To run the full pipeline, open R and execute:

```r
source("run_all.R")
```

This will run:

1. Cohort construction and longitudinal curve fitting.
2. Exploratory functional data analysis and diagnostics.
3. FPS estimation and confounder balancing.
4. Study 1 scalar MSM analysis.
5. Study 2 function-on-function MSM analysis.
6. Bootstrap inference and output generation.

Generated figures, tables, and serialized R objects will be written to:

```text
outputs/
```

---

## References

1. Ciardulli, S., Fontana, N., Vantini, S., & Ieva, F. (2026).  
   *Generalized Propensity Score Weighting for Functional Causal Inference Framework.*  
   arXiv:2608.03200.

2. Zhang, X., Xue, W., & Wang, Q. (2021).  
   *Covariate balancing functional propensity score for functional treatments in cross-sectional observational studies.*  
   Computational Statistics & Data Analysis, 163, 107300.

3. Ramsay, J. O., & Silverman, B. W. (2005).  
   *Functional Data Analysis.*  
   Springer.
