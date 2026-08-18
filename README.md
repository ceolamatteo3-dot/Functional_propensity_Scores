# Causal Inference with Functional Exposures: Evaluating Midlife BMI Trajectories on Metabolic and Cognitive Aging in RAND HRS

This repository implements a **Functional Causal Inference** pipeline using **Functional Propensity Score (FPS) weighting** and **Marginal Structural Models (MSMs)** to evaluate the causal impact of midlife Body Mass Index (BMI) trajectories on late-life health outcomes. Using longitudinal data from the **RAND Health and Retirement Study (HRS 1992–2022)**, this project investigates:
* **Study 1 (Scalar MSM)**: The effect of BMI trajectories ($t \in [52, 64]$) on the risk of late-onset incident Type 2 Diabetes ($t \in (64, 76]$).
* **Study 2 (Function-on-Function MSM)**: The causal effect surface relating midlife BMI trajectories ($s \in [52, 62]$) to subsequent cognitive performance trajectories (`cog27`, $t \in [64, 74]$).

```text
                     ┌────────────────────────────────────────────────────────┐
                     │              1. Longitudinal Smoothing                 │
                     │  - Reconstruct individual continuous trajectories      │
                     │    via natural splines / polynomials                   │
                     └──────────────────────────┬─────────────────────────────┘
                                                │
                                                ▼
                     ┌────────────────────────────────────────────────────────┐
                     │          2. Dimension Reduction via FPCA               │
                     │  - Extract dominant eigenfunctions ϕ_k(t) and scores A_k│
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
Estimating causal effects of continuous longitudinal trajectories requires addressing infinite-dimensional exposures while adjusting for baseline confounding. The Functional Propensity Score (FPS) solves the dual Lagrangian formulation of empirical likelihood balancing weights, eliminating linear dependence between treatment Functional Principal Component (FPC) scores and scalar baseline confounders (demographics, socioeconomic status, cardiovascular comorbidities).  Our design in Study 1 explicitly departs from previous concurrent setups (e.g., Ciardulli et al., 2026) to avoid reverse causation: by strictly separating the exposure window $[52, 64]$ from the subsequent follow-up window $(64, 76]$ and requiring participants to be diabetes-free through age 64, we eliminate post-diagnosis weight contamination and reverse causality. Household-cluster studentized bootstrap algorithms are employed to construct both pointwise and simultaneous uniform confidence bands across continuous domains.  In Study 1 (Incident Late-Onset Diabetes), the estimated effect function $\hat{\mu}(t)$ increases monotonically across the midlife exposure window. Excess body mass immediately prior to the follow-up window (ages 60–64) carries the highest marginal risk contribution ($\hat{\mu}(64) \approx +0.023$, odds multiplier $\text{OR} \approx 1.023$ per unit of centered BMI), achieving statistical significance under both pointwise and simultaneous uniform confidence bands. Early-midlife BMI ($t = 52\text{--}54$) displays a near-zero/slightly negative point estimate ($\hat{\mu}(52) \approx -0.011$), reflecting cohort conditioning on remaining disease-free through age 64: individuals with elevated BMI at age 52 who reach age 64 without diabetes constitute a metabolically resilient subpopulation, while recent late-midlife weight gain acts as the primary acute trigger for late-onset incidence. The integral over the full window is strictly positive ($\int_{52}^{64} \hat{\mu}(t)\,dt > 0$), confirming that a sustained, uniform elevation in BMI across midlife increases overall lifetime risk of incident diabetes.In Study 2 (Cognitive Trajectories cog27), the estimated causal surface $\hat{\mu}(s, t)$ reveals an exposure-timing gradient where higher BMI in early midlife ($s \approx 52\text{--}54$) is associated with lower cognitive scores across older ages, whereas higher BMI in late midlife ($s \approx 60\text{--}62$) exhibits a positive marginal association with cognitive performance in older age (peaking at outcome age $t = 74$). These findings align with literature on the "obesity paradox" in cognitive aging, where midlife adiposity poses long-term vascular/cognitive burdens, while modest body mass in older age provides metabolic/frailty reserve and reflects protection against preclinical neurodegenerative weight loss. Comparing naive unweighted regressions with FPS-weighted models demonstrates substantial confounding bias in unadjusted estimates, which is successfully mitigated through propensity score reweighting.  Plaintext├── R/
│   ├── 00_config.R             # Global parameters, grids, and seeds
│   ├── 01_build_cohort.R       # RAND HRS extraction, inclusion filters, curve fitting
│   ├── 02_eda.R                # FPCA diagnostics, perturbation plots, confounding heatmaps
│   ├── 03_analysis1_scalar.R   # Study 1: Binary MSM, FPS weighting, bootstrap bands
│   ├── 04_analysis2_surface.R  # Study 2: Function-on-function MSM, 2D effect surface & slices
│   ├── fps_functions.R         # Core numerical utilities, dense FPCA, and WLS solvers
│   └── fps_extra.R             # Survey-weighted base-measure extensions
├── outputs/                    # Generated figures, tables, and serialized objects
├── run_all.R                   # Master execution script
└── README.md
The pipeline requires R (≥ 4.1.0) and the packages haven, data.table, tidyselect, and splines. To run the analysis, download the RAND HRS Longitudinal File (1992–2022, v1) (randhrs1992_2022v1.dta) from the Health and Retirement Study portal, place it in the data/ folder, and execute source("run_all.R").  ReferencesCiardulli, S., Fontana, N., Vantini, S., & Ieva, F. (2026). Generalized Propensity Score Weighting for Functional Causal Inference Framework. arXiv:2608.03200.  Zhang, X., Xue, W., & Wang, Q. (2021). Covariate balancing functional propensity score for functional treatments in cross-sectional observational studies. Computational Statistics & Data Analysis, 163, 10730.  Ramsay, J. O., & Silverman, B. W. (2005). Functional Data Analysis. Springer[cite: 1].
