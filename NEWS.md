# spfcICOMP 0.0.0.9009

- Began GitHub and CRAN release hardening without changing the statistical methodology or thesis-result implementation.
- Added repository metadata, README, citation metadata, and CRAN-preparation notes.
- Excluded thesis-analysis artefacts and exploratory outputs from the CRAN source tarball while retaining them for GitHub reproducibility.
- Replaced over-strong public wording about sparsity/selection with feature-screening terminology where appropriate.
- Removed third-party Riboflavin and Golub benchmark copies from the distributable package; the GitHub thesis-analysis layer now obtains Riboflavin from the archived CRAN `hdi` 0.1-10 source (with an installed-`hdi` fast path when available) and Golub from Bioconductor `multtest` when empirical reproduction is requested.
- Passed a local CRAN-style check with 0 errors, 0 warnings, and 0 notes and established a five-platform GitHub Actions R CMD check workflow.
- Hardened validation stratification so numeric binary responses are treated as categorical strata in train-test splitting and cross-validation fold construction.
- Added regression tests for numeric binary-response stratification without changing the statistical methodology or thesis-result implementation.
- Verified 149/149 local tests, a local R CMD check with 0 errors, 0 warnings, and 0 notes, exact R 3.5.0 compatibility, and the five-platform GitHub Actions matrix at commit fc80e696b7a64a2b8ab7c274f0af6c13e0fbb8ff.

# spfcICOMP 0.0.0.9008

- Align the standalone validity-check subspace tolerance with the package unit-test tolerance (1e-7); no methodological or simulation logic changed.

* Relaxed the numerical tolerance in the DGP subspace-equivalence unit test from `1e-8` to `1e-7`; no simulation or estimation logic changed.

# spfcICOMP 0.0.0.9006

- Replaces the collapsed forward continuous simulation with a model-faithful inverse PFC DGP whose nominal d=2 inverse mean has rank two by construction.
- Constructs the inverse loading so `span(Sigma^{-1} Gamma)` equals the sparse true reduction subspace, preserving valid active-variable evaluation under correlated AR(1) residual covariance.
- Defines SNR in the true sufficient coordinates and scales the inverse mean so the requested SNR is exact in population and invariant to p.
- Generates an independent test sample for every continuous replicate and evaluates RMSE/MAE out of sample only.
- Uses a fixed degree-2 continuous response basis in thesis simulation runners, sufficient to represent both d=1 and d=2 without using the true d to tune the fitted basis.
- Adds regression tests guarding against rank collapse and accidental in-sample prediction evaluation.
- Explicitly restricts the current binary PFC simulation to d=1 rather than silently treating a rank-one binary response basis as a multi-dimensional recovery problem.

# spfcICOMP 0.0.0.9005

- Froze the thesis C1F-to-L1 calibration as `icomp_hd_floor`.
- The new default global row penalty is `median(row_norms) * sqrt(2*log(p)/n) * (1 + C1F/(1+C1F))`.
- This separates the high-dimensional multiplicity/noise floor from the bounded covariance-complexity multiplier without introducing a user-tuned LASSO lambda.
- Retained the 0.0.0.9004 `robust_universal` and historical `raw` calibrations for reproducibility.
- Added package tests for the exact frozen penalty formula and default calibration.
- Added high-dimensional-factor and complexity-multiplier diagnostics to simulation summaries.

# spfcICOMP 0.0.0.9004

- Recalibrated C1F-driven row-group sparsity with a scale-compatible, deterministic high-dimensional penalty.
- Removed the 1/p attenuation caused by sum-to-one adaptive-weight normalisation in the thesis selection route.
- Preserved the raw C1F route as a reproducibility option.
- Fixed simulation summaries so SNR and categorical signal-strength levels remain separate factorial cells.
- Added C1F penalty diagnostics to simulation output and smoke-test reporting.

# spfcICOMP 0.0.0.9001

* Added genuine structural-dimension recovery to simulation design runs. The
  true structural dimension is no longer supplied to the estimator in the
  thesis workflow; AIC, BIC, CAIC, ICOMP(IFIM), and CICOMP each select d from a
  candidate grid.
* Added C1F-driven adaptive row-L1 feature selection via
  `selection_rule = "c1f"`, while retaining legacy quantile thresholding for
  backwards compatibility.
* Added dimension recovery indicators, C1F penalty values, and confusion-count
  fields to simulation outputs.
* Updated simulation summaries to report dimension recovery rate and selected-d
  mean/SD by criterion.
* Added a thin `thesis_analysis/` reproducibility workflow built entirely on
  package functions.
