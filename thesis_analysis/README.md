# spfcICOMP thesis-analysis reproducibility layer

This directory contains the research-analysis scripts used to audit and report the
SPFC-ICOMP thesis results. It is intentionally excluded from the CRAN source
tarball by `.Rbuildignore`; the package in `R/` is the methodological
implementation.

## Frozen simulation provenance

The definitive simulation design was completed before the CRAN hardening pass.
Raw completed simulation results remain associated with the v0.0.0.9006
simulation archive, while the definitive audit and aggregation logic is the
v0.0.0.9008 reporting layer. Version 0.0.0.9009 changes release engineering and
documentation only; it does not redefine the statistical methodology.

## External benchmark data

The Riboflavin and Golub benchmark datasets are deliberately **not bundled in
the CRAN package**.

The empirical scripts obtain them at analysis time from established statistical
data packages:

- Riboflavin: the `riboflavin` object from archived CRAN package `hdi` 0.1-10. The analysis script uses an installed `hdi` copy if available and otherwise retrieves the frozen CRAN archive directly.
- Golub leukaemia: `multtest::golub` and `multtest::golub.cl` from the
  Bioconductor package `multtest`.

This keeps third-party benchmark data separate from the MIT-licensed
`spfcICOMP` source package while preserving a reproducible acquisition route.

Install the optional research-data dependencies when needed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("multtest", update = FALSE, ask = FALSE)
```

No installation of `hdi` is required for current R versions. The Riboflavin
analysis script automatically retrieves the frozen `hdi` 0.1-10 source archive
from CRAN when an installed copy is unavailable.

Neither the archived Riboflavin source nor `multtest` is required to install,
load, test, or check the CRAN package itself.

## Definitive real-data settings

Riboflavin:

- continuous response;
- OAS and MEC covariance estimators;
- candidate structural dimensions `1:5`;
- five criteria: AIC, BIC, CAIC, ICOMP(IFIM), CICOMP;
- response basis degree 2;
- five response slices;
- repeated 5-fold cross-validation with 5 repeats;
- seed `20260810`.

Golub leukaemia:

- binary response;
- fixed structural dimension `d = 1`;
- OAS and MEC covariance estimators;
- repeated stratified 5-fold cross-validation with 20 repeats;
- seed `20260810`.

## Run order

Core package and simulation audit:

1. `source("thesis_analysis/00_install_and_test.R")`
2. `source("thesis_analysis/01a_validity_check.R")`
3. `source("thesis_analysis/01_smoke_test.R")`
4. `source("thesis_analysis/02_full_simulation.R")` only when deliberately
   rerunning the full simulation rather than auditing the frozen archive.
5. `source("thesis_analysis/03_audit_and_aggregate.R")`
6. `source("thesis_analysis/04_simulation_figures.R")`

Definitive empirical applications:

7. `source("thesis_analysis/05_riboflavin_analysis.R")`
8. `source("thesis_analysis/06_golub_leukemia_analysis.R")`
9. `source("thesis_analysis/07_real_data_audit_and_aggregate.R")`
10. `source("thesis_analysis/08_real_data_figures.R")`

Do not tune the methodology from the definitive outcome tables. Any future
methodological change should be versioned and validated separately.
