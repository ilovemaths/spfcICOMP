SPFC-ICOMP REAL-DATA RECONCILIATION, EXTERNAL-DATA EDITION

Reconciled empirical workflow for the frozen v0.0.0.9008 thesis analysis, updated for external benchmark-data acquisition under v0.0.0.9009.

Files
-----
05_riboflavin_analysis.R
06_golub_leukemia_analysis.R
07_real_data_audit_and_aggregate.R
08_real_data_figures.R

Run order
---------
source("thesis_analysis/05_riboflavin_analysis.R")
source("thesis_analysis/06_golub_leukemia_analysis.R")
source("thesis_analysis/07_real_data_audit_and_aggregate.R")
source("thesis_analysis/08_real_data_figures.R")

Output directories
------------------
thesis_analysis/real_outputs_v9008/
thesis_analysis/real_tables_v9008/
thesis_analysis/real_figures_v9008/

Key design choices
------------------
Riboflavin:
- Continuous high-dimensional regression
- OAS and MEC
- d = 1,...,5
- AIC, BIC, CAIC, ICOMP_IFIM, CICOMP
- Polynomial response basis degree 2, matching the frozen simulation
- C1F ICOMP-HD feature screening
- Repeated 5-fold CV, 5 repeats
- Dimension selected within each training fold

Golub leukaemia:
- Binary high-dimensional classification
- OAS and MEC
- d = 1 fixed because the validated binary inverse-response basis has rank 1
- Information criteria reported at d = 1, not treated as competing d selectors
- Ridge-logistic reduced classifier
- C1F ICOMP-HD feature screening
- Repeated stratified 5-fold CV, 20 repeats
- Accuracy, balanced accuracy, sensitivity, specificity, precision, F1 and AUC


CRAN/GitHub release note:
Riboflavin and Golub data are no longer bundled in spfcICOMP. The v0.0.0.9009 research scripts load Riboflavin from the archived CRAN hdi 0.1-10 source (or an installed hdi copy when available) and Golub from Bioconductor package multtest. Numerical analysis settings and frozen v9008 output provenance are unchanged.
