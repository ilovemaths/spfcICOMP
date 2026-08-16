# Thesis-method update in version 0.0.0.9001

This development snapshot makes the package simulation workflow consistent with
the intended SPFC-ICOMP thesis design.

## Structural dimension

`run_spfc_design_row()` no longer passes the true simulation dimension to the
SPFC estimator. It calls `spfc_select_dimension()` on a candidate grid and
reports one result row for each requested information criterion. Candidate SPFC
fits are reused across criteria.

Default thesis criteria:

- AIC
- BIC
- CAIC
- ICOMP(IFIM)
- CICOMP

The simulation output contains `fitted_d` and `dimension_correct`, allowing
Monte Carlo estimates of P(d_hat = d_true), mean selected d, and SD(selected d).

## Feature selection

`spfc_select_variables()` now supports:

- `selection_rule = "c1f"`: recommended thesis route. C1F from the fitted-model
  covariance is used as the data-adaptive scalar L1 term. Adaptive row weights
  distribute that total penalty across rows of the SPFC loading matrix. The
  row-group proximal update sets some row norms exactly to zero, so selected
  model cardinality is data-determined.
- `selection_rule = "quantile"`: legacy route retained for reproducibility.
- `selection_rule = "fixed"`: explicit user-supplied threshold.

The output includes raw importance, adaptive weight, row penalty, shrunk
importance, selected status, and the C1F value as an attribute.

## Simulation subspace metric

Subspace distance is calculated when the selected dimension equals the true
dimension. Incorrect dimension selection is handled by the separate
`dimension_correct` outcome rather than by silently truncating unequal
subspaces.

## Criterion definitions

The existing package definitions of AIC, BIC, CAIC, ICOMP(IFIM),
ICOMP(Misspec), and CICOMP were intentionally not changed. They use the
maximised reduced-model likelihood and its estimated parameter covariance,
consistent with the package's current model-covariance interpretation.


## 0.0.0.9002 categorical C1F fix

The ridge-logistic reduced model now stores its penalised observed Hessian and
its inverse covariance approximation. `extract_ifim()` uses this matrix for
C1F-based variable selection with binary categorical responses. Dimension
selection remains based on the ordinary logistic companion model, so the
existing likelihood-based information criteria are unchanged.


## 0.0.0.9005 frozen C1F calibration

Before the definitive Monte Carlo experiment, the thesis variable-selection rule was frozen as

`lambda_IC-HD = s_Gamma * sqrt(2 log(p) / n) * [1 + C1F/(1 + C1F)]`,

where `s_Gamma` is the median positive row norm of the fitted SPFC loading matrix. Adaptive row weights are rescaled to have median one before the row-group proximal update. This provides an explicit high-dimensional multiplicity/noise floor and a bounded C1F complexity multiplier, with no cross-validated or manually tuned L1 parameter. Version 0.0.0.9004's `robust_universal` rule and the historical `raw` route remain available only for reproducibility.

## 0.0.0.9006 simulation-validity correction

The first 0.0.0.9005 definitive simulation was audited before thesis reporting.
Two objective validity defects were identified: (i) `rowSums(X %*% B)` made the
nominal d=2 forward response algebraically equal to a single linear index, and
(ii) simulation RMSE/MAE were evaluated on the training observations.

Version 0.0.0.9006 replaces that generator with a model-faithful inverse PFC
DGP. For d=2, `f(Y)` contains the population-orthonormal functions Y and
`(Y^2 - 1)/sqrt(2)`. With residual covariance Sigma, sparse true reduction basis
B, and G = B' Sigma B, the inverse loading is `Gamma = Sigma B G^{-1/2}`.
Consequently `span(Sigma^{-1} Gamma) = span(B)`, so the sparse truth used for
subspace and variable-selection evaluation is mathematically aligned with the
PFC reduction target. The inverse mean is multiplied by `sqrt(SNR)`, giving the
requested SNR exactly in the sufficient coordinates and avoiding a p-dependent
signal inflation.

Each continuous replicate now includes an independent test sample of size n.
All fitting and selection are performed on training data only, and RMSE/MAE are
computed exclusively on the independent test sample. The thesis simulation
uses a fixed degree-2 response basis for both d=1 and d=2. The current binary
simulation is explicitly restricted to d=1 until a validated multiclass basis
is implemented.
