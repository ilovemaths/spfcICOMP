# Simulation validity correction for v0.0.0.9006

The v0.0.0.9005 continuous generator used `signal <- rowSums(X %*% B)`. For d=2 this equals `X %*% (b1 + b2)`, so the nominal two-dimensional forward model collapsed to one linear index. In addition, RMSE and MAE were evaluated on the same observations used to fit the reduced model.

Version 0.0.0.9006 corrects both issues before the replacement definitive simulation is run.

## Model-faithful inverse PFC DGP

The corrected generator follows the PFC inverse model

    X | Y = Gamma A f(Y) + epsilon,
    epsilon ~ N_p(0, Sigma),

with AR(1) residual covariance Sigma. The true sparse reduction basis B has its non-zero rows among the first s predictors. Let G = B' Sigma B. The inverse-mean loading matrix is

    Gamma = Sigma B G^{-1/2}.

Therefore

    span(Sigma^{-1} Gamma) = span(B),

so the true sufficient reduction subspace evaluated in simulation is exactly the sparse B supplied by the generator, even under correlated residual noise.

For d=1,

    f(Y) = Y.

For d=2,

    f(Y) = [Y, (Y^2 - 1)/sqrt(2)]',   Y ~ N(0,1).

These basis functions are population-centred, uncorrelated, and have unit variance. The inverse mean is scaled by sqrt(SNR), so in the true sufficient coordinates B'X,

    Var{B' E(X|Y)} = SNR * G,
    Var(B' epsilon) = G.

Thus the requested SNR is preserved exactly in population and does not increase artificially with p.

## Honest prediction

Each replicate generates an independent test sample of size n from the same population, using the same Sigma, B, Gamma and SNR. Covariance estimation, dimension selection, C1F feature selection and downstream-model fitting use training data only. RMSE and MAE are computed only on the independent test sample.

## Response basis

The thesis simulation uses the fixed degree-2 basis f(Y) = (Y, Y^2)' after centring through `build_fy(..., poly_degree = 2)` for both nominal d=1 and d=2 scenarios. This avoids using the true d to choose the fitted basis and ensures that the candidate model has sufficient response-basis rank to represent d=2.
