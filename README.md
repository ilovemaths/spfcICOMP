# spfcICOMP

<!-- badges: start -->
[![R-CMD-check](https://github.com/ilovemaths/spfcICOMP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ilovemaths/spfcICOMP/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/spfcICOMP)](https://CRAN.R-project.org/package=spfcICOMP)
[![DOI](https://zenodo.org/badge/1335997201.svg)](https://doi.org/10.5281/zenodo.22310094)
<!-- badges: end -->

`spfcICOMP` implements Shrinkage Principal Fitted Components (SPFC) for
high-dimensional sufficient dimension reduction, structural-dimension
selection, feature screening, regression, and classification.

The package currently provides:

- regularised covariance estimation, including Oracle Approximating Shrinkage
  (OAS) and Maximum Entropy Covariance (MEC);
- Principal Fitted Components estimation for continuous and categorical
  responses;
- structural-dimension selection using AIC, BIC, CAIC, ICOMP(IFIM), and
  CICOMP;
- C1F-calibrated original-feature screening;
- reduced-space regression and classification;
- simulation and benchmarking utilities.

## Release status

Version 0.1.0 is published on CRAN and archived on Zenodo. The release
passed 149 test expectations, the five-platform GitHub Actions matrix,
the declared R 3.5.0 compatibility audit, and CRAN incoming pretests
with no errors or warnings.

- CRAN: https://CRAN.R-project.org/package=spfcICOMP
- GitHub release: https://github.com/ilovemaths/spfcICOMP/releases/tag/v0.1.0
- Version-specific Zenodo DOI: https://doi.org/10.5281/zenodo.22310095
- Concept DOI for all versions: https://doi.org/10.5281/zenodo.22310094

## Installation

Install the released version from CRAN:

```r
install.packages("spfcICOMP")
library(spfcICOMP)
```

## Development version

The development version can be installed from GitHub with either `pak` or
`remotes`:

```r
pak::pak("ilovemaths/spfcICOMP")
# or
remotes::install_github("ilovemaths/spfcICOMP")
```

## Basic example

```r
library(spfcICOMP)

set.seed(123)
X <- matrix(rnorm(40 * 10), nrow = 40, ncol = 10)
y <- X[, 1] - 0.5 * X[, 2] + rnorm(40)

fit <- spfc_fit(
  X = X,
  y = y,
  d = 1,
  ytype = "continuous",
  cov_method = "mec",
  nslices = 5,
  poly_degree = 2
)

fit
summary(fit)
head(fitted(fit))
```

Automatic structural-dimension selection is available through
`spfc_select_dimension()`:

```r
dsel <- spfc_select_dimension(
  X = X,
  y = y,
  d_grid = 1:3,
  cov_method = "mec",
  ytype = "continuous"
)

dsel$selected
```

## Benchmark data and thesis reproducibility

Third-party Riboflavin and Golub gene-expression datasets are not bundled in
the CRAN package. The research scripts under `thesis_analysis/` obtain the
benchmarks from their established statistical-data packages when required:

- Riboflavin from the archived CRAN `hdi` 0.1-10 source via an archive-safe loader;
- Golub leukaemia from Bioconductor package `multtest`.

These optional research-data packages are not required for the normal
installation, tests, vignettes, or examples of `spfcICOMP`. See
`thesis_analysis/README.md` for the frozen empirical-analysis settings and run
order.

## Citation

To obtain the package citation in R, run:

```r
citation("spfcICOMP")
```

For the exact 0.1.0 software archive, use
https://doi.org/10.5281/zenodo.22310095. The stable concept DOI for the
software project is https://doi.org/10.5281/zenodo.22310094.

## Methodological references

- Cook, R. D. and Forzani, L. (2008). Principal Fitted Components for
  dimension reduction in regression. *Statistical Science*, 23(4), 485-501.
  doi:10.1214/08-STS275.
- Chen, Y., Wiesel, A., Eldar, Y. C. and Hero, A. O. (2010). Shrinkage
  algorithms for MMSE covariance estimation. *IEEE Transactions on Signal
  Processing*, 58(10), 5016-5029. doi:10.1109/TSP.2010.2053029.
- Bozdogan, H. (2000). Akaike's Information Criterion and recent developments
  in information complexity. *Journal of Mathematical Psychology*, 44(1),
  62-91. doi:10.1006/jmps.1999.1277.
- Olorede, K. O. and Yahya, W. B. (2019). A new covariance estimator for
  sufficient dimension reduction in high-dimensional and undersized sample
  problems. *arXiv*. doi:10.48550/arXiv.1909.13017.

## Licence

MIT.
