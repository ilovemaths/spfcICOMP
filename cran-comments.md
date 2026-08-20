## Initial CRAN submission

This is the pre-submission record for the first `spfcICOMP` CRAN release. It will be updated with the exact v0.1.0 release-candidate results immediately before submission.

## Current development checks

For development version 0.0.0.9009:

* 149 tests pass.
* Local Windows CRAN-style check: 0 errors | 0 warnings | 0 notes.
* GitHub Actions R CMD check succeeds on macOS release, Windows release, Ubuntu release, Ubuntu R-devel, and Ubuntu oldrel-1.
* Exact R 3.5.0 compatibility audit passes.
* Validation stratification hardening is verified at commit `fc80e696b7a64a2b8ab7c274f0af6c13e0fbb8ff`.
* Both public-repository commits passed all five GitHub Actions matrix jobs.

The exact v0.1.0 source tarball will be rebuilt and rechecked before submission.

## External benchmark data

Third-party Riboflavin and Golub benchmark datasets are not included in the CRAN source package. Optional thesis-reproducibility scripts, excluded from the source tarball, obtain Riboflavin from the archived CRAN `hdi` 0.1-10 source and Golub from Bioconductor package `multtest` when requested.

## Reverse dependencies

This is an initial CRAN submission, so no CRAN reverse dependencies are expected before publication.

## Release-freeze policy

No statistical methodology will be changed during final CRAN packaging. Release-engineering, documentation, portability, and policy-compliance changes will be revalidated with the full test suite and R CMD checks.
