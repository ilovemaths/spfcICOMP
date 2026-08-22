## Initial CRAN submission

This is the pre-submission record for the first `spfcICOMP` CRAN release.

## Exact release-candidate checks

For release candidate version 0.1.0:

* 149 of 149 test expectations pass, with no failures, skips, warnings, or errors.
* The exact source tarball is `spfcICOMP_0.1.0.tar.gz`, containing 97 entries and measuring 77,678 bytes.
* Local Windows R 4.6.0 `R CMD check --as-cran --no-manual`: 0 errors | 0 warnings | 0 notes.
* The authoritative `00check.log` ends with `Status: OK`.
* All required source-package files are present and no forbidden repository, thesis-analysis, generated-output, reverse-dependency, or local-build artefacts are included.
* Standalone citation evaluation passes and the British-English spelling audit reports no flags.
* Exact R 3.5.0 compatibility was verified during development hardening.

The earlier public development commits and the exact 0.1.0 release-candidate commit `e48f48e26f3e607690a9809129b21ec7665bd7d1` passed the configured five-platform GitHub Actions matrix on macOS release, Windows release, Ubuntu release, Ubuntu R-devel, and Ubuntu oldrel-1. The release-candidate workflow run completed successfully on 21 August 2026. The checked source-package content has not changed since local tarball verification.

## External benchmark data

Third-party Riboflavin and Golub benchmark datasets are not included in the CRAN source package. Optional thesis-reproducibility scripts, excluded from the source tarball, obtain Riboflavin from the archived CRAN `hdi` 0.1-10 source and Golub from Bioconductor package `multtest` when requested.

## Reverse dependencies

This is an initial CRAN submission, so no CRAN reverse dependencies are expected before publication.

## Release-freeze policy

No statistical methodology will be changed during final CRAN packaging. Release-engineering, documentation, portability, and policy-compliance changes will be revalidated with the full test suite and R CMD checks.
