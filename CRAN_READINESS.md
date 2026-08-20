# spfcICOMP CRAN Readiness Record

## Baseline

- Source baseline: 0.0.0.9008.
- CRAN/GitHub hardening version: 0.0.0.9009.
- Statistical methodology and frozen thesis-result implementation: unchanged.
- Default branch: `main`.
- GitHub Actions: standard five-platform R CMD check workflow configured.

## Verified checks completed through 20 August 2026

The current v0.0.0.9009 development state has:

- 149 package tests passing;
- local `devtools::check(cran = TRUE, manual = TRUE)` with **0 errors, 0 warnings, 0 notes**;
- successful GitHub Actions checks on macOS release, Windows release, Ubuntu release, Ubuntu R-devel, and Ubuntu oldrel-1.
- exact R 3.5.0 compatibility audit passed;
- numeric binary-response stratification hardened and verified at commit `fc80e696b7a64a2b8ab7c274f0af6c13e0fbb8ff`;
- both the initial public snapshot and the subsequent hardening commit passed all five GitHub Actions matrix jobs.

These checks are development-hardening evidence. The exact future v0.1.0 source tarball must still be rebuilt and checked again immediately before CRAN submission.

## Resolved: benchmark-data redistribution separation

Third-party Riboflavin and Golub benchmark copies have been removed from the distributable package.

The excluded `thesis_analysis/` research layer obtains the benchmarks only when the empirical analyses are deliberately reproduced:

- Riboflavin from the archived CRAN `hdi` 0.1-10 source via an archive-safe loader;
- Golub leukaemia from Bioconductor package `multtest`.

Neither external data package is required for package installation, examples, vignettes, unit tests, or CRAN checks. This keeps third-party benchmark data outside the MIT-licensed `spfcICOMP` distribution.

## Remaining CRAN release audit

### 1. Public API review

The current exports were reviewed as stable user-facing, advanced research, and plotting interfaces. No export was identified as clearly accidental. The existing API is therefore retained for the initial 0.1.0 release to avoid unnecessary late-stage code and documentation churn.

### 2. Examples and documentation

Add concise, fast executable examples to the principal user-facing functions. Examples must not run heavy Monte Carlo simulations, require optional benchmark data, or write outside temporary directories.

### 3. Dependency, URL, spelling, and portability checks

Before release, run `urlchecker::url_check()`, run a package spelling check, review declared imports/suggests, confirm source-package and installed-package sizes, inspect licence and citation metadata, and verify no generated research output enters the source tarball.

### 4. Exact release-candidate checks

Only after the audit is complete: bump 0.0.0.9009 to 0.1.0; build the exact source tarball; check that exact tarball with `R CMD check --as-cran`; rerun the GitHub Actions matrix on the release commit; update `cran-comments.md`; and submit the exact checked source tarball.

After CRAN acceptance, tag the accepted commit `v0.1.0`, create a GitHub release, and advance the development version.
