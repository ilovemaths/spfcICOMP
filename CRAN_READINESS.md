# spfcICOMP CRAN Readiness Record

## Baseline

- Source baseline: 0.0.0.9008.
- CRAN/GitHub hardening version: 0.0.0.9009.
- Initial CRAN release candidate: 0.1.0.
- Statistical methodology and frozen thesis-result implementation: unchanged.
- Default branch: `main`.
- GitHub Actions: standard five-platform R CMD check workflow configured.

## Verified checks completed through 21 August 2026

The v0.0.0.9009 development-hardening state has:

- 149 package tests passing;
- local `devtools::check(cran = TRUE, manual = TRUE)` with **0 errors, 0 warnings, 0 notes**;
- successful GitHub Actions checks on macOS release, Windows release, Ubuntu release, Ubuntu R-devel, and Ubuntu oldrel-1;
- exact R 3.5.0 compatibility audit passed;
- numeric binary-response stratification hardened and verified at commit `fc80e696b7a64a2b8ab7c274f0af6c13e0fbb8ff`;
- both public development commits passed all five GitHub Actions matrix jobs.

The exact v0.1.0 local release candidate has:

- 149 of 149 test expectations passing, with no failures, skips, warnings, or errors;
- a source tarball named `spfcICOMP_0.1.0.tar.gz`, containing 97 entries and measuring 77,678 bytes;
- all required package files present and no forbidden repository, thesis-analysis, generated-output, reverse-dependency, or local-build artefacts;
- a Windows R 4.6.0 `R CMD check --as-cran --no-manual` result of **0 errors, 0 warnings, 0 notes**, with authoritative `Status: OK`;
- standalone citation evaluation passing and no spelling flags under the declared British-English language metadata.

The exact 0.1.0 release-candidate commit `e48f48e26f3e607690a9809129b21ec7665bd7d1` was pushed to `main` and passed the configured five-platform GitHub Actions matrix on 21 August 2026. All local and remote release-candidate verification gates are complete.

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

The principal user-facing functions contain concise, fast examples that do not require optional benchmark data or write outside temporary directories. The examples, donttest examples, tests, and vignette rebuild completed successfully in the pre-bump source-tarball check.

### 3. Dependency, URL, spelling, and portability checks

The dependency declaration was reviewed and the only strong dependencies are R and MASS. Manual extended-timeout checks returned HTTP 200 for the repository and issue-tracker URLs. British English is declared and the reviewed spelling audit returns no flags. The pre-bump source tarball contained 96 entries, was 77,511 bytes, contained all required package files, and contained no excluded thesis, output, GitHub, reverse-dependency, or local-build artefacts. Licence and standalone citation metadata were also verified.

### 4. Exact release-candidate checks

The exact local 0.1.0 source tarball passed the full release-candidate test, content, citation, spelling, and `R CMD check --as-cran` gates. Commit `e48f48e26f3e607690a9809129b21ec7665bd7d1` also passed the configured five-platform GitHub Actions matrix. The checked package content has not changed since these verifications. No release tag has been created and CRAN submission has not started.

After CRAN acceptance, tag the accepted commit `v0.1.0`, create a GitHub release, and advance the development version.
