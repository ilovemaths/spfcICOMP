# Install the local spfcICOMP source package and run its tests.
# Run this script from the package root.

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("testthat", quietly = TRUE)) {
  install.packages("testthat", repos = "https://cloud.r-project.org")
}

remotes::install_local(".", upgrade = "never", dependencies = TRUE)
library(spfcICOMP)

testthat::test_dir("tests/testthat", reporter = "summary")
packageVersion("spfcICOMP")
