library(tinytest)

# Locate analysis/scripts/sim_study.R relative to the package root
# regardless of whether tests are run via tinytest::test_package()
# (installed package, working directory is the test directory) or
# pkgload::load_all() + run_test_dir() (working directory is the
# project root). Real assertions here replace the previous
# placeholder `expect_true(TRUE)` (whitepaper section 5, framing
# option (c), "the test suite is a placeholder").
#
# analysis/ is intentionally .Rbuildignore'd (it is report-driver
# code, not package API), so it is absent from the tarball R CMD
# check tests against; sim_study.R is unreachable from that sandbox
# by any relative path or here::here(). Skip the file rather than
# error in that case -- the CI "Run tests" step separately runs
# tinytest::run_test_dir() against the full source checkout, where
# analysis/ exists and these assertions do execute.
find_sim_script <- function() {
  candidates <- c(
    "analysis/scripts/sim_study.R",
    "../../analysis/scripts/sim_study.R",
    "../../../analysis/scripts/sim_study.R"
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) {
    return(normalizePath(hit[1]))
  }
  if (requireNamespace("here", quietly = TRUE)) {
    cand <- here::here("analysis", "scripts", "sim_study.R")
    if (file.exists(cand)) {
      return(normalizePath(cand))
    }
  }
  NULL
}

sim_script <- find_sim_script()
if (is.null(sim_script)) {
  exit_file(paste(
    "sim_study.R not reachable from working directory", getwd(),
    "-- analysis/ is not shipped in the built package (R CMD check",
    "sandbox); skipping"
  ))
}
source(sim_script)

set.seed(20260820)

## generate_data() ------------------------------------------------

dat <- generate_data(
  n_per_arm = 40, p_visits = 4, miss_rate = 0,
  miss_mech = "MCAR", delta = 0.25
)

expect_true(
  all(c("id", "time", "trt", "y", "baseline") %in% names(dat)),
  info = "generate_data() returns the expected columns"
)
expect_equal(
  nrow(dat), 80 * 5,
  info = "generate_data() produces n_total * (p_visits + 1) rows"
)
expect_equal(
  sum(is.na(dat$y)), 0L,
  info = "no missingness generated when miss_rate = 0"
)

dat_miss <- generate_data(
  n_per_arm = 40, p_visits = 4, miss_rate = 0.30,
  miss_mech = "MCAR", delta = 0.25
)
expect_true(
  sum(is.na(dat_miss$y)) > 0,
  info = "MCAR missingness with miss_rate = 0.30 produces some NAs"
)
expect_equal(
  sum(is.na(dat_miss$y[dat_miss$baseline])), 0L,
  info = "baseline (time = 0) values are never missing"
)

## truth_for_method() ---------------------------------------------

expect_equal(
  truth_for_method("SMA-slope", delta = 0.25, p_visits = 4),
  0.25,
  info = "SMA-slope estimand is delta (whitepaper M2)"
)
expect_equal(
  truth_for_method("MMRM", delta = 0.25, p_visits = 4),
  0.25,
  info = "MMRM (fit with continuous time) estimand is delta"
)
expect_equal(
  truth_for_method("SMA-change", delta = 0.25, p_visits = 4),
  0.25 * (4 + 1) / 2,
  info = paste(
    "SMA-change estimand is delta * (p_visits + 1) / 2",
    "(whitepaper M2)"
  )
)
expect_equal(
  truth_for_method("SMA-change", delta = 0.25, p_visits = 8),
  0.25 * (8 + 1) / 2,
  info = "SMA-change estimand scales with p_visits"
)

## fit_sma_slope() and fit_sma_change() recover their own truth ---

dat_large <- generate_data(
  n_per_arm = 400, p_visits = 4, miss_rate = 0,
  miss_mech = "MCAR", delta = 0.25
)

slope_fit <- fit_sma_slope(dat_large)
expect_true(
  is.finite(slope_fit$est),
  info = "fit_sma_slope() returns a finite estimate"
)
expect_true(
  abs(slope_fit$est - 0.25) < 0.05,
  info = paste(
    "fit_sma_slope() recovers delta = 0.25 within Monte Carlo",
    "tolerance at n_per_arm = 400 (single large-sample check,",
    "not a full simulation)"
  )
)

change_fit <- fit_sma_change(dat_large)
change_truth <- truth_for_method(
  "SMA-change", delta = 0.25, p_visits = 4
)
expect_true(
  is.finite(change_fit$est),
  info = "fit_sma_change() returns a finite estimate"
)
expect_true(
  abs(change_fit$est - change_truth) < 0.15,
  info = paste(
    "fit_sma_change() recovers its own estimand,",
    "delta * (p_visits + 1) / 2, not a raw delta = 0.25",
    "(regression check for whitepaper M2)"
  )
)

## fit_mmrm() converges and targets the slope estimand ------------

if (requireNamespace("mmrm", quietly = TRUE)) {
  mmrm_fit <- fit_mmrm(dat_large)
  expect_true(
    isTRUE(mmrm_fit$converged),
    info = paste(
      "fit_mmrm() converges on a clean complete-data set;",
      "regression check for the mmrm::us() formula bug",
      "(whitepaper M1)"
    )
  )
  expect_true(
    is.finite(mmrm_fit$est),
    info = "fit_mmrm() returns a finite estimate when converged"
  )
  expect_true(
    abs(mmrm_fit$est - 0.25) < 0.1,
    info = paste(
      "fit_mmrm() (continuous time trt:time coefficient) recovers",
      "delta = 0.25 directly, not a visit-contrast multiple of it"
    )
  )
  expect_true(
    isTRUE(mmrm_fit$pval >= 0 && mmrm_fit$pval <= 1),
    info = paste(
      "fit_mmrm() pval lies in [0, 1]; regression check that the",
      "Pr(>|t|) column (5) is read from mmrm::summary()$coefficients,",
      "not the t value column (4), which is unbounded and produced",
      "~50% Type I error at the null in the production-scale run",
      "(remediation 2026-08-20, second pass)"
    )
  )
  expect_true(
    mmrm_fit$pval < 0.01,
    info = paste(
      "fit_mmrm() pval is small for a strong true effect",
      "(delta = 0.25, n_per_arm = 400); a t value miscoded as a",
      "p-value would not reliably satisfy this"
    )
  )
} else {
  cat("mmrm not installed; skipping fit_mmrm() tests\n")
}

## check_simulation_quality() catches broken pipelines -------------

good_results <- data.frame(
  method = c("SMA-slope", "MMRM"),
  convergence = c(1, 1),
  n_valid = c(100, 100),
  delta = c(0.25, 0.25),
  miss_mech = c("MCAR", "MCAR"),
  coverage = c(0.95, 0.95),
  mcse_coverage = c(0.01, 0.01)
)
expect_true(
  isTRUE(check_simulation_quality(good_results)),
  info = "check_simulation_quality() passes a well-formed result set"
)

zero_conv_results <- good_results
zero_conv_results$convergence[2] <- 0
zero_conv_results$n_valid[2] <- 0
expect_error(
  check_simulation_quality(zero_conv_results),
  info = paste(
    "check_simulation_quality() stops on zero-convergence rows,",
    "the exact failure mode of the mmrm::us() bug (whitepaper M1,",
    "m9)"
  )
)
