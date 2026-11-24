library(tinytest)

# Tests for the supplementary methods added in response to the
# referee report of 2026-11-17 (fit_sma_slope_wls, fit_mmrm_visit,
# Kenward-Roger option of fit_mmrm, missingness_summary, new truths).
# Locating sim_study.R follows test_basic.R.
cands <- c(
  "analysis/scripts/sim_study.R",
  "../../analysis/scripts/sim_study.R",
  "../../../analysis/scripts/sim_study.R"
)
hit <- cands[file.exists(cands)]
if (length(hit) == 0) {
  exit_file("sim_study.R not reachable (R CMD check sandbox); skipping")
}
source(hit[1])

set.seed(20261117)

## truth_for_method() for the supplementary methods ---------------

expect_equal(
  truth_for_method("SMA-slope-W", delta = 0.25, p_visits = 8),
  0.25,
  info = "weighted slope summary targets delta"
)
expect_equal(
  truth_for_method("MMRM-KR", delta = 0.25, p_visits = 8),
  0.25,
  info = "continuous-time MMRM with KR df targets delta"
)
expect_equal(
  truth_for_method("MMRM-visit", delta = 0.25, p_visits = 8),
  0.25 * 8,
  info = "categorical-visit MMRM final-visit contrast targets delta * p"
)

## generate_data() default is unchanged by the mar_slope argument --

set.seed(1)
d_default <- generate_data(30, 4, 0.2, "MAR")
set.seed(1)
d_explicit <- generate_data(30, 4, 0.2, "MAR", mar_slope = 0.04)
expect_identical(
  d_default, d_explicit,
  info = "mar_slope = 0.04 reproduces the production generator"
)

## missingness_summary() --------------------------------------------

d_complete <- generate_data(20, 4, 0, "MCAR")
expect_equal(
  unname(missingness_summary(d_complete)), c(0, 0, 0),
  info = "no missingness reported for complete data"
)
d_mar <- generate_data(100, 8, 0.2, "MAR")
ms <- missingness_summary(d_mar)
expect_true(
  ms["final_missing"] >= ms["obs_missing"],
  info = "under monotone dropout the final visit is the most missing"
)

## fit_sma_slope_wls() ----------------------------------------------

d_big <- generate_data(300, 4, 0, "MCAR", delta = 0.25)
ols <- fit_sma_slope(d_big)
wls <- fit_sma_slope_wls(d_big)
expect_equal(
  wls$est, ols$est, tolerance = 1e-10,
  info = "with complete data all weights are equal, so WLS = OLS"
)
expect_true(
  abs(wls$est - 0.25) < 0.05,
  info = "weighted slope summary recovers delta in a large sample"
)

d_mar_big <- generate_data(300, 8, 0.2, "MAR", delta = 0.25)
wls_mar <- fit_sma_slope_wls(d_mar_big)
ols_mar <- fit_sma_slope(d_mar_big)
expect_true(
  is.finite(wls_mar$est) && is.finite(wls_mar$se),
  info = "weighted slope summary is finite under MAR dropout"
)
expect_true(
  wls_mar$se < ols_mar$se,
  info = paste(
    "under monotone dropout the weighted summary has a smaller",
    "model-based SE than the unweighted one (Major 3 mechanism)"
  )
)
expect_true(
  wls_mar$ci_lower < wls_mar$est && wls_mar$est < wls_mar$ci_upper,
  info = "weighted summary CI brackets its estimate"
)

## fit_mmrm(df_method = "Kenward-Roger") and fit_mmrm_visit() ------

if (requireNamespace("mmrm", quietly = TRUE)) {
  d_m <- generate_data(150, 4, 0, "MCAR", delta = 0.25)
  kr <- fit_mmrm(d_m, df_method = "Kenward-Roger")
  sat <- fit_mmrm(d_m)
  expect_true(
    isTRUE(kr$converged) && is.finite(kr$est),
    info = "fit_mmrm() converges with Kenward-Roger df"
  )
  expect_equal(
    kr$est, sat$est, tolerance = 1e-6,
    info = "df method changes inference, not the REML point estimate"
  )
  vis <- fit_mmrm_visit(d_m)
  expect_true(
    isTRUE(vis$converged) && is.finite(vis$est),
    info = "fit_mmrm_visit() converges on complete data"
  )
  # With complete data, a saturated-mean MMRM with unstructured
  # covariance reproduces the per-visit ANCOVA estimate exactly.
  wide <- merge(
    d_m[d_m$time == 4, c("id", "trt", "y")],
    setNames(d_m[d_m$time == 0, c("id", "y")], c("id", "y0")),
    by = "id"
  )
  ancova_p <- unname(coef(lm(y ~ trt + y0, data = wide))["trt"])
  expect_equal(
    vis$est, ancova_p, tolerance = 1e-4,
    info = paste(
      "fit_mmrm_visit() returns the final-visit treatment",
      "difference (equal to final-visit ANCOVA with complete data)"
    )
  )
  expect_true(
    vis$pval >= 0 && vis$pval <= 1,
    info = "fit_mmrm_visit() reads the Pr(>|t|) column"
  )
}
