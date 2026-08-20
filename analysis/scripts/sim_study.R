generate_data <- function(n_per_arm, p_visits, miss_rate,
                         miss_mech = "MCAR",
                         beta0 = 50, beta1 = -0.5,
                         delta = 0.25,
                         sigma_b0 = 5, sigma_b1 = 0.3,
                         rho_01 = -0.3, sigma_e = 3) {
  n_total <- 2 * n_per_arm
  trt <- rep(c(0, 1), each = n_per_arm)
  times <- seq(0, p_visits, by = 1)

  b0 <- rnorm(n_total, 0, sigma_b0)
  b1 <- rho_01 * (sigma_b1 / sigma_b0) * b0 +
    sqrt(1 - rho_01^2) * sigma_b1 * rnorm(n_total)

  dat <- expand.grid(id = 1:n_total, time = times) |>
    dplyr::mutate(
      trt = trt[id],
      b0_i = b0[id],
      b1_i = b1[id],
      y = beta0 + b0_i +
        (beta1 + b1_i) * time +
        delta * trt * time +
        rnorm(dplyr::n(), 0, sigma_e),
      baseline = (time == 0)
    )

  if (miss_rate > 0 && miss_mech == "MCAR") {
    post_bl <- which(dat$time > 0)
    drop <- rbinom(length(post_bl), 1, miss_rate)
    dat$y[post_bl[drop == 1]] <- NA
  }

  # MAR dropout is calibrated so that the marginal probability of
  # dropping out by the final visit is approximately `miss_rate`. The
  # constant per-visit hazard consistent with that marginal rate under
  # independent hazards would be `h_base` solved from
  # `1 - (1 - h_base)^p_visits = miss_rate`, but the informative
  # perturbation below (subjects with worse-than-average observed
  # outcomes are more likely to drop) is asymmetric under the logistic
  # link and systematically pulls the achieved marginal rate below
  # `h_base`'s nominal target. A multiplicative correction
  # (`mar_hazard_scale = 2`) was calibrated empirically (Monte Carlo
  # check over n_per_arm in {30, 50, 100}, p_visits in {4, 8},
  # miss_rate in {0.10, 0.20}, 60 replicates per cell) to bring the
  # achieved marginal rate within about 2-3 percentage points of the
  # nominal `miss_rate` (whitepaper m3); it is not an exact analytic
  # calibration and is disclosed as such in the manuscript Methods.
  if (miss_rate > 0 && miss_mech == "MAR") {
    mar_hazard_scale <- 2
    h_base <- 1 - (1 - miss_rate)^(1 / p_visits)
    intercept <- qlogis(h_base) + log(mar_hazard_scale)
    for (id_i in unique(dat$id)) {
      rows_i <- which(dat$id == id_i)
      for (j in 2:length(rows_i)) {
        if (is.na(dat$y[rows_i[j - 1]])) {
          dat$y[rows_i[j]] <- NA
          next
        }
        prev_y <- dat$y[rows_i[j - 1]]
        p_drop <- plogis(intercept + 0.04 * (beta0 - prev_y))
        p_drop <- min(p_drop, 1)
        if (runif(1) < p_drop) {
          dat$y[rows_i[j:length(rows_i)]] <- NA
          break
        }
      }
    }
  }

  dat
}

fit_sma_slope <- function(dat) {
  bl <- dat |>
    dplyr::filter(time == 0) |>
    dplyr::select(id, trt, baseline_y = y)

  slopes <- dat |>
    dplyr::filter(time > 0, !is.na(y)) |>
    dplyr::group_by(id) |>
    dplyr::filter(dplyr::n() >= 2) |>
    dplyr::summarise(
      slope = coef(lm(y ~ time))[2],
      .groups = "drop"
    )

  ana <- dplyr::inner_join(slopes, bl, by = "id")
  if (nrow(ana) < 6) {
    return(list(
      est = NA, se = NA, pval = NA,
      ci_lower = NA, ci_upper = NA, converged = TRUE
    ))
  }

  fit <- lm(slope ~ trt + baseline_y, data = ana)
  s <- summary(fit)$coefficients
  trt_row <- which(rownames(s) == "trt")

  list(
    est = s[trt_row, 1],
    se = s[trt_row, 2],
    pval = s[trt_row, 4],
    ci_lower = s[trt_row, 1] - qt(0.975, fit$df.residual) *
      s[trt_row, 2],
    ci_upper = s[trt_row, 1] + qt(0.975, fit$df.residual) *
      s[trt_row, 2],
    converged = TRUE
  )
}

# SMA-change: change-score ANCOVA. The Design specifications
# subsection of report.Rmd also describes an unadjusted two-sample
# t-test comparator; that variant is not implemented here (whitepaper
# m1). `fit_sma_change_ttest()` below adds it so text and code agree
# without changing the estimand handled by `fit_sma_change()`.
fit_sma_change <- function(dat) {
  bl <- dat |>
    dplyr::filter(time == 0) |>
    dplyr::select(id, trt, baseline_y = y)

  post_means <- dat |>
    dplyr::filter(time > 0, !is.na(y)) |>
    dplyr::group_by(id) |>
    dplyr::summarise(post_mean = mean(y), .groups = "drop")

  ana <- dplyr::inner_join(post_means, bl, by = "id") |>
    dplyr::mutate(change = post_mean - baseline_y)

  if (nrow(ana) < 6) {
    return(list(
      est = NA, se = NA, pval = NA,
      ci_lower = NA, ci_upper = NA, converged = TRUE
    ))
  }

  fit <- lm(change ~ trt + baseline_y, data = ana)
  s <- summary(fit)$coefficients
  trt_row <- which(rownames(s) == "trt")

  list(
    est = s[trt_row, 1],
    se = s[trt_row, 2],
    pval = s[trt_row, 4],
    ci_lower = s[trt_row, 1] - qt(0.975, fit$df.residual) *
      s[trt_row, 2],
    ci_upper = s[trt_row, 1] + qt(0.975, fit$df.residual) *
      s[trt_row, 2],
    converged = TRUE
  )
}

# Unadjusted two-sample t-test on the change score, as promised (but
# not previously implemented) by the Design specifications text.
fit_sma_change_ttest <- function(dat) {
  bl <- dat |>
    dplyr::filter(time == 0) |>
    dplyr::select(id, trt, baseline_y = y)

  post_means <- dat |>
    dplyr::filter(time > 0, !is.na(y)) |>
    dplyr::group_by(id) |>
    dplyr::summarise(post_mean = mean(y), .groups = "drop")

  ana <- dplyr::inner_join(post_means, bl, by = "id") |>
    dplyr::mutate(change = post_mean - baseline_y)

  if (nrow(ana) < 6 || dplyr::n_distinct(ana$trt) < 2) {
    return(list(
      est = NA, se = NA, pval = NA,
      ci_lower = NA, ci_upper = NA, converged = TRUE
    ))
  }

  tt <- t.test(change ~ trt, data = ana, var.equal = TRUE)
  est <- unname(tt$estimate[2] - tt$estimate[1])

  list(
    est = est,
    se = unname(diff(tt$conf.int) / (2 * qt(0.975, tt$parameter))),
    pval = tt$p.value,
    ci_lower = tt$conf.int[1],
    ci_upper = tt$conf.int[2],
    converged = TRUE
  )
}

# MMRM is fit with `trt * time` (time entered continuously in the
# fixed-effects formula) rather than `trt * visit`, so that the
# `trt:time` coefficient is a per-unit-time treatment effect directly
# comparable to `delta`, the same estimand SMA-slope targets
# (whitepaper M2, remediation option (b)). The unstructured covariance
# term still groups on the categorical `visit` factor, which is
# unaffected by the fixed-effects parameterization. Note the
# covariance term is unqualified `us(visit | id)`, not
# `mmrm::us(visit | id)`; the namespace-qualified form is not
# recognized by mmrm's formula parser and silently fails on every
# call (whitepaper M1).
fit_mmrm <- function(dat) {
  ana <- dat |>
    dplyr::filter(time > 0, !is.na(y)) |>
    dplyr::mutate(
      visit = factor(time),
      id = factor(id)
    )

  bl <- dat |>
    dplyr::filter(time == 0) |>
    dplyr::select(id, baseline_y = y) |>
    dplyr::mutate(id = factor(id))

  ana <- dplyr::inner_join(ana, bl, by = "id")

  if (nrow(ana) < 10) {
    return(list(
      est = NA, se = NA, pval = NA,
      ci_lower = NA, ci_upper = NA, converged = FALSE
    ))
  }

  tryCatch({
    fit <- mmrm::mmrm(
      y ~ trt * time + baseline_y + us(visit | id),
      data = ana,
      method = "Satterthwaite"
    )

    s <- summary(fit)$coefficients
    idx <- which(rownames(s) == "trt:time")

    if (length(idx) == 0) {
      return(list(
        est = NA, se = NA, pval = NA,
        ci_lower = NA, ci_upper = NA,
        converged = FALSE
      ))
    }

    est_val <- s[idx, 1]
    se_val <- s[idx, 2]
    df_val <- s[idx, 3]
    pval <- s[idx, 4]

    list(
      est = est_val,
      se = se_val,
      pval = pval,
      ci_lower = est_val - qt(0.975, df_val) * se_val,
      ci_upper = est_val + qt(0.975, df_val) * se_val,
      converged = TRUE
    )
  }, error = function(e) {
    list(
      est = NA, se = NA, pval = NA,
      ci_lower = NA, ci_upper = NA, converged = FALSE
    )
  })
}

# Per-method true value of the estimated treatment-effect coefficient
# under the generating model (whitepaper M2). SMA-slope and MMRM
# (fit with continuous time, see `fit_mmrm()`) both target the
# per-unit-time slope effect `delta`. SMA-change targets the
# difference in mean post-baseline change, which under linear
# divergence from baseline equals `delta * (p_visits + 1) / 2`
# (average of `delta * t` over `t = 1, ..., p_visits`).
truth_for_method <- function(method, delta, p_visits) {
  dplyr::case_when(
    method == "SMA-slope" ~ delta,
    method == "SMA-change" ~ delta * (p_visits + 1) / 2,
    method == "MMRM" ~ delta,
    TRUE ~ NA_real_
  )
}

run_one_rep <- function(n_per_arm, p_visits, miss_rate,
                        miss_mech, delta = 0.25, ...) {
  dat <- generate_data(
    n_per_arm = n_per_arm,
    p_visits = p_visits,
    miss_rate = miss_rate,
    miss_mech = miss_mech,
    delta = delta,
    ...
  )

  sma_slope <- fit_sma_slope(dat)
  sma_change <- fit_sma_change(dat)
  mmrm_res <- fit_mmrm(dat)

  out <- dplyr::bind_rows(
    tibble::tibble(
      method = "SMA-slope",
      est = sma_slope$est,
      se = sma_slope$se,
      pval = sma_slope$pval,
      ci_lower = sma_slope$ci_lower,
      ci_upper = sma_slope$ci_upper,
      converged = sma_slope$converged
    ),
    tibble::tibble(
      method = "SMA-change",
      est = sma_change$est,
      se = sma_change$se,
      pval = sma_change$pval,
      ci_lower = sma_change$ci_lower,
      ci_upper = sma_change$ci_upper,
      converged = sma_change$converged
    ),
    tibble::tibble(
      method = "MMRM",
      est = mmrm_res$est,
      se = mmrm_res$se,
      pval = mmrm_res$pval,
      ci_lower = mmrm_res$ci_lower,
      ci_upper = mmrm_res$ci_upper,
      converged = mmrm_res$converged
    )
  )
  out$truth <- truth_for_method(out$method, delta, p_visits)
  out
}

run_simulation <- function(n_per_arm, p_visits, miss_rate,
                           miss_mech, n_reps = 1000,
                           delta = 0.25, rng_dir = NULL, ...) {
  # Morris et al. (2019) §4.1: the RNG seed is set ONCE by the caller;
  # this function does not call set.seed(). Per-replicate RNG states
  # are captured. If `rng_dir` is supplied, the full list of states is
  # written to a sidecar RDS in that directory and the sidecar path is
  # returned as a column on every output row (surviving
  # `purrr::pmap_dfr` row-binding, unlike an attribute, which
  # `pmap_dfr` drops; whitepaper m4). If `rng_dir` is NULL, the states
  # are attached only as an attribute of the returned object, which is
  # sufficient for interactive use but not for the multi-condition
  # `pmap_dfr` loop in `report.Rmd`.
  rng_states <- vector("list", n_reps)

  reps <- purrr::map_dfr(
    seq_len(n_reps),
    function(i) {
      rng_states[[i]] <<- .Random.seed
      run_one_rep(
        n_per_arm = n_per_arm,
        p_visits = p_visits,
        miss_rate = miss_rate,
        miss_mech = miss_mech,
        delta = delta,
        ...
      )
    },
    .id = "rep"
  )

  # Monte Carlo SE formulas per Morris, White & Crowther (2019)
  # Table 6. Degenerate cases (n_valid < 1 or < 2) yield NA_real_
  # rather than NaN/Inf. Bias, MSE, and coverage are scored against
  # each method's own `truth` (whitepaper M2), not a single shared
  # `delta`.
  out <- reps |>
    dplyr::group_by(method, truth) |>
    dplyr::summarise(
      n_valid = sum(!is.na(est)),
      bias = if (n_valid >= 1) {
        mean(est, na.rm = TRUE) - dplyr::first(truth)
      } else NA_real_,
      mcse_bias = if (n_valid >= 1) {
        stats::sd(est, na.rm = TRUE) / sqrt(n_valid)
      } else NA_real_,
      emp_se = if (n_valid >= 2) {
        stats::sd(est, na.rm = TRUE)
      } else NA_real_,
      mcse_emp_se = if (n_valid >= 2) {
        emp_se / sqrt(2 * (n_valid - 1))
      } else NA_real_,
      mse = if (n_valid >= 1) {
        mean((est - dplyr::first(truth))^2, na.rm = TRUE)
      } else NA_real_,
      mcse_mse = if (n_valid >= 1) {
        sqrt(
          stats::var((est - dplyr::first(truth))^2,
                      na.rm = TRUE) / n_valid
        )
      } else NA_real_,
      mean_se = if (n_valid >= 1) {
        mean(se, na.rm = TRUE)
      } else NA_real_,
      mcse_mean_se = if (n_valid >= 1) {
        stats::sd(se, na.rm = TRUE) / sqrt(n_valid)
      } else NA_real_,
      power = if (n_valid >= 1) {
        mean(pval < 0.05, na.rm = TRUE)
      } else NA_real_,
      mcse_power = if (n_valid >= 1) {
        sqrt(power * (1 - power) / n_valid)
      } else NA_real_,
      coverage = if (n_valid >= 1) {
        mean(ci_lower <= dplyr::first(truth) &
               ci_upper >= dplyr::first(truth),
             na.rm = TRUE)
      } else NA_real_,
      mcse_coverage = if (n_valid >= 1) {
        sqrt(coverage * (1 - coverage) / n_valid)
      } else NA_real_,
      convergence = mean(converged, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      n_per_arm = n_per_arm,
      p_visits = p_visits,
      miss_rate = miss_rate,
      miss_mech = miss_mech,
      delta = delta
    )

  if (!is.null(rng_dir)) {
    dir.create(rng_dir, showWarnings = FALSE, recursive = TRUE)
    rng_file <- file.path(
      rng_dir,
      sprintf(
        "rng_n%d_p%d_miss%02d_%s_delta%s.rds",
        n_per_arm, p_visits, round(100 * miss_rate), miss_mech,
        gsub("\\.", "", format(delta))
      )
    )
    saveRDS(rng_states, rng_file)
    out$rng_state_file <- rng_file
  } else {
    attr(out, "rng_states") <- rng_states
  }

  out
}

# Render-blocking sanity checks (whitepaper M3, m9): stop before any
# table or inline number is produced from a simulation run that
# contains impossible or clearly-broken output. This does not
# validate scientific correctness, only that the pipeline did not
# silently fail the way the `mmrm::us` bug (M1) previously did.
check_simulation_quality <- function(results, conv_floor = 0.5) {
  stopifnot(
    "results has zero rows" = nrow(results) > 0,
    "convergence column missing" = "convergence" %in%
      names(results)
  )

  bad_conv <- results |>
    dplyr::filter(convergence < conv_floor)
  if (nrow(bad_conv) > 0) {
    stop(
      "check_simulation_quality: convergence below floor (",
      conv_floor, ") in ", nrow(bad_conv), " method-condition ",
      "cells; inspect fit_mmrm() and the missingness settings ",
      "for the affected conditions before rendering."
    )
  }

  all_na <- results |>
    dplyr::filter(n_valid == 0)
  if (nrow(all_na) > 0) {
    stop(
      "check_simulation_quality: ", nrow(all_na), " method-",
      "condition cells produced zero valid replicates."
    )
  }

  # SMA-slope coverage should sit near the nominal 95% up to Monte
  # Carlo error, in every condition with a nonzero effect (delta > 0)
  # and complete or MCAR data (the setting under which the estimator
  # is unbiased by construction).
  slope_complete <- results |>
    dplyr::filter(
      method == "SMA-slope", delta > 0, miss_mech == "MCAR"
    )
  if (nrow(slope_complete) > 0) {
    bad_cov <- slope_complete |>
      dplyr::filter(
        abs(coverage - 0.95) > 5 * mcse_coverage + 0.03
      )
    if (nrow(bad_cov) > 0) {
      stop(
        "check_simulation_quality: SMA-slope coverage under MCAR ",
        "departs from nominal 0.95 by more than 5 MCSEs + 0.03 in ",
        nrow(bad_cov), " condition(s)."
      )
    }
  }

  invisible(TRUE)
}
