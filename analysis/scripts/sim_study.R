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

  if (miss_rate > 0 && miss_mech == "MAR") {
    for (id_i in unique(dat$id)) {
      rows_i <- which(dat$id == id_i)
      for (j in 2:length(rows_i)) {
        if (is.na(dat$y[rows_i[j - 1]])) {
          dat$y[rows_i[j]] <- NA
          next
        }
        prev_y <- dat$y[rows_i[j - 1]]
        p_drop <- plogis(-3 + 0.04 * (beta0 - prev_y))
        p_drop <- min(p_drop, miss_rate * 2)
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
      y ~ trt * visit + baseline_y +
        mmrm::us(visit | id),
      data = ana,
      method = "Satterthwaite"
    )

    last_visit <- levels(ana$visit)[
      length(levels(ana$visit))
    ]
    coef_name <- paste0("trt:visit", last_visit)

    if (!(coef_name %in% names(coef(fit)))) {
      coef_name <- "trt"
    }

    s <- summary(fit)$coefficients
    idx <- which(rownames(s) == coef_name)

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

  dplyr::bind_rows(
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
}

run_simulation <- function(n_per_arm, p_visits, miss_rate,
                           miss_mech, n_reps = 1000,
                           delta = 0.25, ...) {
  # Morris et al. (2019) §4.1: the RNG seed is set ONCE by the caller;
  # this function does not call set.seed(). Per-replicate RNG states
  # are captured and attached to the returned summary for diagnostic
  # reproducibility of any failing rep.
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
  # rather than NaN/Inf.
  out <- reps |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      n_valid = sum(!is.na(est)),
      bias = if (n_valid >= 1) {
        mean(est, na.rm = TRUE) - delta
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
        mean((est - delta)^2, na.rm = TRUE)
      } else NA_real_,
      mcse_mse = if (n_valid >= 1) {
        sqrt(
          stats::var((est - delta)^2, na.rm = TRUE) / n_valid
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
        mean(ci_lower <= delta & ci_upper >= delta,
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
      miss_mech = miss_mech
    )

  attr(out, "rng_states") <- rng_states
  out
}
