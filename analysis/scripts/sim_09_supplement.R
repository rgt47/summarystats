# Supplementary analyses requested by the referee report of
# 2026-11-17. Each replicate data set of the production run
# (sim_09.rds) is regenerated exactly by restoring the RNG state that
# run_simulation() saved before that replicate (rng_states/ sidecars),
# and the new methods are fitted to the identical data. The replay
# was checked to reproduce the stored SMA-slope bias and empirical SE
# to 12 significant digits.
#
# Parts (environment variable SUPP_PART):
#   wls     All 60 cells. SMA-slope (replay check), the
#           inverse-variance-weighted slope summary SMA-slope-W, and
#           achieved missingness (Majors 3 and 4).
#   kr      All p = 8 cells. Continuous-time MMRM with Kenward-Roger
#           degrees of freedom, MMRM-KR (Major 2).
#   visit   All p = 8 cells. Categorical-visit MMRM with a final-visit
#           contrast and Kenward-Roger degrees of freedom, MMRM-visit
#           (Major 5).
#   strong  n = 100, p = 8, MAR cells regenerated with a strongly
#           informative dropout coefficient (mar_slope = 0.2 instead
#           of 0.04); SMA-slope, SMA-slope-W, MMRM (Major 4).
#
# Usage, from the project root (wall-clock times on a 4+4 core
# laptop with SUPP_CORES = 6 and SUPP_REPS = 1000 are logged in
# analysis/report/referee-remediation-log-2026-11-17.md):
#   SUPP_PART=wls Rscript analysis/scripts/sim_09_supplement.R

source(file.path(here::here(), "analysis", "scripts", "sim_study.R"))

part <- Sys.getenv("SUPP_PART", "wls")
n_cores <- as.integer(Sys.getenv("SUPP_CORES", "6"))
supp_reps <- as.integer(Sys.getenv("SUPP_REPS", "1000"))

derived <- file.path(here::here(), "analysis", "data", "derived_data")
main <- readRDS(file.path(derived, "sim_09.rds"))

cells <- main$results |>
  dplyr::distinct(
    n_per_arm, p_visits, miss_rate, miss_mech, delta, rng_state_file
  ) |>
  dplyr::mutate(
    rng_state_file = file.path(
      derived, "rng_states", basename(rng_state_file)
    )
  )

cells <- switch(
  part,
  wls = cells,
  kr = dplyr::filter(cells, p_visits == 8),
  visit = dplyr::filter(cells, p_visits == 8),
  strong = dplyr::filter(
    cells, n_per_arm == 100, p_visits == 8, miss_mech == "MAR"
  ),
  stop("unknown SUPP_PART: ", part)
)

# Optional cell filters, for reduced runs on a busy host, e.g.
#   SUPP_N=30 SUPP_P=8 SUPP_MISS=MCAR0 SUPP_DELTA=0 SUPP_TAG=null30
# SUPP_MISS takes labels such as MCAR0, MCAR10, MAR20 (comma list).
# SUPP_TAG is appended to the output file name so that reduced runs
# do not overwrite one another; the report binds all tagged files.
env_num <- function(var) {
  v <- Sys.getenv(var, "")
  if (v == "") NULL else as.numeric(strsplit(v, ",")[[1]])
}
f_n <- env_num("SUPP_N")
f_p <- env_num("SUPP_P")
f_delta <- env_num("SUPP_DELTA")
f_miss <- Sys.getenv("SUPP_MISS", "")
tag <- Sys.getenv("SUPP_TAG", "")
cells <- cells |>
  dplyr::mutate(
    miss_label = paste0(miss_mech, round(100 * miss_rate))
  ) |>
  dplyr::filter(
    is.null(f_n) | n_per_arm %in% f_n,
    is.null(f_p) | p_visits %in% f_p,
    is.null(f_delta) | delta %in% f_delta,
    f_miss == "" | miss_label %in% strsplit(f_miss, ",")[[1]]
  )
stopifnot("no cells selected" = nrow(cells) > 0)

mar_slope <- if (part == "strong") 0.2 else 0.04

fit_methods <- function(dat) {
  pack <- function(method, f) {
    tibble::tibble(
      method = method, est = f$est, se = f$se, pval = f$pval,
      ci_lower = f$ci_lower, ci_upper = f$ci_upper,
      converged = f$converged
    )
  }
  switch(
    part,
    wls = dplyr::bind_rows(
      pack("SMA-slope", fit_sma_slope(dat)),
      pack("SMA-slope-W", fit_sma_slope_wls(dat))
    ),
    kr = pack("MMRM-KR", fit_mmrm(dat, df_method = "Kenward-Roger")),
    visit = pack("MMRM-visit", fit_mmrm_visit(dat)),
    strong = dplyr::bind_rows(
      pack("SMA-slope", fit_sma_slope(dat)),
      pack("SMA-slope-W", fit_sma_slope_wls(dat)),
      pack("MMRM", fit_mmrm(dat))
    )
  )
}

run_cell <- function(i) {
  cell <- cells[i, ]
  states <- readRDS(cell$rng_state_file)
  n_use <- min(supp_reps, length(states))
  miss <- matrix(NA_real_, n_use, 3)
  reps <- purrr::map_dfr(seq_len(n_use), function(r) {
    assign(".Random.seed", states[[r]], envir = globalenv())
    dat <- generate_data(
      n_per_arm = cell$n_per_arm,
      p_visits = cell$p_visits,
      miss_rate = cell$miss_rate,
      miss_mech = as.character(cell$miss_mech),
      delta = cell$delta,
      mar_slope = mar_slope
    )
    miss[r, ] <<- missingness_summary(dat)
    out <- fit_methods(dat)
    out$truth <- truth_for_method(out$method, cell$delta, cell$p_visits)
    out
  })
  summarise_reps(reps) |>
    dplyr::mutate(
      n_per_arm = cell$n_per_arm,
      p_visits = cell$p_visits,
      miss_rate = cell$miss_rate,
      miss_mech = as.character(cell$miss_mech),
      delta = cell$delta,
      n_reps = n_use,
      ach_obs_missing = mean(miss[, 1]),
      ach_final_missing = mean(miss[, 2]),
      ach_subj_lt2 = mean(miss[, 3])
    )
}

t0 <- Sys.time()
res_list <- parallel::mclapply(
  seq_len(nrow(cells)), run_cell, mc.cores = n_cores
)
failed <- vapply(res_list, inherits, logical(1), what = "try-error")
if (any(failed)) {
  stop("cells failed: ", paste(which(failed), collapse = ", "),
       "\n", res_list[[which(failed)[1]]])
}
results <- dplyr::bind_rows(res_list)
elapsed <- difftime(Sys.time(), t0, units = "mins")

stopifnot(
  "supplement produced no rows" = nrow(results) > 0,
  "a method-cell has zero valid replicates" = all(results$n_valid > 0)
)

saveRDS(
  list(
    results = results,
    part = part,
    mar_slope = mar_slope,
    supp_reps = supp_reps,
    source = "replay of sim_09.rds RNG states",
    elapsed_min = as.numeric(elapsed),
    session_info = sessionInfo(),
    generated = Sys.time()
  ),
  file.path(
    derived,
    paste0("sim_09_supp_", part, if (tag != "") paste0("_", tag), ".rds")
  )
)
cat("part", part, "done in", round(as.numeric(elapsed), 1), "min\n")
