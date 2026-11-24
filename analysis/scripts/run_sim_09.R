# Production driver for the main simulation reported in
# analysis/report/report.Rmd. Writes
# analysis/data/derived_data/sim_09.rds and the per-condition RNG
# state sidecars under analysis/data/derived_data/rng_states/.
#
# This code was moved out of the report's `run-sim` chunk (referee
# report 2026-11-17, Major 7): the knitr cache key did not include
# N_REPS or the contents of sim_study.R, so a cached render could
# silently reuse stale output. The report now only reads sim_09.rds.
# The code below is otherwise unchanged from the chunk that produced
# the committed sim_09.rds (generated 2026-08-20, n_reps = 1000).
#
# Usage, from the project root. Runtime is not verified: the earlier
# manuscript said 1.5 to 2 hours on one core, but per-fit timings
# taken on 2026-11-17 suggest several hours.
#   N_REPS=1000 Rscript analysis/scripts/run_sim_09.R

source(file.path(here::here(), "analysis", "scripts", "sim_study.R"))

# Morris et al. (2019) section 4.1: RNGkind is pinned and the seed
# is set once at the start of the program.
RNGkind("L'Ecuyer-CMRG")
set.seed(20260310)

n_reps <- as.integer(Sys.getenv("N_REPS", "1000"))

design_grid <- expand.grid(
  n_per_arm = c(30, 50, 100),
  p_visits = c(4, 8),
  miss_rate = c(0, 0.10, 0.20),
  miss_mech = c("MCAR", "MAR"),
  delta = c(0, 0.25)
) |>
  # MCAR and MAR are identical at 0% missingness; keep one copy.
  dplyr::filter(!(miss_rate == 0 & miss_mech == "MAR"))

rng_dir <- file.path(
  here::here(), "analysis", "data", "derived_data", "rng_states"
)

results <- purrr::pmap_dfr(
  design_grid,
  function(n_per_arm, p_visits, miss_rate, miss_mech, delta) {
    run_simulation(
      n_per_arm = n_per_arm,
      p_visits = p_visits,
      miss_rate = miss_rate,
      miss_mech = miss_mech,
      delta = delta,
      n_reps = n_reps,
      rng_dir = rng_dir
    )
  }
)

check_simulation_quality(results)

sim_09_path <- file.path(
  here::here(), "analysis", "data", "derived_data", "sim_09.rds"
)
saveRDS(
  list(
    results = results,
    session_info = sessionInfo(),
    n_reps = n_reps,
    seed = 20260310,
    generated = Sys.time()
  ),
  sim_09_path
)
