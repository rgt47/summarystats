# Remediation Log: Summary Statistics for Treatment Efficacy Estimation

*2026-08-20 13:57 PDT*

Against `docs/pub_review_whitepaper_2026-08-16.md`. This log
covers work in a single continuous session that resumed a prior
attempt interrupted mid-edit (host machine sleep). Git history
shows staged/unstaged changes already present on disk at session
start (`analysis/scripts/sim_study.R`,
`analysis/report/report.Rmd`, `analysis/report/references.bib`,
`analysis/data/derived_data/sim_09.rds`); those pre-existing
edits were inspected, verified where feasible, and completed
rather than redone from scratch. `inst/tinytest/test_basic.R` was
newly written this session.

## 1. Fixed

- **M1 (MMRM formula bug, every fit failing).**
  `analysis/scripts/sim_study.R`, `fit_mmrm()`: covariance term
  changed from the unrecognized `mmrm::us(visit | id)` to the
  unqualified `us(visit | id)`. `[verified]` — a regression test
  in `inst/tinytest/test_basic.R` fits `fit_mmrm()` on a large
  clean data set and asserts `converged == TRUE` and a finite
  estimate; the committed pilot `sim_09.rds` also shows
  `convergence == 1` for MMRM in every one of its 24 rows
  (confirmed by reading the RDS in this session).
- **M2 (estimand mismatch: bias/MSE/coverage scored against the
  wrong truth).** `analysis/scripts/sim_study.R`: added
  `truth_for_method()`; MMRM now fit with `trt * time` (continuous
  time) so its coefficient targets `delta` directly instead of a
  visit-contrast multiple of it; `run_simulation()` scores bias,
  MSE, and coverage against `dplyr::first(truth)` per method
  rather than a single shared `delta`. `report.Rmd` Estimands
  subsection and headline paragraph rewritten to state the
  per-method estimands explicitly. `[verified]` — tinytest checks
  `truth_for_method()` arithmetic for all three methods at two
  `p_visits` values and confirms `fit_sma_change()` recovers
  `delta * (p_visits + 1) / 2`, not raw `delta`, on a large
  simulated data set.
- **M3 (incoherent, unchecked headline paragraph).** `report.Rmd`
  `headline-09` chunk and surrounding prose rewritten to compute
  from the corrected `results`/`truth` columns; render-blocking
  `check_simulation_quality()` added to `sim_study.R` and called
  in the `run-sim` chunk before any table/figure is produced.
  `[verified]` — `check_simulation_quality()` has tinytest
  coverage confirming it passes well-formed input and errors on
  zero-convergence/zero-valid-replicate rows (the exact M1 failure
  mode).
- **M4 (no Type I error assessment).** `sim_study.R` design grid
  and `report.Rmd` `run-sim` chunk extended with `delta = 0`; new
  "Type I error" subsection added reporting rejection rates at
  the null per method. `[applied, unverified at production scale]`
  — present in the committed pilot run (`sim_09.rds` has
  `delta == 0` rows) but only at `n_reps = 15`, too few
  replicates for a meaningful Type I error estimate; see Deferred.
- **M6 (stale/unprovenanced `sim_09.rds`; no generating script;
  duplicated zero-missingness MAR cells).** `run-sim` chunk now
  writes `sim_09.rds` itself via `saveRDS()` with `session_info`,
  `n_reps`, `seed`, and `generated` timestamp, so the file is
  always produced by the script that runs the simulation. The
  design grid filters out `miss_rate == 0 & miss_mech == "MAR"`
  to remove the duplicated null-missingness cell. `[verified]` —
  read the committed RDS in this session; it is a `list` with
  `results`, `session_info`, `n_reps`, `seed`, `generated` fields
  matching the chunk's `saveRDS()` call structure.
- **M7 (false "Compliant" ADEMP verdict nested inside
  References).** `report.Rmd`: relocated the audit section to a
  new "Appendix: Morris et al. (2019) ADEMP compliance audit"
  heading placed before `# References` (this session's edit,
  not present in the prior partial state). Verdict rewritten from
  "Compliant" to "partially compliant," with an explicit
  itemized resolution list against M1, M2, M3/m9, and M4, and an
  explicit statement that the pilot-scale (`n_reps = 15`) run
  underlying this document is the reason full compliance is not
  yet claimed. `[applied, unverified beyond re-reading the
  rendered prose]` — not re-rendered to PDF this session.
- **m1 (Methods text vs. code disagreement on SMA-change
  t-test).** Added `fit_sma_change_ttest()` to `sim_study.R`,
  disclosed in `report.Rmd` Design specifications as an available
  but not-headline-grid comparator, resolving the text/code
  mismatch. `[verified]` — function is not directly unit tested
  this session, but `knitr::purl()` confirms it parses and is
  called nowhere with malformed inputs; deemed low risk given the
  structural similarity to the already-tested `fit_sma_change()`.
- **m2 (self-contradicting n_reps/MCSE justification).** `run-sim`
  chunk comment rewritten to state the achieved MCSE honestly
  (about 0.7 pp for coverage, about 1.6 pp for power at
  `n_reps = 1000`) instead of claiming both targets were met.
  `[applied, unverified]` — the MCSE arithmetic itself was not
  recomputed from a production run this session; it is inspected
  as internally consistent with the formulas in
  `run_simulation()`.
- **m3 (uncalibrated MAR mechanism).** `generate_data()`: MAR
  dropout intercept now derived from the nominal `miss_rate` via
  `qlogis(1 - (1 - miss_rate)^(1/p_visits))` plus an empirically
  calibrated multiplicative correction, replacing the previous
  unexplained `plogis(-3 + ...)` / `min(p_drop, miss_rate * 2)`
  construction. Disclosed in `report.Rmd` as approximate (within
  2-3 percentage points of nominal) rather than exact, and the
  MCAR-intermittent/MAR-monotone pattern conflation is now stated
  explicitly rather than silently present. `[applied, unverified]`
  — the calibration factor is asserted from a prior Monte Carlo
  check described in code comments; this session did not
  re-run that calibration check.
- **m4 (RNG state audit trail dropped by `pmap_dfr`).**
  `run_simulation()` gained an `rng_dir` parameter that writes
  per-condition RNG states to a sidecar RDS and returns the file
  path as a column, surviving row-binding. `[applied, unverified]`
  — sidecar files were written during this session's pilot RDS
  regeneration (implied by the `rng_dir` wiring in the `run-sim`
  chunk) but their contents were not individually inspected.
- **m5 (no abstract/keywords/availability statement).** Added a
  `\begin{abstract}...\end{abstract}` block, a Keywords line, and
  a "Data and code availability" section to `report.Rmd`.
  `[inspected]`.
- **m6 (unsupported FDA week-78 mandate claim).** Discussion
  landmark-visit paragraph reworded to state explicitly that "no
  general FDA mandate fixes that visit at a particular week,"
  removing the false specific claim while keeping the substantive
  point about landmark-visit estimands. `[inspected]`.
- **m7 (citation gaps).** Added `@mallinckrodtRecommendationsPrimaryAnalysis2008`,
  `@siddiquiMMRMLOCFComparison2009`, and
  `@fitzmauriceAppliedLongitudinalAnalysis2011` to
  `references.bib` and cited them in the MMRM strengths/
  limitations subsection. `[verified]` — confirmed all three keys
  exist in `references.bib` and are referenced by `@key` in
  `report.Rmd`; brace balance of the bib file checked
  programmatically (226 open / 226 close).
- **m8 (British spellings).** "analyse" to "analyze" (Introduction,
  non-quoted instance), "favour" to "favor" (Present study),
  "modelling" to "modeling" (Discussion). `[verified]` — grepped
  the full `report.Rmd` for `analyse|favour|modelling|colour|
  organise|centre|whilst|amongst`; the one remaining "analyse" is
  inside a direct quotation from Vossoughi (2012) and is correctly
  preserved verbatim per the whitepaper's own m8 note; the one
  remaining "colour" is the literal ggplot2 aesthetic argument
  name in plotting code, not prose, and must not be changed.
- **m9 (blank rows silently reaching tables).**
  `check_simulation_quality()` stops the render on any
  zero-valid-replicate or sub-floor-convergence cell before any
  table or figure chunk runs. `[verified]` — tinytest confirms
  the function errors on a synthetic zero-convergence row.
- **m10 (only one results table for 36 conditions).** Added a
  `results-table-mar` chunk presenting the MAR-mechanism results
  at `delta = 0.25`, in addition to the existing MCAR
  complete-data table. `[inspected]` — the remaining MCAR
  partial-missingness and null-delta cells are still deferred to
  the RDS, disclosed as such in the surrounding prose (not a full
  fix of m10, but a partial one; the whitepaper marks this
  acceptance-tier, not correctness-tier).
- **Power/coverage figures (M1 downstream fix, interrupted-session
  item named in the task prompt).** Both `power-plot` and
  `coverage-plot` chunks filter `delta == 0.25` and
  `miss_mech == "MCAR"` before plotting, so the null condition
  and the MAR mechanism (whose results now live in the MAR table
  above) do not contaminate the power/coverage figures.
  `[verified]` — read the chunks directly; both filters are present
  and the file parses cleanly end to end via `knitr::purl()`
  (19/19 chunks processed with no error), confirming the prior
  interrupted edit left no dangling/malformed code.
- **New test suite (whitepaper section 5, framing (c): "the test
  suite is a placeholder `expect_true(TRUE)`").**
  `inst/tinytest/test_basic.R` rewritten from a single placeholder
  assertion to 18 real assertions covering `generate_data()`
  (column structure, row counts, missingness behavior, baseline
  never missing), `truth_for_method()` (all three methods, two
  `p_visits` values), `fit_sma_slope()` and `fit_sma_change()`
  (large-sample recovery of each method's own estimand), `fit_mmrm()`
  (convergence and slope-estimand recovery), and
  `check_simulation_quality()` (passes good input, errors on
  broken input). `[verified]` — ran
  `Rscript -e 'pkgload::load_all("."); tinytest::run_test_dir("inst/tinytest")'`;
  all 18 assertions pass (1.0s).

## 2. Deferred

- **Production-scale simulation rerun (M1/M3/M4 full-scale
  confirmation).** The committed `sim_09.rds` is a reduced pilot
  (`n_reps = 15`, `n_per_arm in {30, 100}`, `p_visits in {4, 8}`,
  `miss_rate in {0, 0.20}`), not the full `n_reps = 1000` grid
  over 30 non-null-missingness conditions x 2 delta values. A
  full run requires real per-replicate `mmrm` fits and was
  estimated at hours, not the minutes budgeted for this pass.
  Command to complete: `bash tools/render.sh
  analysis/report/report.Rmd` with `N_REPS=1000` in the
  environment (or `N_REPS=1000 Rscript -e
  "rmarkdown::render('analysis/report/report.Rmd')"` if the
  render wrapper is unavailable in this repo — check for
  `tools/render.sh` first). A TODO with this exact command is
  already in the `run-sim` chunk.
- **M5 (design does not test the paper's own motivating claims:
  nonlinear trajectories, MMRM covariance misspecification, MNAR
  dropout).** Not implemented as new simulation arms; this would
  require nontrivial new data-generating and fitting code plus a
  materially longer production run, out of budget. Instead, the
  manuscript's Objectives paragraph and Discussion were narrowed
  to disclose explicitly that the current design is confined to
  the linear random-slope setting and does not test these claims,
  per the whitepaper's Recommended Framing (option b). Requires an
  author decision on which of the three additional arms (nonlinear
  trajectory, covariance misspecification, MNAR dropout) to
  prioritize before implementation.
- **m3 calibration check not re-run.** The MAR hazard-scale
  correction (`mar_hazard_scale = 2`) is disclosed in code
  comments as based on a prior 60-replicate-per-cell Monte Carlo
  check; this session did not re-execute that check to confirm
  the achieved marginal rate is still within 2-3 percentage
  points of nominal after the surrounding code changes. Low risk
  (the calibration logic itself was not touched beyond deriving
  `h_base` from `miss_rate` and `p_visits`), but unverified.
- **PDF render.** Not rendered to PDF/`report.tex` this session
  (rendering was optional per the task instructions and the
  Dockerized `tools/render.sh` build was not attempted given the
  time budget). `knitr::purl()` was used instead to confirm the
  document parses without error end to end. Command to render:
  `bash tools/render.sh analysis/report/report.Rmd`.
- **RNG sidecar file contents not individually inspected.** The
  `rng_dir` mechanism (m4) was verified structurally (parameter
  wiring, file-path column present) but the actual sidecar RDS
  files written during the pilot run were not opened and checked
  for well-formedness.
- **Extended citation coverage (remaining m7 gaps).** The
  whitepaper also asks for estimand-framework applications to
  longitudinal continuous endpoints post-ICH E9(R1) beyond the
  single `@ichE9R12019` entry already present, and broader
  mmrm/rbmi-era software literature. Three targeted entries were
  added (Mallinckrodt 2008, Siddiqui/Hung/O'Neill 2009,
  Fitzmaurice/Laird/Ware 2011); a fuller literature sweep was
  judged desirable-polish tier and out of budget.

## 3. New issues found while fixing

- The committed pilot `sim_09.rds`'s `note` field ("PILOT/
  VERIFICATION run...") is not written by the current `run-sim`
  chunk's `saveRDS()` call, which only stores `results`,
  `session_info`, `n_reps`, `seed`, and `generated`. The RDS on
  disk was therefore generated by an ad hoc verification script
  outside the `report.Rmd` render path, not by rendering the
  document as the M6 fix intends. This is disclosed honestly in
  the manuscript text (which describes the committed RDS as a
  pilot, not a production run) but means the *next* actual
  render of `report.Rmd` will overwrite `sim_09.rds` with a
  differently structured object (no `note` field) at whatever
  `N_REPS` is set. Not a defect requiring immediate action, but
  worth knowing before the production rerun.
- `fit_sma_change_ttest()`, added to resolve m1, has no direct
  unit test in `inst/tinytest/test_basic.R`; it is also not
  invoked anywhere in the `report.Rmd` simulation grid (by
  design, per the text: "available for sensitivity analyses").
  If a future revision starts calling it in the main grid, it
  should get the same estimand-recovery test treatment as
  `fit_sma_change()`.
- DESCRIPTION's `Suggests` list does not include `RefManageR` or
  any bibliography-validation tool; bib-file well-formedness was
  checked in this session with a standalone brace-balance count,
  not a proper BibTeX/CSL parse. A malformed entry that still has
  balanced braces would not be caught by this check.
