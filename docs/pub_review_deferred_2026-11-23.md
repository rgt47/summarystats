# Deferred-Item Pass: Summary Statistics for Treatment Efficacy
*2026-11-23 12:45 PST*

Source logs: `docs/pub_review_remediation_2026-08-20.md` (sections 2,
3, and the "Outstanding / deferred after this pass" list of section
4), answering `docs/pub_review_whitepaper_2026-08-16.md`.

Later work checked before acting: commits after 2026-08-20
(`7fca2d3`, `634535d`, `7319a1b`, `7d8c8bc`, `a983d68`, `7f763ae`,
`461618e`) and the uncommitted 2026-11-17 referee pass
(`analysis/report/referee-remediation-log-2026-11-17.md`,
`analysis/scripts/run_sim_09.R`, `analysis/scripts/sim_09_supplement.R`,
`inst/tinytest/test-supplement.R`, four `sim_09_supp_*.rds` files,
and edits to `report.Rmd`, `references.bib`, `sim_study.R`). That
pass is outside the `pub_review_remediation_*` series but resolves
several items below, so it is cited as the resolving source.

## 1. Inventory

| # | Source (08-20 log) | Item | Class |
|---|---|---|---|
| 1 | Sec. 2 | Production-scale (`n_reps = 1000`) rerun | RESOLVED-SINCE |
| 2 | Sec. 2 (M5) | Nonlinear-trajectory, covariance-misspecification, MNAR arms | AUTHOR-DECISION |
| 3 | Sec. 2 (M5), via 11-17 | Categorical-visit MMRM and strong-MAR arms: implemented, not run | LONG-RUN |
| 4 | Sec. 2 and Sec. 4 (m3) | MAR calibration check not re-run; achieved rate not audited | RESOLVED-SINCE |
| 5 | Sec. 2 | PDF render | RESOLVED-SINCE |
| 6 | Sec. 2 (m4) | RNG sidecar contents not inspected | RESOLVED-SINCE |
| 7 | Sec. 2 (m7) | Extended citation coverage | AUTHOR-DECISION (partly resolved) |
| 8 | Sec. 3 | Pilot RDS `note` field differs from render output | NO-LONGER-APPLIES |
| 9 | Sec. 3 and Sec. 4 | `fit_sma_change_ttest()` has no unit test | DO-NOW |
| 10 | Sec. 3 | Bib checked only by brace balance | DO-NOW |
| 11 | Sec. 4 | MMRM Type I inflation (up to 0.099) not checked against Kenward-Roger | LONG-RUN (partly resolved) |
| 12 | Sec. 4 | Whether to switch the main MMRM analysis to Kenward-Roger | AUTHOR-DECISION |

Totals: 12 items; 4 resolved-since; 1 no-longer-applies; 2 done now;
2 long-run; 3 author-decision.

Notes on the RESOLVED-SINCE and NO-LONGER-APPLIES rows:

- **Item 1.** Done in the second pass of the same 08-20 log
  (section 4) and committed in `7fca2d3` / `634535d`. Confirmed today
  by reading `analysis/data/derived_data/sim_09.rds`: `n_reps = 1000`,
  180 rows (60 cells x 3 methods). [verified]
- **Item 4.** The 11-17 pass measured achieved missingness by
  replaying the stored data sets (final-visit missingness 0.196 to
  0.203 at nominal 0.10 MAR and 0.362 to 0.380 at nominal 0.20 MAR),
  added Table `tab:achieved` to `report.Rmd`, removed the
  'within 2-3 percentage points of nominal' claim, and corrected the
  `generate_data()` comment (`sim_study.R` lines 37 to 40 now say the
  scale factor does not calibrate final-visit missingness). Confirmed
  today by grep: no '2-3' or 'percentage points of nominal' text
  remains in `report.Rmd` or `sim_study.R`. [inspected] Recalibrating
  MCAR and MAR to a common quantity is a separate author decision
  carried in the 11-17 log; it is not re-listed here.
- **Item 5.** Rendered in the 08-20 second pass and again on
  2026-11-17 through `tools/render.sh`. `report.pdf` (Nov 17 18:30)
  is newer than `report.Rmd` (18:29), `references.bib` (18:28) and
  `sim_study.R` (17:42), so no re-render is needed. [inspected]
- **Item 6.** The 08-20 second pass checked that each of the 60
  sidecars has 1000 states, and the 11-17 replay reproduced the
  stored SMA-slope bias and empirical SE to 12 significant digits.
  Confirmed today with a script: 60 files, each a list of 1000
  `integer(7)` L'Ecuyer-CMRG states, all distinct within file, and
  all 60 files referenced by `sim_09.rds$results$rng_state_file`
  exist. [verified]
- **Item 8.** The pilot RDS has been overwritten by the production
  run; `names(sim_09.rds)` is now `results, session_info, n_reps,
  seed, generated` (no `note`). The driver also moved out of the
  report into `analysis/scripts/run_sim_09.R`, and the report only
  reads the RDS, so a render no longer overwrites it. [verified]

## 2. Done now

- **Item 9: unit tests for `fit_sma_change_ttest()`.**
  `inst/tinytest/test_basic.R` gained five assertions on the
  `n_per_arm = 400` complete-data set: the estimate, SE and p-value
  equal those of the unadjusted `lm(change ~ trt)` fit (which checks
  the treated-minus-control sign convention of
  `tt$estimate[2] - tt$estimate[1]` and the SE back-calculated from
  the confidence interval), the estimate is within 0.2 of
  `delta * (p_visits + 1) / 2`, and the function returns `NA` when
  only one arm remains. [verified]
  `Rscript -e "pkgload::load_all('.'); tinytest::run_test_dir('inst/tinytest')"`:
  41 of 41 pass (25 in `test_basic.R`, previously 20; 16 in
  `test-supplement.R`).
- **Item 10: proper bibliography parse.** Parsed
  `analysis/report/references.bib` with
  `pandoc -f biblatex -t csljson`: 31 entries parsed, 0 warnings,
  31 unique keys, no duplicate keys, no key that pandoc failed to
  parse, and no `@key` in `report.Rmd` that is absent from the bib.
  [verified] This is a parse check, not a check that each record's
  content is correct; the 11-17 pass did a Crossref DOI check of the
  whole file.

No manuscript, figure, or result file was changed in this pass.

## 3. Long runs

All commands run from the project root
(`/Users/zenn/prj/res/09-summary-stats-efficacy/summarystats`).
Wall times below are extrapolated from the reduced runs logged on
2026-11-17 on a heavily loaded host; they were not re-timed today.

- **Item 11: full Kenward-Roger check of the continuous-time MMRM
  over all `p = 8` cells (30 cells x 1000 replays).**
  `SUPP_PART=kr SUPP_P=8 SUPP_CORES=4 Rscript analysis/scripts/sim_09_supplement.R`.
  The only existing measurement is one cell (n = 30, p = 8,
  complete, null) at 300 replays in 12.5 min with 5 workers under
  load, which extrapolates to several hours of wall time for the
  full set; the 11-17 log estimates 1.5 to 2 CPU-hours unloaded.
  Afterwards: delete `analysis/data/derived_data/sim_09_supp_kr_null30.rds`
  so that `read_supp()` does not bind it alongside the full file,
  re-render with `bash tools/render.sh analysis/report/report.Rmd`,
  and remove the TODO at the 'Kenward-Roger inference for MMRM'
  paragraph (around `report.Rmd` line 1206); the inline values read
  the new file directly. Check that the Type I error prose
  (13 of 30 null cells above 0.05 + 1.96 MCSE under Satterthwaite)
  still reads correctly against the KR rates.
- **Item 3: categorical-visit MMRM and strong-MAR arms.**
  `SUPP_PART=visit SUPP_CORES=4 Rscript analysis/scripts/sim_09_supplement.R`
  (30 cells, categorical-visit unstructured MMRM with KR; slower
  per fit than the continuous-time model, so expect at least the
  `kr` run time) and
  `SUPP_PART=strong SUPP_CORES=4 Rscript analysis/scripts/sim_09_supplement.R`
  (3 cells; well under an hour). Afterwards: write the Results prose
  for both (none exists yet), remove the 'implemented but not yet
  run' TODO (around `report.Rmd` lines 453, 583 and 1211), and
  re-render.
- **Related, from the 11-17 log (not a 08-20 item):** the full
  weighted-slope run `SUPP_PART=wls SUPP_CORES=4 Rscript
  analysis/scripts/sim_09_supplement.R` (60 cells x 1000; the reduced
  runs took about 10 min per 15 cells x 200 to 300 replays, which
  extrapolates to roughly 2 to 3 hours under load). Delete the three
  tagged `sim_09_supp_wls_*.rds` files afterwards.

The script's default `SUPP_CORES` is 6; on a shared host pass a
smaller value explicitly.

## 4. Author decisions

- **Item 2: which of the untested motivating claims to simulate.**
  Options: (a) add a nonlinear-trajectory arm (tests the paper's
  claim that slopes misrepresent curved decline), (b) add a
  covariance-misspecification arm (tests the Vossoughi robustness
  claim the paper attributes), (c) add an MNAR arm, or (d) keep the
  current scope, which the manuscript already discloses (Present
  study and Discussion say the simulation cannot test these
  arguments). Evidence: the headline result is that SMA-slope is
  unbiased but loses precision under MAR; none of the three arms
  exists in code. Recommendation: (a) first, since a nonlinear
  truth directly threatens the slope estimand that the paper
  recommends and is the cheapest to add to `generate_data()`;
  defer (b) and (c).
- **Item 7: remaining citation gaps.** The 11-17 pass added
  Kenward-Roger 2009, Wang-Ogburn-Rosenblum 2019, Lu 2010,
  Donohue et al. 2023, Raket 2022, Liang-Zeger 2000 and Mallinckrodt
  et al. 2020, and recorded the mmrm version. Not added: the FDA
  2023 guidance on covariate adjustment in randomized trials, Ye et
  al. (2023) on covariate adjustment, and rbmi-era reference-based
  imputation literature. Options: add them in the Discussion with
  verified records, or leave the literature as is. Recommendation:
  add the FDA 2023 guidance (it bears directly on the baseline
  adjustment used by SMA-change) after verifying the record; the
  others are optional for a methods note of this scope. Not done
  here because the entries need external verification and a choice
  of where the argument goes.
- **Item 12: Kenward-Roger as the primary MMRM inference.**
  Options: keep Satterthwaite (matches `sim_09.rds`; disclose the
  inflation, as now) or switch `fit_mmrm()`'s default to
  Kenward-Roger and rerun the main grid
  (`N_REPS=1000 Rscript analysis/scripts/run_sim_09.R`, several
  hours) so that power comparisons at p = 8 are between tests of
  similar size. Evidence: one-cell KR null rate 0.060 (MCSE 0.014,
  300 replays) versus Satterthwaite 0.083 (1000 replays) in the same
  cell; unpaired and suggestive only. Size-adjusted power is not
  computable from `sim_09.rds` because it stores only cell
  summaries. Recommendation: run Item 11 first; if KR rates are
  near nominal across the p = 8 cells, switch the main analysis and
  save replicate-level p-values in the rerun.

## 5. New issues found

- **CI will fail on a clean checkout until the 11-17 files are
  committed.** `report.Rmd` (around line 1078) stops if the
  supplementary `wls` or `kr` results are missing, and
  `sim_09_supp_*.rds`, `run_sim_09.R`, `sim_09_supplement.R` and
  `test-supplement.R` are untracked. This is noted in the 11-17 log
  but not in the `pub_review` series; it needs the author's commit.
  [inspected]
- **Stale code comment.** `analysis/scripts/sim_study.R` lines 179
  to 183 say the t-test variant 'is not implemented here', then
  point to `fit_sma_change_ttest()` immediately below, which does
  implement it. Cosmetic; left unchanged to keep this pass to the
  deferred items. [inspected]
- **Superseded render copies.** `analysis/report/share/` holds two
  superseded 11-17 renders (`report-2026-11-17-1821-*` and
  `-1823-*`) and their `MANIFEST.md` rows alongside the final
  `-1829-*` copy. Not deleted here (no deletion of result files).
  [inspected]
