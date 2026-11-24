# Manuscript revision history: report.Rmd
*2026-11-17 18:04 PST*

This file holds revision-history material that was removed from
`analysis/report/report.Rmd` on 2026-11-17, in response to Major 7
of `analysis/report/referee-report-2026-11-17.md`. The referee
noted that internal remediation narrative is not meaningful to a
journal reader. The text below is preserved verbatim. The
Discussion note and chunk comment are the committed wording;
the appendix is the voice-revised working-tree wording.

## Discussion note on the estimand correction (removed)

**A note on the estimand correction.** An earlier version of
this manuscript scored the bias, MSE, and coverage of all
three methods against a single value, $\delta = 0.25$, even
though SMA-change and MMRM, as originally coded, target
different functionals of the generating model. That version's
headline paragraph consequently asserted that all three
methods "recovered the true treatment slope effect
$\delta = 0.25$," while its own results table showed
SMA-change bias as large as 0.9, a self-contradiction. The
Estimands subsection and the `truth_for_method()` function in
`sim_study.R` now define and score each method against its own
true value: $\delta$ for SMA-slope and for MMRM (which is fit
with continuous time, so its `trt:time` coefficient targets
$\delta$ directly, rather than a covariance-model-dependent
visit-specific contrast), and $\delta \cdot (p+1)/2$ for
SMA-change. We note that this correction does not change which
method has the best operating characteristics under the linear
random-slope model. It changes what "bias" means for two of
the three methods, from an artifact of the reporting code to a
correctly computed near-zero quantity.

## run-sim chunk comment (removed with the chunk)

The simulation driver moved to `analysis/scripts/run_sim_09.R`.
The chunk comment read:

```r
# Morris et al. (2019) §4.1: RNGkind is pinned and the seed is set ONCE
# at the start of the program. Caller controls RNG state; the
# `run_simulation()` worker does not call `set.seed()`.
RNGkind("L'Ecuyer-CMRG")
set.seed(20260310)

# n_reps: production target is 1000, for which the achieved MCSE is
# about 0.7 percentage points for coverage at 0.95 (n_reps >= 475
# needed for <= 1 pp) but about 1.6 percentage points for power at
# p = 0.5 (n_reps >= 2500 would be needed for <= 1 pp there). This
# text previously claimed both targets were met at n_reps = 1000,
# which is only true for the coverage target; the power MCSE
# statement is corrected here rather than silently left wrong
# (whitepaper m2). Raising n_reps to 2500 for the power-focused cells
# is future work (see Future research).
#
# RESOLVED (remediation 2026-08-20, second pass): this chunk was
# re-rendered at the production n_reps = 1000 over the full
# corrected design grid (30 non-null-missingness conditions x 2
# delta values). That run also exposed and fixed a second bug in
# `fit_mmrm()` (see `sim_study.R`): the p-value was read from
# column 4 of `summary(mmrm_fit)$coefficients`, which is the t
# value, not column 5 (`Pr(>|t|)`); this produced an MMRM Type I
# error near 0.50 at every design condition in the first
# production-scale attempt, undetected in the earlier n_reps = 15
# pilot because the pilot's headline table only reports
# delta = 0.25 cells and its wide Monte Carlo SE masked the
# problem in the delta = 0 cells that were computed. The fix
# (`pval <- s[idx, 5]`) is regression-tested in
# `inst/tinytest/test_basic.R`. `sim_09.rds` below now holds the
# real production-scale output. To regenerate it, run:
#   N_REPS=1000 Rscript -e \
#     "rmarkdown::render('analysis/report/report.Rmd')"
# or, preferably, `bash tools/render.sh analysis/report/report.Rmd`
# with `N_REPS=1000` set in the environment; expect roughly
# 1.5-2 hours on a single host core.
```

## Appendix: ADEMP compliance audit (removed)


### Morris et al. (2019) ADEMP compliance audit

This audit was previously nested under the References
heading, where it rendered between the References title and
the citeproc bibliography (whitepaper M7). It is relocated
here as an appendix. We have re-run it against the corrected
simulation described in this manuscript. We do not reproduce
the pre-remediation 'Compliant' verdict, which was false. At
the time that verdict was written, every MMRM fit had
silently failed (M1). In addition, two of three methods were
scored against the wrong estimand (M2). As a result, the study
did not meet Morris et al.'s reporting standard in any
substantive sense. This was so even though the checklist below
was satisfied syntactically. The original audit document is
retained unmodified at `docs/morris-audit-2026-04-17.md` for
the historical record. The verdict below supersedes it.

**Verdict as of this remediation pass (2026-08-20, second
pass): compliant.** The ADEMP structural elements (aims,
data-generating mechanism, estimands, methods, performance
measures) are present. After the M2 correction, the estimands
subsection correctly states that the three methods target
different functionals. It also states that each method is
scored against its own true value. The ADEMP-recommended null
(`\(\delta = 0\)`) condition is now included (M4). MMRM
converges in every reported cell rather than failing
wholesale (M1). Every reported performance measure carries a
Morris Table 6 Monte Carlo SE. The currently committed
`sim_09.rds` is the production-scale run, with
`n_reps = 1000` across the full 60-condition design grid.
Every inline number, table, and figure in this document
therefore comes from that run. The reported Monte Carlo SEs
are about 0.7 percentage points for coverage. For power they
are about 1.4-1.6 percentage points in the least favorable
cells. These SEs are therefore informative at the precision
that Morris et al.'s reporting standard presumes. This
production run also surfaced a second coding defect not
visible in the pilot. `fit_mmrm()` read the t statistic
rather than the p-value from `mmrm`'s coefficient table
(column 4 instead of column 5). This produced an MMRM Type I
error near 0.50 in an initial production-scale attempt. The
fix is regression-tested (`inst/tinytest/test_basic.R`), and
the run was redone. See the `run-sim` chunk in the Methods for
detail.

**Audit-time gaps and their resolution:**

- *Sim chunks `eval=FALSE`*: Resolved. The `run-sim` chunk is
  `eval=TRUE, cache=TRUE` and runs at render time. Cached
  output is reused on later renders.
- *No `set.seed()` in `sim_study.R`*: Resolved by setting the
  seed once at the caller side (Morris et al. 2019 §4.1.1) in
  the `run-sim` chunk: `RNGkind('L'Ecuyer-CMRG'); set.seed(20260310)`.
- *No Monte Carlo SE machinery*: Resolved. The
  `run_simulation()` function in `sim_study.R` returns Morris
  Table 6 MCSEs for bias, empirical SE, MSE, mean model-based
  SE, power, and coverage. All reported numbers carry the
  corresponding MCSE.
- *No render-blocking sanity checks*: Resolved by
  `check_simulation_quality()`, added during this remediation
  pass. It stops the render if any method-condition cell has
  zero convergence or zero valid replicates. It also stops the
  render if SMA-slope coverage under MCAR departs from nominal
  0.95 by more than 5 MCSEs plus 0.03 (whitepaper M3, m9).
- *MMRM total fit failure (M1)*: Resolved. The covariance term
  in `fit_mmrm()` is now the unqualified `us(visit | id)`. It
  replaces the unrecognized `mmrm::us(visit | id)`. Convergence
  is verified nonzero (at least 0.998 in every reported cell)
  in the production-scale run. The tinytest suite in
  `inst/tinytest/test_basic.R` provides a regression check.
- *Estimand mismatch (M2)*: Resolved in code and text via
  `truth_for_method()` and the corrected Estimands subsection.
- *No Type I error assessment (M4)*: Resolved. `delta = 0` is
  now part of the design grid and reported in the Type I error
  subsection.
- *Production-scale run not yet executed*: Resolved. The
  committed `sim_09.rds` is the `n_reps = 1000` production run
  over the full design grid.
- *MMRM p-value read from the wrong `mmrm::summary()` coefficient
  column (found during the production-scale run)*: Resolved.
  `fit_mmrm()` now reads `Pr(>|t|)` (column 5) rather than
  `t value` (column 4). The fix is regression-tested in
  `inst/tinytest/test_basic.R`.

Pilot runs use `N_REPS=100` for fast smoke-test rendering. (A
smaller ad hoc pilot was used to verify this remediation.)
Production submission requires `N_REPS=1000` (the default)
over the full design grid. This is the setting used to
generate the results reported in this document.

