# Referee Remediation Log: report.Rmd
*2026-11-17 18:10 PST*

Input: `analysis/report/referee-report-2026-11-17.md`. The edits build
on the uncommitted voice revision (`voice-revision-log-2026-11-17.md`)
and are left in the working tree. Nothing was committed, staged, or
pushed.

Evidence labels: **[verified]** means I ran code or a lookup and saw
the result; **[applied, unverified]** means the edit was made but not
checked by execution beyond a successful render.

**Supplementary runs actually executed.** The host load average was
150 to 700 during this session (many concurrent agents), so every run
was reduced. All runs replay the stored RNG states of `sim_09.rds`,
so the supplementary methods see the identical data sets; the replay
was checked to reproduce stored SMA-slope bias and empirical SE to 12
significant digits. mmrm 0.3.19 was used here; production used
0.3.17.

| Output file | Command (from project root) | Cells | Reps | Wall time |
|---|---|---|---|---|
| `sim_09_supp_wls_p8d025.rds` | `SUPP_PART=wls SUPP_P=8 SUPP_DELTA=0.25 SUPP_REPS=300 SUPP_CORES=5 SUPP_TAG=p8d025` | 15 | 300 | 11.3 min |
| `sim_09_supp_wls_p4d025.rds` | `SUPP_PART=wls SUPP_P=4 SUPP_DELTA=0.25 SUPP_REPS=200 SUPP_CORES=4 SUPP_TAG=p4d025` | 15 | 200 | 9.6 min |
| `sim_09_supp_wls_p8d0.rds` | `SUPP_PART=wls SUPP_P=8 SUPP_DELTA=0 SUPP_REPS=250 SUPP_CORES=4 SUPP_TAG=p8d0` | 15 | 250 | 10.6 min |
| `sim_09_supp_kr_null30.rds` | `SUPP_PART=kr SUPP_N=30 SUPP_P=8 SUPP_MISS=MCAR0 SUPP_DELTA=0 SUPP_REPS=300 SUPP_CORES=5 SUPP_TAG=null30` | 1 | 300 | 12.5 min |

Each command is prefixed to `Rscript analysis/scripts/sim_09_supplement.R`.
Two runs went slightly over the 10-minute guideline because of host
load, not because of their size. Key results [verified]: at n = 100,
p = 8, nominal MAR 20 percent, power was 0.31 (SMA-slope), 0.69
(SMA-slope-W, same data sets) and 0.75 (MMRM, main run); at MAR 10
percent, 0.83 (SMA-slope-W) versus 0.82 (MMRM). Weighting halved the
empirical SE (0.196 to 0.098). The weighted null rejection rates were
0.024 to 0.084 at 250 reps (MCSE 0.014), against 0.028 to 0.096 for
the unweighted summary on the same data. Kenward-Roger (one cell
only, 300 replicates, compared unpaired with the 1000-replicate
Satterthwaite rate; suggestive only) at n = 30,
p = 8, complete, null: 0.060 (MCSE 0.014, first 300 data sets) versus
Satterthwaite 0.083 in the main run (all 1000). Achieved MAR
missingness: final visit 0.196 to 0.203 at nominal 0.10 and 0.362 to
0.380 at nominal 0.20; subjects with fewer than two post-baseline
visits 5 to 20 percent. Smoke-test outputs (30 and 10 reps) were
deleted so that they are not bound into the report.

**Render.** `bash tools/render.sh analysis/report/report.Rmd` succeeded
twice (22 pages; no unresolved references; only float-placement
warnings). It updated `report.pdf`, `report.tex`, two figure PDFs, and
`share/MANIFEST.md`, and staged
`share/report-2026-11-17-1821-49bb875-wip.pdf` and
`share/report-2026-11-17-1823-49bb875-wip.pdf`. Both are superseded by
the final render, `share/report-2026-11-17-1829-49bb875-wip.pdf`
(22 pages, no unresolved references), made after all edits. They can
be deleted, along with their `MANIFEST.md` rows.

**Untracked inputs the user must add.** The report now stops (by
design) if the supplementary `wls` or `kr` results are missing. The
four files `analysis/data/derived_data/sim_09_supp_*.rds` are
untracked, as are `run_sim_09.R`, `sim_09_supplement.R` and
`test-supplement.R`. CI rendering will fail until they are committed.
The Data availability section describes them as version-controlled,
which becomes true only on commit. Inline R values were
also checked by purling the document and evaluating every inline
expression [verified].

## 1. Fixed

Completeness: Majors 1, 6 and 7 are fully addressed. Majors 2, 3, 4
and 5 are **partially** addressed: the prose and code defects are
fixed, but part of each remedy (full Kenward-Roger rerun or
size-adjusted power; full-scale weighted-slope run; recalibrated MAR
arm and strong-MAR results; categorical-visit MMRM results and new
scenarios) is deferred below.

### Major comments

- **Major 1 (misattributed citations).** `references.bib`: the three
  records were corrected against Crossref and PubMed, and the cite
  keys were kept so the prose did not have to change keys.
  `ardMixedModelRepeated2011` is now Donohue and Aisen, J Nutr Health
  Aging 2012;16(4):360-364 (PMID 22499459).
  `ardSimulationStudyComparing2015` is now Chen, Ni, Fleisher, Zhou,
  Aisen and Mohs, Alzheimer's Dement TRCI 2018;4:46-53 (doi
  10.1016/j.trci.2017.12.002). The old DOI returned 404 at Crossref.
  `liuMixedModelsCovariateInteractions2021` is now Schuler, Ther Innov
  Regul Sci 2022;56:145-154 (doi 10.1007/s43441-021-00348-y). The old
  DOI resolved to an unrelated 2021 Pharm Stat paper; the old PMID was
  already Schuler's. The attributed findings were rewritten from the
  abstracts. The Introduction AD paragraph no longer says that slope
  models are 'more intuitive and clinically meaningful', because the
  Donohue-Aisen abstract does not say so. (Content claims for both
  papers come from abstracts: Donohue-Aisen read verbatim from PubMed,
  Chen et al. from a Crossref abstract summary. Full texts were not
  read.) Chen et al. are now
  described as finding the slope model more powerful in most
  scenarios while still preferring a relaxed MMRM in general. The
  abstract's 'Ard' became 'Chen et al.'. A whole-bib DOI check over
  Crossref found no other mismatches. Two missing DOIs were added
  (Mallinckrodt 2008, Siddiqui 2009) and the mmrm version used in
  production (0.3.17, from `sim_09.rds` session info) was recorded.
  [verified: Crossref and PubMed E-utilities lookups]
- **Major 2 (Type I error prose contradicts the output).**
  `report.Rmd` Type I error subsection rewritten from the data, with
  inline R only. It now says the SMA rates are consistent with nominal
  up to Monte Carlo error (not conservative). It reports MMRM ranges
  by n at p = 8, the count of null cells above 0.05 + 1.96 MCSE
  (13 of 30), the mean-SE to empirical-SE ratio, and the coverage at
  p = 8. It warns that p = 8 power comparisons are between tests of
  unequal size. The over-parameterization of the 36-parameter
  unstructured covariance against a 4-parameter truth is stated (also
  Minor 4). `mean_se` was added to both tables (Minor 1).
  Kenward-Roger: `fit_mmrm()` gained a `df_method` argument (default
  unchanged); see the supplementary status above for the reduced run.
  [verified: values recomputed from `sim_09.rds`; rendered]
- **Major 3 (MAR power loss unreported).** New Headline paragraph
  reports SMA-slope versus MMRM power at n = 100, p = 8 under complete,
  MCAR 20, MAR 10 and MAR 20, the doubling of the SMA-slope empirical
  SE, and the fact that SMA-slope stays unbiased with nominal coverage
  (bias at most 0.022, coverage 0.944 to 0.963 across MAR cells). The
  mechanism is explained. Added `fit_sma_slope_wls()` (inverse-variance
  weights from a pooled residual variance and a moment estimate of the
  slope variance) and a table comparing it with SMA-slope on the same
  replayed data sets. The adequacy bullet ('<20%') was replaced, and
  the list is now labeled as the authors' judgment (Minor 11). The
  power figure now shows MAR panels (Minor 7). [verified: numbers from
  `sim_09.rds` and the supplementary run; figure rendered]
- **Major 4 (missingness not comparable; calibration claim wrong).**
  Methods rewritten: MCAR rate is per observation; MAR rate roughly
  matches observations missing but final-visit missingness is about
  twice nominal. The '2 to 3 points' claim was removed. A new achieved
  missingness table reads the supplementary output. The confounding of
  MCAR versus MAR contrasts is stated. The misleading code comment in
  `generate_data()` was corrected, and the dropout coefficient was
  exposed as `mar_slope` (default 0.04, unchanged; tested to reproduce
  the production generator exactly). Recalibrating both arms to a
  common quantity was not done (Deferred). [verified: achieved rates
  from replayed data]
- **Major 5 (estimand argument; simulated MMRM not the criticized
  one).** The abstract, the Interpretability item, and the Discussion
  estimand paragraph no longer say the MMRM estimand is conditional on
  the working covariance. The Discussion now states that the estimand
  is fixed by population and intercurrent-event strategy, and that the
  slope estimand depends on design under dropout. An 'Intercurrent
  events' bullet (hypothetical strategy for dropout) was added to
  ADEMP. SMA-change drift under MAR (bias -0.16 to -0.21, 14 to 19
  percent of 1.125) is now interpreted in Results and Discussion. The
  Design section now says the simulated MMRM is a linear-trend model,
  not the categorical-visit landmark MMRM. Added `fit_mmrm_visit()`
  (`trt * visit + baseline * visit`, unstructured, KR, final visit as
  reference level; truth delta * p). The Present study section now
  says the simulation cannot test the robustness and interpretability
  arguments. [verified for code by tinytest; prose rendered]
- **Major 6 (unsupported claims).** 'Exact p-values' qualified to
  complete balanced normal data (Introduction; Discussion rewritten).
  'Convergence failures are common' replaced by a statement tied to
  the mmrm package and to our own convergence result. The SAS-versus-R
  'regulatory risk' claims were removed; the Software item now cites
  what the installed mmrm vignettes say (close agreement; linear KR
  matches SAS for unstructured), which I read locally. 'Robust to
  covariance misspecification' is now attributed to Vossoughi et al.
  with an explicit note that our simulation does not test it. The
  Senn et al. (2000) intention-to-treat claim was deleted from the
  Introduction and the Discussion paragraph was rewritten without
  attributing it, because I could not read the source. [inspected]
- **Major 7 (revision narrative; cache hazard).** All 'whitepaper',
  'earlier version', and 'remediation' text was removed. The
  Discussion estimand-correction note, the run-sim chunk comment, and
  the ADEMP audit appendix moved to `docs/manuscript_revision_history.md`.
  The simulation driver moved verbatim to
  `analysis/scripts/run_sim_09.R`; the report now only reads
  `sim_09.rds` (with `stopifnot` checks) and the global `cache = TRUE`
  was removed. `run_simulation()` was refactored to call a new
  `summarise_reps()`; old and new code gave `identical()` output on a
  test cell. [verified: identity check; replay of stored RNG states
  reproduced stored SMA-slope bias and empirical SE to 12 digits]

### Minor comments

- 1: `mean_se` added to both tables. MSE is stated as kept in the RDS.
  [applied, unverified]
- 2: Runtime statements reconciled. The single remaining statement
  says the run takes at least several hours and that it has not been
  timed on an unloaded host; the inherited '1.5 to 2 hours' figure
  is not verified. [applied]
- 3: '30 non-null-missingness conditions' rewritten as 30
  (n, p, missingness) cells per delta. [applied]
- 4: Over-parameterization note added to the generating-model text.
  [applied]
- 5: Equation now uses delta times a treatment indicator g, with
  random effects indexed by i and g. [applied, unverified (LaTeX)]
- 6: Laird-Ware 'two-stage' wording now separates the two senses.
  [applied]
- 7: Power and coverage figures show all five missingness settings;
  coverage has +/- 1.96 MCSE bars. [applied, unverified (render)]
- 8: Abstract now says 'Monte Carlo standard errors (Morris et al.,
  2019)'. [applied]
- 9: DOIs added; mmrm version recorded; Liu entry moot after Major 1.
  Fitzmaurice entry left (already has `edition = {2}`). [applied]
- 10: Kenward and Roger (2009) cited, with the parameterization point
  stated from an index summary of its abstract, not from the paper.
  [verified citation via Crossref; content inspected only]
- 11: '<20%' threshold removed; list labeled as authors' judgment.
  [applied]
- 12: Future research item 2 now cites Wang, Ogburn and Rosenblum
  (2019) and narrows the gap to unequal-precision slope summaries.
  [applied]

### Missing and questionable references (section 6)

Added after Crossref verification: Kenward and Roger 2009; Wang,
Ogburn and Rosenblum 2019; Lu 2010; Donohue et al. 2023; Raket 2022.
Liang and Zeger 2000 (Sankhya B 62:134-148) and Mallinckrodt et al.
2020 (CRC) were verified by web search of bibliographic indexes and
publisher pages. All are cited in the Discussion or Future research.
Uncited software claims were removed rather than cited (row 9).

**Tests.** `inst/tinytest/test-supplement.R` (16 assertions) covers
the new truths, the unchanged `mar_slope` default, `missingness_summary()`,
WLS = OLS on complete data, smaller WLS SE under MAR, KR convergence and
REML-estimate invariance, and `fit_mmrm_visit()` equal to final-visit
ANCOVA on complete data. `tinytest::run_test_dir('inst/tinytest')`:
36 of 36 pass [verified].

## 2. Deferred

- **Full-scale supplementary runs (Majors 2, 3, 4).** The reduced
  runs above have MCSEs of about 0.02 to 0.03 for power. Run the full
  1000-replicate versions on a quiet host, then delete the tagged
  reduced files so they are not double-counted by `read_supp()`:
  `SUPP_PART=wls Rscript analysis/scripts/sim_09_supplement.R` and
  `SUPP_PART=kr SUPP_P=8 Rscript analysis/scripts/sim_09_supplement.R`
  (roughly 0.15 s per MMRM fit unloaded; about 1.5 to 2 CPU-hours for
  `kr`). The manuscript text will then read the new numbers directly.
- **Size-adjusted power (Major 2).** Not computed: `sim_09.rds` keeps
  only cell summaries, not per-replicate p-values. It needs a rerun
  that saves replicate-level output. Switching the main analysis to
  Kenward-Roger also needs the full rerun
  (`N_REPS=1000 Rscript analysis/scripts/run_sim_09.R` after changing
  `fit_mmrm()`'s default); this is an author decision.
- **Categorical-visit MMRM and strongly informative MAR (Majors 4, 5;
  suggestion 2).** Implemented; `fit_mmrm_visit()` is unit-tested on
  complete data only. Both parts were run end to end through the script as
  throwaway smoke tests (one cell each; `visit`: n = 30, p = 8, MAR 20,
  15 reps, converged 15/15; `strong`: n = 100, p = 8, MAR 20, 10 reps,
  all methods valid) and the outputs were deleted. Note that
  `mar_slope = 0.2` also raises the amount of missingness (final-visit
  missingness about 0.52 in the smoke run), so the strong-MAR arm is
  again confounded by amount; calibrating it is part of the Major 4
  deferral. The manuscript Methods now say
  both are implemented but not yet run. Commands:
  `SUPP_PART=visit Rscript analysis/scripts/sim_09_supplement.R` and
  `SUPP_PART=strong Rscript analysis/scripts/sim_09_supplement.R`.
  The manuscript carries a visible TODO, and the Results prose for
  them still needs to be written once results exist.
- **Recalibrating MCAR and MAR to a common missingness quantity
  (Major 4).** Requires changing the generator and rerunning the whole
  grid. Disclosed in Methods instead; an author decision.
- **Nonlinear-truth and covariance-misspecification scenarios (Major
  5; suggestion 2).** Not built. The Present study and Discussion now
  say the simulation cannot test these arguments.
- **Dawson stratified test (Major 3).** Not implemented; the weighted
  summary was the cheaper of the two remedies the referee named.
- **Senn et al. (2000) intention-to-treat claim (Major 6).** Removed,
  not verified. If the author wants it back, it must be quoted from
  the source.
- **Suggestion 4 (table replacing the adequacy bullets)** and
  **suggestion 1 (reframing around the dropout result)** were done
  only partially: the dropout result now leads the Results and
  Discussion, but the bullet list remains.
- **Contemporariness (section 3).** cLDA, splines, progression models,
  Wang et al., Kenward-Roger 2009 and Mallinckrodt et al. 2020 are now
  cited, but only briefly. The 2023 FDA covariate-adjustment guidance
  and Ye et al. (2023) were not added (not verified in this session).
- **Zotero.** `analysis/report/zotero-sync/missing-references-for-zotero-import.bib`
  still holds the three wrong records. The Zotero library itself was
  not touched and needs the same corrections, or a future sync could
  reintroduce them.
- **Venue fit (section 2).** The referee suggests a tutorial or review
  venue. That is an author decision.

## 3. New issues found while fixing

- **Narrative citations dropped author names.** With the numeric
  `statistics-in-medicine.csl`, every narrative `@key` rendered as a
  bare number at the start of a sentence (for example '19 compared
  MMRM ...'). All 28 narrative citations were rewritten as
  'Author [@key]' (names checked against the corrected bib), and the
  affected paragraphs were rewrapped to 78 columns. [verified in the
  rendered PDF]
- **Any render previously reran the full simulation.** No `cache/`
  directory existed, so a render (including CI's `render-report.yml`)
  would have rerun the 1.5 to 2 hour simulation and overwritten
  `sim_09.rds` at whatever `N_REPS` was set. Fixed by Major 7.
- **SMA-slope under MAR is unbiased; only precision is lost.** Bias
  at most 0.022 and coverage 0.944 to 0.963 across MAR cells. The
  referee's framing ('loses most of its power') is right, but the
  manuscript now states that validity is retained.
- **`fit_sma_slope()` silently drops subjects with fewer than two
  post-baseline visits** (5 to 20 percent under MAR). The weighted
  summary does the same. This is now stated in the text but is
  itself a complete-case restriction worth an author decision.
