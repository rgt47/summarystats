# Referee Report: 'Summary statistics for treatment efficacy estimation in clinical trials: when simplicity rivals complexity'
*2026-11-17 17:17 PST*

Manuscript: `analysis/report/report.Rmd` (summarystats compendium).
Target venue assumed: *Statistics in Medicine* (the CSL in the YAML).

## 1. Summary of the manuscript

The manuscript reviews about 30 years of literature on the two-stage
summary-measures approach. In that approach, each subject's repeated
measurements are reduced to a slope or a change score, which is then
analyzed by ANCOVA. The manuscript contrasts this approach with
mixed-model repeated measures (MMRM). It reports a Monte Carlo study
under a random-intercept-and-slope generating model. The design has
n = 30, 50, or 100 per arm, p = 4 or 8 post-baseline visits, 0, 10, or
20 percent missingness under MCAR (intermittent) or MAR (monotone),
and delta = 0 or 0.25, with 1000 replicates per cell. The study
compares SMA-slope ANCOVA, SMA-change ANCOVA, and an MMRM with
continuous time and unstructured covariance. The authors explicitly
do not claim that the simulation is novel. They say the simulation
verifies known theory once each method is scored against its own
estimand. They also argue, on ICH E9(R1) grounds, that the slope
estimand does not depend on a covariance model, whereas the MMRM
visit-specific estimand does.

## 2. Overall assessment and recommendation

**Recommendation: major revision. In its current form the manuscript
does not clear the bar for a top-tier methods journal.**

The simulation code is careful in several respects. It scores each
method against its own estimand, includes a null condition, reports
Morris et al. Monte Carlo SEs, and adds render-blocking sanity
checks. However, the manuscript has four problems that prevent
acceptance:

1. Three citations carry wrong authors, years, or venues (Major 1).
2. The Results prose contradicts the simulation output on Type I
   error (Major 2).
3. The most informative result in the authors' own data goes
   unreported. Under monotone dropout, the unweighted slope summary
   loses most of its power relative to MMRM (Major 3). This result
   undercuts the manuscript's own list of settings in which summary
   statistics are adequate.
4. The estimand argument in the abstract and Discussion is
   technically incorrect as stated (Major 5).

The contribution is a review plus a confirmatory simulation of a
setting already covered by Frison and Pocock (1997), Vossoughi et al.
(2012), and Chen et al. (2018). This fits a tutorial or review venue
(*Pharmaceutical Statistics*, *Therapeutic Innovation & Regulatory
Science*, or a *Statistics in Medicine* tutorial) better than a
primary methods paper.

## 3. Significance and novelty

- **What is new.** The manuscript adds no new theory and no new
  method. The simulation repeats an established comparison (linear
  random-slope truth, unstructured MMRM versus slope summary). The
  manuscript concedes this (Present study, lines 362 to 378). The
  remaining claimed contribution is the ICH E9(R1) decision-support
  framing. However, the Discussion does not develop it into a
  decision framework. Only one paragraph (lines 1014 to 1020)
  addresses estimands, and that paragraph is flawed (Major 5).
- **Overlapping work.** Frison and Pocock (1997) derive the optimal
  linear summary under known covariance. Vossoughi et al. (2012)
  report a Monte Carlo comparison of SMA, LMM, and UMA for linear
  trend data. The paper cited as `ardSimulationStudyComparing2015`
  appears to be Chen et al. (2018, *Alzheimer's & Dementia: TRCI*
  4:46 to 53). That paper compares a slope model with MMRM on
  ADAS-Cog trajectories and finds a power advantage for the slope
  model under near-linearity, with similar Type I error. This is
  essentially the present simulation's conclusion.
- **Contemporariness.** The literature engaged with stops at about
  2012, apart from the misattributed 2021 item. A current referee
  would expect engagement with several bodies of work:
  1. Constrained longitudinal data analysis (cLDA), from Liang and
     Zeger (2000) and Lu (2010).
  2. Spline and progression-model alternatives to MMRM in
     Alzheimer's disease, including Donohue et al. (2023) and Raket
     (2022).
  3. Model-robust covariate adjustment in randomized trials, from
     Wang, Ogburn and Rosenblum (2019), the 2023 FDA covariate
     adjustment guidance, and Ye et al. (2023).
  4. Estimand-aligned estimators and intercurrent-event strategies
     (Mallinckrodt et al., 2020).
  5. The improved Kenward-Roger approximation (Kenward and Roger,
     2009), which mmrm implements directly.

  As it stands, the manuscript reads as situated around 2012.

## 4. Major comments

**Major 1. Misattributed citations (correctness; inspected via
literature search).**
Three bibliography entries do not correspond to the cited works as
described.

- `ardMixedModelRepeated2011` ('Mixed Model of Repeated Measures
  versus Slope Models in Alzheimer's Disease Clinical Trials', Ard
  and Edland, *J Nutr Health Aging* 2011;15(6):437 to 441). A paper
  with this exact title exists, but the indexed record lists M. C.
  Donohue as first author. That record gives 2012, *J Nutr Health
  Aging* 16(4):360 to 364. The bib entry's DOI suffix ('-012-')
  itself suggests 2012.
- `ardSimulationStudyComparing2015` ('A Simulation Study Comparing
  Slope Model with Mixed-Model Repeated Measure...', Ard, Raghavan
  and Edland, *Alzheimer's Research & Therapy* 2015). No 2015 paper
  with this title was found in that journal. The matching paper is
  Chen et al., *Alzheimer's & Dementia: TRCI* 2018;4:46 to 53, doi
  10.1016/j.trci.2017.12.002.
- `liuMixedModelsCovariateInteractions2021` (Liu and Pang,
  *Pharmaceutical Statistics*, doi 10.1002/pst.2165). The work with
  this title is by Alejandro Schuler (arXiv:2108.06621). It appears
  to have been published in *Therapeutic Innovation & Regulatory
  Science* 2022;56:145 to 154 (doi 10.1007/s43441-021-00348-y). I
  did not confirm the journal version's title directly.

These are not stylistic slips. Author names and venues are wrong,
and the manuscript cites these papers repeatedly ('Ard' appears in
the abstract). *Remedy:* replace each entry with the verified record.
Re-check every attributed finding (for example, 'slope models provide
a more intuitive and clinically meaningful way...', line 322) against
the actual paper. Run a DOI-resolution check over the whole `.bib`.

**Major 2. Type I error prose contradicts the simulation output
(correctness; verified by reading `sim_09.rds`, n_reps = 1000).**
Lines 758 to 766 state that 'the SMA methods are slightly
conservative relative to nominal 0.05 throughout'. They also state
that MMRM is 'close to nominal in most conditions but modestly
inflated (up to about 0.10)' only at n = 30, p = 8. The saved
results show otherwise.

- SMA-slope rejection rates range from 0.031 to 0.061, and
  SMA-change rates from 0.039 to 0.070. Many cells exceed 0.05.
  These rates are consistent with nominal size up to Monte Carlo
  error (MCSE about 0.007), not with conservatism.
- MMRM is inflated in every p = 8 cell at n = 30 (0.076 to 0.099)
  and at n = 50 (0.057 to 0.081). At n = 100 it is still 0.059 to
  0.070.
- At p = 8 and delta = 0.25, MMRM coverage is 0.897 to 0.939 for
  n = 30 and 50. Its mean model-based SE underestimates the
  empirical SE by about 5 to 19 percent (12 to 19 percent at
  n = 30; for example, 0.154 versus 0.189 at n = 30, MAR 20
  percent).

A direct consequence concerns the power comparisons. Where MMRM
'wins' at p = 8 (for example, 0.457 versus 0.418 at n = 30,
complete data), the comparison is between tests of unequal size, so
it is not interpretable. *Remedy:* rewrite the Type I error
paragraph from the data. Refit MMRM with Kenward-Roger degrees of
freedom (`method = 'Kenward-Roger'` in `mmrm`). Alternatively, report
size-adjusted power. State the SE underestimation explicitly.

**Major 3. The main substantive result of the simulation is
unreported (inspected and verified).**
Under MAR (monotone) dropout, the unweighted SMA-slope estimator
loses a large share of its power, while MMRM does not. At n = 100
and p = 8, SMA-slope power is as follows:

| Missingness | SMA-slope | MMRM |
|---|---|---|
| Complete | 0.890 | 0.888 |
| MCAR 20 percent | 0.791 | 0.822 |
| MAR 10 percent | 0.437 | 0.825 |
| MAR 20 percent | 0.281 | 0.747 |

The empirical SE of SMA-slope roughly doubles under MAR (0.195
versus 0.091 under MCAR 20 percent). The mechanism is
straightforward. Under monotone dropout, many subjects have only two
or three post-baseline observations. Their OLS slopes have much
larger variance than full-trajectory slopes, and the second-stage
ANCOVA weights them equally. In my check of the generator (100
replicates, n = 100 per arm), 5 to 20 percent of subjects had fewer
than two post-baseline observations under MAR. Those subjects are
silently dropped by `fit_sma_slope()`.

This result contradicts the 'Missing data rates are low to moderate
(<20%)' bullet (line 342) and the Discussion's favorable framing of
summary statistics under missingness (lines 967 to 976). The headline
paragraph reports only a pooled power range and so conceals it.
*Remedy:*

- Report and discuss this result prominently.
- Add an inverse-variance-weighted slope summary, or a
  Dawson-stratified test (Dawson 1994; 1998), as a further method.
  The literature the manuscript cites already prescribes this fix.
- Revise the adequacy conditions accordingly.

**Major 4. Missingness rates are not comparable across mechanisms,
and the MAR calibration claim appears inaccurate (verified on the
generator; 100 replicates per setting).**
For MCAR, `miss_rate` is the per-observation missing probability. My
check gave 9.9 and 19.7 to 19.9 percent of post-baseline
observations missing, and the same fraction at the final visit. For
MAR, the code comment says the hazard is calibrated so that dropout
by the final visit is about `miss_rate`. The achieved rates were:

- Final-visit missingness: 0.20 to 0.21 at nominal 0.10, and 0.36 to
  0.38 at nominal 0.20.
- Post-baseline observations missing: 0.12 to 0.13 and 0.22 to 0.24.

The Methods claim (lines 532 to 535) that achieved missingness was
within 2 to 3 points of nominal holds only for the observation-level
fraction. It does not hold for the dropout-by-final-visit target
that the code says it calibrates. The MCAR and MAR arms therefore
differ simultaneously in mechanism, pattern (the authors disclose
this), and amount. Any MCAR-versus-MAR contrast is confounded.

Separately, the informative component is weak. The dropout logit
coefficient is 0.04 per unit of response, against a between-subject
SD of about 5. *Remedy:* define 'missing data rate' once, calibrate
both arms to the same quantity, and report achieved rates in a table.
Also report results under a more strongly informative MAR mechanism.

**Major 5. The estimand argument is incorrect as stated, and the
simulated 'MMRM' is not the MMRM the paper critiques (inspected).**

- The abstract (lines 77 to 80) and Discussion (lines 1017 to 1020)
  say the MMRM visit-specific estimand 'is conditional on the
  working covariance model'. This is not correct. The estimand, a
  difference in marginal means at a visit, is defined by the
  population and the intercurrent-event strategy, not by the working
  covariance. With a saturated mean model and an unstructured
  covariance, the MMRM estimator is consistent for it under MAR. The
  working covariance affects efficiency, and under MAR it can affect
  validity when the covariance is misspecified. That is a property
  of the estimator, not of the estimand.
- Conversely, the SMA-slope estimand is not free of design
  dependence. Under dropout, a per-subject OLS slope estimand depends
  on the visit schedule and on which visits each subject contributes.
  The manuscript specifies no intercurrent-event strategy for any
  method, as ICH E9(R1) would require.
- The SMA-change bias under MAR (-0.16 to -0.21 at p = 8, MAR 20
  percent; about 15 to 19 percent of the true value 1.125) is a
  symptom of this. The post-baseline mean averages over fewer and
  earlier visits for dropouts, so the target drifts. This is
  estimand drift, not estimator noise. It is reported in the
  headline only as 'largest absolute bias' and is not interpreted.
- The simulated MMRM (`y ~ trt * time + baseline_y + us(visit | id)`)
  is a linear-trend mixed model with a correctly specified mean. It
  is not the categorical-visit MMRM with a landmark contrast that the
  Introduction and Discussion criticize. The simulation therefore
  cannot speak to the interpretability, landmark, or misspecification
  arguments made for summary statistics.

*Remedy:*

- Correct the estimand statements.
- State an intercurrent-event strategy for each method.
- Add the standard `trt * visit` MMRM with its final-visit contrast
  as a fourth method, scored against delta times p.
- Add at least one nonlinear-truth scenario and one
  covariance-misspecification scenario. Without these, the
  manuscript's robustness claims rest on citation alone.

**Major 6. Several claims in the Introduction and Discussion are
unsupported or overstated (inspected).**

- 'Exact p-values from well-understood distributions' (lines 134 and
  958). This holds only for complete, balanced data under normality.
  With unequal numbers of observations per subject, the slope
  summaries are heteroscedastic and the t/F reference is
  approximate.
- 'Convergence failures are common in R packages...' (lines 118 to
  121). This is uncited, and it is contradicted by the authors' own
  result: MMRM convergence was at least 0.998, including n = 30,
  p = 8, with 36 covariance parameters.
- 'Different software implementations yield discrepant p-values' and
  'regulatory risk when results from SAS and R disagree' (lines 115
  to 117 and 963 to 965). These need a specific citation. The mmrm
  package documentation reports agreement with SAS PROC MIXED to
  numerical tolerance.
- 'The summary measure approach is robust to covariance
  misspecification' (lines 951 to 953). The present simulation never
  varies the covariance, and its Results do not support this
  statement.
- The Senn et al. (2000) intention-to-treat claim (lines 189 to 191
  and 967 to 973) should be checked against the source and quoted. I
  did not verify it.

*Remedy:* cite, qualify, or delete each claim.

**Major 7. Material that does not belong in a submitted manuscript,
and a reproducibility hazard (inspected).**

- The manuscript text refers to 'whitepaper M1/M2/M6/m10', 'an
  earlier version of this manuscript', and 'this remediation pass'
  (lines 395 to 399, 511, 664, 817 to 823, 920 to 939, and 1069 to
  1155). The ADEMP audit appendix narrates internal revision history.
  None of this is meaningful to a journal reader, and it should be
  moved to a changelog.
- The `run-sim` chunk is `cache = TRUE`. The knitr cache key does not
  include the `N_REPS` environment variable or the contents of the
  sourced `sim_study.R`. Changing either can silently reuse a stale
  cache, while the Reproducibility note implies otherwise. *Remedy:*
  add `cache.extra = list(Sys.getenv('N_REPS'),
  tools::md5sum(sim_script))`, or read results from `sim_09.rds`
  produced by a standalone script.

## 5. Minor comments

1. Lines 413 to 417 and 672 to 673: 'mean model-based SE' and MSE are
   listed as performance measures, but neither table reports them.
   Add `mean_se` beside `emp_se`, since that comparison is what
   reveals the MMRM SE underestimation in Major 2.
2. Runtime is stated as '1.5-2 hours' (chunk comment, line 597) and
   as 'on the order of an hour or more' (line 916). Reconcile them.
3. Line 547: '30 non-null-missingness conditions' is confusing. The
   count is (n, p, missingness) cells per delta value. Say so.
4. The generating model (lines 427 to 437) has no serial correlation,
   so its true covariance is the 2-parameter-plus-correlation random
   coefficients structure. Note that the unstructured MMRM is
   heavily over-parameterized here (10 or 36 parameters). This bears
   directly on the Type I error result.
5. In the data-generating equation, delta_g is indexed by group, but
   the code uses delta times the treatment indicator. The subscript
   i for random effects should be nested in g. Minor notational
   tidying.
6. Line 226 to 228: 'developed further in the two-stage
   random-effects lineage summarized by Fitzmaurice et al.' Laird and
   Ware is not a two-stage summary-measures method. The phrasing
   conflates two different senses of 'two-stage'.
7. The power figure shows MCAR only. Add MAR panels, which is where
   the methods diverge (Major 3). The coverage figure should show
   Monte Carlo error bands.
8. The abstract's 'Morris-Table-6 Monte Carlo standard errors' is
   jargon for an abstract. Write 'Monte Carlo standard errors
   (Morris et al., 2019)'.
9. Bibliography: `liuMixedModels...` lacks volume and pages (moot
   after Major 1). `mallinckrodtRecommendations2008` and
   `siddiquiMMRMLOCF2009` lack DOIs. `mmrmPackage2024` should record
   the version used. `fitzmaurice...2011` should list edition 2 in
   the style the CSL expects.
10. Lines 255 to 264: the Kenward-Roger parameterization dependence
    was addressed by Kenward and Roger (2009). Cite it and note that
    mmrm offers it.
11. Line 342: the '<20%' threshold has no source. Either cite one or
    present it as the authors' judgment.
12. The Future research items 1, 2, and 7 are reasonable. Item 2
    (sandwich SEs for summary statistics) is not 'unexplored'. HC
    standard errors for ANCOVA in randomized trials are well
    studied, for example by Wang, Ogburn and Rosenblum (2019).

## 6. Missing and questionable references

| Claim location | Issue | Suggested source |
|---|---|---|
| Lines 321, 980; abstract ('Ard') | Misattributed (Major 1) | Donohue MC, et al. Mixed model of repeated measures versus slope models in Alzheimer's disease clinical trials. *J Nutr Health Aging*. 2012;16(4):360-364. (Author list beyond first author not verified.) |
| Lines 324-331, 980 | Misattributed (Major 1) | Chen et al. A simulation study comparing slope model with mixed-model repeated measure to assess cognitive data in clinical trials of Alzheimer's disease. *Alzheimer's & Dementia: TRCI*. 2018;4:46-53. doi:10.1016/j.trci.2017.12.002 |
| Lines 293-298, 1006 | Misattributed (Major 1) | Schuler A. Mixed models for repeated measures should include time-by-covariate interactions to assure power gains and robustness against dropout bias relative to complete-case ANCOVA. *Ther Innov Regul Sci*. 2022;56:145-154. doi:10.1007/s43441-021-00348-y (journal version to confirm; arXiv:2108.06621) |
| Lines 255-264 | KR parameterization issue already addressed | Kenward MG, Roger JH. An improved approximation to the precision of fixed effects from restricted maximum likelihood. *Comput Stat Data Anal*. 2009;53:2583-2595. |
| Lines 157-169, 1014-1020 | Baseline-as-response alternative ignored | Liang KY, Zeger SL. Longitudinal data analysis of continuous and discrete responses for pre-post designs. *Sankhya B*. 2000;62:134-148; Lu K. On efficiency of constrained longitudinal data analysis versus longitudinal analysis of covariance. *Biometrics*. 2010;66:891-896. |
| Lines 978-986 | Modern AD trajectory models absent | Donohue MC, et al. Natural cubic splines for the analysis of Alzheimer's clinical trials. *Pharm Stat*. 2023;22:508-519; Raket LL. Progression models for repeated measures. *Stat Med*. 2022;41:5537-5557. |
| Lines 1014-1020 | Estimand framing lacks a modern source | Mallinckrodt C, Molenberghs G, Lipkovich I, Ratitch B. *Estimands, Estimators and Sensitivity Analysis in Clinical Trials*. CRC Press; 2020. |
| Lines 1034-1038 | 'Unexplored' is incorrect | Wang B, Ogburn EL, Rosenblum M. Analysis of covariance in randomized trials: more precision and valid confidence intervals, without model assumptions. *Biometrics*. 2019;75:1391-1400. |
| Lines 115-121, 963-965 | Uncited software-discrepancy and convergence claims | Cite a specific comparison, or the mmrm validation documentation (Sabanés Bové et al.). |

The bibliographic details in rows 4 to 8 are from my own knowledge
of the literature and were not re-verified by search in this review.
Confirm them before use.

## 7. Suggestions for strengthening the paper

1. Reframe the simulation around the dropout result in Major 3. That
   result is new to the manuscript, practically important, and has a
   known remedy (weighted or stratified summaries). Showing that the
   remedy closes the gap would be a real contribution.
2. Add three scenarios: the categorical-visit MMRM with a landmark
   contrast, a nonlinear (plateau or delayed-effect) truth, and a
   generating covariance with serial correlation. Together these test
   the claims the Introduction makes.
3. Report size-adjusted power, or use Kenward-Roger throughout, so
   that power comparisons are between tests of equal size.
4. Replace the 'settings where summary statistics are adequate'
   bullet list with a table derived from the simulation output.
5. Move the remediation history and audit narrative to `NEWS.md` or
   `docs/`, and read published numbers from a versioned `sim_09.rds`
   rather than a knitr cache.

## 8. Scope of this review

- **Verified (executed):** I read `sim_09.rds` (n_reps = 1000,
  generated 2026-08-20) and tabulated every cell. This is the basis
  for the Type I error, coverage, SE, power, and SMA-change bias
  figures in Majors 2, 3 and 5. I also ran `generate_data()` from
  `analysis/scripts/sim_study.R` (100 replicates per setting,
  n = 100 per arm, seed 1) to measure achieved missingness and the
  fraction of subjects with fewer than two post-baseline observations
  (Majors 3 and 4).
- **Inspected (read):** the full `report.Rmd`, `sim_study.R`, and
  `references.bib`. I checked the truth for each method (delta,
  delta(p+1)/2, delta) against the generating model; it is correct.
  I also checked the MCSE formulas against Morris et al. Table 6
  (correct as coded) and the `fit_mmrm()` column index (column 5,
  correct).
- **Not checked:** I did not re-render the manuscript or rerun the
  production simulation, so I did not confirm that the cache
  reproduces `sim_09.rds`. I did not read the full texts of the
  Senn et al. (2000), Vossoughi et al. (2012), Frison and Pocock,
  Dawson, or Weinberg and Lagakos papers, so the attributed findings
  are inspected for plausibility only. The quoted Vossoughi phrases
  were not verified against the source.
- **Searches run (WebSearch):** the Ard and Edland 2011 title; the
  Ard et al. 2015 simulation title; the Liu and Pang 2021 title and
  the Schuler DOI; the Donohue 2012 *J Nutr Health Aging* record; and
  the Mallinckrodt 2004 duloxetine figures (82.2 percent and 12.4
  percent of 202 comparisons, confirmed). No PubMed-specific tool was
  used.
