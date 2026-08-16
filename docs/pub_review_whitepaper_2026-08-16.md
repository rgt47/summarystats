# Publication Review White Paper: Summary Statistics for Treatment Efficacy Estimation

*Review date: 2026-08-16 10:10 PDT*

Workspace: `~/prj/res/09-summary-stats-efficacy/summarystats`
Manuscript reviewed: `analysis/report/report.Rmd` (rendered
artifacts `report.tex`, `report.pdf`, and
`share/report-2026-08-15-1635-7dbf322.pdf`). Supporting code:
`analysis/scripts/sim_study.R`. Stored results:
`analysis/data/derived_data/sim_09.rds`.

Epistemic status conventions used throughout: "verified" means I
ran code and observed the result; "inspected" means I read the
source or rendered output and confirmed the claim; "inferred"
means consistent with surrounding evidence but not directly
checked; "unverified" means not checked.

## 1. Summary of the work under review

The single manuscript, "Summary statistics for treatment efficacy
estimation in clinical trials: when simplicity rivals
complexity," combines (a) a narrative synthesis of the summary
measures literature (Matthews 1990; Frison and Pocock 1992, 1997;
Senn 2000, 2006; Dawson 1994, 1998; Vossoughi 2012; Ard 2011,
2015) against MMRM, (b) a list of design settings claimed to
favor summary statistics, and (c) a Monte Carlo simulation
comparing SMA-slope (per-subject OLS slope plus ANCOVA on
baseline), SMA-change (post-baseline mean minus baseline plus
ANCOVA), and MMRM (unstructured covariance, Satterthwaite) under
a random-intercept-and-slope generating model across a 3 x 2 x 3
x 2 factorial (n per arm, visit count, missingness rate, MCAR vs
MAR), with ADEMP reporting and Morris Table 6 Monte Carlo SEs.
The literature review is competent and well organized. The
simulation, however, is invalid as rendered: every MMRM fit in
the study failed, and two of the three methods are evaluated
against an estimand they do not target.

## 2. Major issues

### M1. Every MMRM fit failed; the paper's central comparison never happened

Location: `analysis/scripts/sim_study.R`, `fit_mmrm()` (line
152); rendered evidence in `analysis/report/report.tex` lines
582 to 640 and `analysis/data/derived_data/sim_09.rds`.

The model formula is written as
`y ~ trt * visit + baseline_y + mmrm::us(visit | id)`. The
`mmrm` package parses the covariance term by name from the
formula, and the namespace-qualified call `mmrm::us(...)` is not
recognized. Every call throws "Covariance structure must be
specified in formula" and the surrounding `tryCatch` silently
records `converged = FALSE`. I verified this by running
`fit_mmrm`'s formula on generated data: the namespaced form
errors, and the unqualified `us(visit | id)` form fits cleanly
(verified). The consequences, all inspected in the rendered
`report.tex`:

- The results table (Table 2) contains entirely blank MMRM rows
  in all six complete-data cells.
- The inline sentence reads "MMRM convergence was at least 0.000
  in every condition examined," a number that directly
  contradicts the surrounding prose and the paper's thesis.
- The stored results object `sim_09.rds` shows `n_valid = 0` and
  `convergence = 0` for MMRM in every condition (verified by
  reading the RDS).
- The power and coverage figures can contain at most two of the
  three methods (inferred from the data; figures not re-rendered
  by this review).

The paper's stated objective (3) is a comparison of operating
characteristics between summary statistics and MMRM. That
comparison does not exist anywhere in the results. A referee
would treat this as grounds for rejection in the submitted form.
Remediation: change the formula to `us(visit | id)` (or
construct the formula with the covariance term unqualified and
attach the `mmrm` namespace at fit time), rerun the full grid,
and add a hard assertion (e.g., `stopifnot(convergence > 0.9)`
per condition) so that silent wholesale failure can never again
reach a rendered document.

### M2. Estimand mismatch: bias, MSE, and coverage are computed against the wrong truth for two of three methods

Location: `analysis/scripts/sim_study.R`, `run_simulation()`
summary block (bias defined as `mean(est) - delta` with
`delta = 0.25` for all methods); `report.Rmd` Estimands bullet
(ADEMP section) and Design specifications.

The generating model gives treatment a slope effect, so the
group difference at time t is `delta * t`. The three methods
estimate three different quantities:

- SMA-slope estimates `delta` (per-unit-time effect). Correct.
- SMA-change estimates the treatment difference in mean
  post-baseline change, `delta * (p + 1) / 2`, which is 0.625
  for p = 4 and 1.125 for p = 8, not 0.25.
- MMRM, as coded, extracts the `trt:visit{last}` interaction
  coefficient. With treatment contrasts and visit 1 as the
  reference, this is `delta * (p - 1)`, not `delta` and not
  the effect at the final visit (`delta * p`), which would
  require the main effect plus the interaction. The fallback to
  the `trt` main effect when the interaction name is absent
  targets `delta * 1`, a third quantity.

Yet all three are scored against 0.25. The rendered Table 2
confirms the arithmetic: SMA-change "bias" is 0.37 to 0.40 for
p = 4 (0.625 - 0.25 = 0.375) and 0.83 to 0.89 for p = 8
(1.125 - 0.25 = 0.875), with "coverage" as low as 0.434
(inspected). These are not biases; they are the difference
between two different estimands. The headline paragraph then
reports, verbatim, that all three procedures "recovered the true
treatment slope effect delta = 0.25 with bias below 0.939," a
sentence that is self-refuting as rendered (inspected,
`report.tex` lines 584 to 591). Power comparisons remain
interpretable (all three are valid tests of the global null),
but every bias, MSE, and coverage number for SMA-change and MMRM
is meaningless as reported.

Remediation: define the estimand per method and score each
estimator against its own true value (rescale SMA-change by
`2 / (p + 1)`, or equivalently score it against
`delta * (p + 1) / 2`); for MMRM either (a) test the visit-p
contrast `trt + trt:visit{p}` against `delta * p`, or (b) fit
MMRM with continuous time and compare slopes directly. State
explicitly in the Estimands subsection that the methods target
different quantities and how comparability is achieved.

### M3. The reported inline results were never sanity-checked; the headline paragraph is incoherent

Location: `report.Rmd` chunk `headline-09` and the "Headline
findings" paragraph; rendered text in `report.tex` lines 582
to 607.

Beyond M1 and M2, the paragraph as rendered states: bias "below
0.939" while claiming recovery of 0.25; largest absolute bias
1.005 "at 0 percent missing, MAR mechanism" (a duplicated null
condition, see M6); `min(convergence)` printed as 0.000 with no
authorial reaction. The inline-computation design is good
practice, but the absence of any assertion layer or human review
between simulation output and submitted prose is itself a
methodological failure a referee will notice immediately.
Remediation: add automated checks (coverage of SMA-slope within
Monte Carlo error of 0.95, convergence above a floor, bias of
each method near zero against its own estimand) that stop the
render on violation, then rewrite the paragraph from the
corrected results.

### M4. No Type I error assessment

Location: `report.Rmd` Parameter values (Table 1) and
Performance metrics; `sim_study.R` (single `delta = 0.25`).

The simulation runs a single nonzero effect size. There is no
null (`delta = 0`) condition, so Type I error control is never
demonstrated for any method, despite the introduction's claim
that MMRM misspecification "can inflate Type I error" and
despite the paper's recommendation of summary statistics as a
primary analysis. Morris et al. (2019), which the paper adopts
as its reporting standard, treats null scenarios as a core
component of method-comparison studies. A referee for
Biometrics or Statistics in Medicine will require this.
Remediation: add `delta = 0` to the factorial and report
rejection rates with MCSEs; consider a second nonzero effect
size so power curves are not anchored at a single point.

### M5. The simulation design does not test any of the paper's own motivating claims

Location: `report.Rmd` Introduction sections "MMRM: strengths
and practical limitations" and "When are summary statistics
adequate?"; Methods, Data generating model.

The narrative case for summary statistics rests on (a)
robustness to covariance misspecification (Vossoughi 2012), (b)
MMRM small-sample degrees-of-freedom fragility, (c) nonlinear
progression as the setting where MMRM retains an advantage, and
(d) MNAR dropout as a limitation. The simulation varies none of
these. The generating model is exactly the random-slope linear
model under which SMA-slope is known to be nearly efficient
(Frison and Pocock 1997), data are at worst MAR at 20 percent,
and MMRM is fit with the correct-by-construction unstructured
covariance. The design therefore cannot produce a surprise: it
re-verifies a 30-year-old theoretical result and tests nothing
the paper actually argues about. This is the difference between
a publishable simulation and a confirmatory exercise.
Remediation: add at least (i) a nonlinear (e.g., quadratic or
plateau) trajectory arm where the slope summary is misspecified,
(ii) a covariance-misspecification arm for MMRM (e.g., fit
AR(1) or compound symmetry when the truth is random slopes),
and (iii) an MNAR dropout arm. These directly map to the claims
in the Discussion and would let the limitations paragraph rest
on evidence instead of citation.

### M6. The stored derived-data object contradicts the rendered results and has no generating script

Location: `analysis/data/derived_data/sim_09.rds`; `report.Rmd`
Results ("see `analysis/data/derived_data/` for the full
results data frame").

The manuscript directs readers to the derived-data directory
for the full results. The object there has 12 rows (4
conditions x 3 methods) with `n_valid = 500` per cell, while
the rendered document reports 108 scenario-method combinations
at approximately 1,000 replications (verified by reading the
RDS; inspected in `report.tex`). No script in the repository
writes `sim_09.rds` (inspected: `sim_study.R` returns a data
frame but never saves; the `run-sim` chunk relies on the knitr
cache, not the RDS). The pointer is therefore to a stale,
unprovenanced artifact from a different run. Relatedly, the
factorial includes `miss_rate = 0` under both MCAR and MAR,
which are identical conditions; the "MAR, 0 percent missing"
cell cited in the headline is a duplicate of the MCAR null.
Remediation: have the `run-sim` chunk (or a standalone script)
write the full results to a dated RDS with session info; drop
the duplicated zero-missingness MAR cells or collapse them; add
a data README row documenting provenance.

### M7. The Morris compliance section asserts "Compliant" for a study whose results falsify the claim, and it sits inside the References section

Location: `report.Rmd` lines 769 to 796 (`# References`
followed by `## Morris et al. (2019) ADEMP Compliance`);
`docs/morris-audit-2026-04-17.md`.

The subsection declares "Verdict at submission: Compliant."
A study in which 100 percent of one method's fits silently
failed, MCSE-bearing tables contain blank rows, and the
n_reps justification is self-contradicting (see m2 below) is
not compliant with Morris et al.'s reporting standard in any
meaningful sense; the checklist was satisfied syntactically
while the substance failed. Additionally, because the
subsection is nested under the `References` heading, the
compliance text renders between the References title and the
citeproc bibliography (inspected in `report.tex`). Remediation:
move the audit material to an appendix or supplementary
material, re-run the audit after M1 to M6 are fixed, and state
the verdict against the corrected results.

## 3. Minor issues

### m1. Methods text and code disagree on SMA-change

`report.Rmd` Design specifications promises both a two-sample
t-test and a baseline-adjusted ANCOVA for the change score;
`fit_sma_change()` implements only the ANCOVA (inspected).
Align text and code.

### m2. Self-contradicting n_reps justification

The `run-sim` chunk comment states that a power MCSE of at
most 1 percentage point at p = 0.5 requires n_reps of at least
2,500, then sets production n_reps to 1,000 (inspected). Either
raise n_reps or restate the achieved MCSE target honestly
(about 1.6 points at p = 0.5).

### m3. MAR mechanism is uncalibrated and conflated with pattern

The MAR dropout probability
`plogis(-3 + 0.04 * (beta0 - prev_y))`, capped at
`2 * miss_rate`, does not produce a marginal missingness rate
equal to the nominal 10 or 20 percent, so "missing data rate"
means different things across the MCAR and MAR arms
(inspected; not empirically quantified). MCAR is also
intermittent while MAR is monotone dropout, confounding
mechanism with pattern. Calibrate the intercept to hit the
nominal marginal rate and state the pattern difference, or make
both monotone.

### m4. Per-replicate RNG state capture does not survive to the saved output

`run_simulation()` attaches `rng_states` as an attribute of its
return value, but `purrr::pmap_dfr` row-binding in the
`run-sim` chunk drops attributes, so the audit-trail claim
("RNG states stored") is not delivered end to end (inspected).
Store states in a list-column or write them to a sidecar file.

### m5. No abstract, keywords, or data/code availability statement

The manuscript has no abstract and no availability statement;
both are required by every plausible target journal
(inspected).

### m6. Unsupported factual claims in the Discussion

"The FDA-mandated co-primary endpoint at week 78 in
Alzheimer's trials" is asserted without citation and is not a
general FDA mandate; the Mallinckrodt 82 percent concordance
figure is cited but the surrounding inference (LOCF divergence
"cemented MMRM's position") would benefit from the DIA working
group literature (unverified claims; see m7).

### m7. Citation coverage gaps

Missing anchors a referee will expect: the DIA/PhRMA MMRM
position papers (Mallinckrodt et al. 2008; Siddiqui, Hung,
O'Neill 2009), a longitudinal-data text (Fitzmaurice, Laird,
Ware or Diggle et al.) for the two-stage random-effects
lineage, estimand-framework applications to longitudinal
continuous endpoints post ICH E9(R1), and the recent
mmrm/rbmi-era software literature beyond the single `mmrm`
manual entry. Twenty-one entries is thin for a paper whose
first half is a literature synthesis.

### m8. British spellings in prose

"analyse" (Introduction, summary measures paragraph), "favour"
(Present study), "modelling" (Discussion) violate the US
English house style; the Vossoughi quotation's "analyse" is
correctly preserved as verbatim (inspected).

### m9. Blank MMRM rows would survive into any table reusing the current pipeline

`knitr::kable` renders NA as empty strings via
`knitr.kable.NA = ""`, which is how six blank rows reached a
submitted PDF without an error (inspected). After M1, add a
completeness assertion before tabulation.

### m10. Only one results table for 36 conditions

MAR and partial-missingness results are deferred to the (stale,
see M6) RDS. A paper about robustness to missing data must show
the missing-data results, at least in a supplement.

## 4. What remains to be done

Ordered by priority for submission readiness.

Required for correctness:

1. Fix the `mmrm::us` formula bug, rerun the full simulation,
   and confirm nonzero MMRM convergence in every condition (M1).
2. Redefine per-method estimands and rescore bias, MSE, and
   coverage against each method's own true value (M2).
3. Add null (`delta = 0`) conditions and report Type I error
   with MCSEs (M4).
4. Regenerate all inline numbers, tables, and figures from the
   corrected run; add render-blocking sanity assertions (M3,
   m9).
5. Write the full results to a provenanced RDS from a script in
   the repository; remove or document `sim_09.rds`; deduplicate
   the zero-missingness MAR cells (M6).

Required for acceptance:

6. Extend the design to test the paper's motivating claims:
   nonlinear trajectories, MMRM covariance misspecification,
   MNAR dropout (M5).
7. Calibrate the MAR mechanism to the nominal missingness rate
   and disentangle mechanism from pattern (m3).
8. Present missing-data results in tables or a supplement (m10).
9. Add an abstract, keywords, and a data and code availability
   statement (m5).
10. Move the ADEMP audit out of the References section, re-run
    it against corrected results, and report an honest verdict
    (M7).
11. Fill citation gaps (m7) and remove or support the FDA
    week-78 claim (m6).

Desirable polish:

12. Align the SMA-change description with the code or implement
    the promised t-test comparator (m1).
13. Resolve the n_reps versus MCSE-target contradiction (m2);
    consider n_reps = 2,500 for power cells.
14. Deliver the per-replicate RNG state audit trail end to end
    (m4).
15. Correct British spellings (m8).
16. Report hardware, run time, and package versions
    (sessionInfo) for the production run.

## 5. Recommended framing

Plausible framings for this paper:

(a) Methodological comparison via new simulation evidence: a
head-to-head operating-characteristics study of summary
measures versus MMRM. The problem is that the literature
already covers the favorable settings thoroughly: Frison and
Pocock (1992, 1997) established near-optimal efficiency of
adjusted summaries under linear divergence; Vossoughi (2012)
already ran the SMA-versus-LMM simulation including covariance
misspecification; Ard (2015) already compared slope models with
MMRM in the AD setting. A simulation confined to the linear
random-slope world (the current design) replicates known
results and would be judged incremental at best.

(b) Practice-oriented synthesis with decision guidance: a
review plus targeted simulation aimed at trial statisticians
deciding when a summary-measure primary analysis is defensible
under the ICH E9(R1) estimands framework. The estimands angle
is the genuinely underexploited gap: the literature the paper
surveys predates E9(R1), and the observation that the slope
difference is a cleanly defined, covariance-model-free estimand
while the MMRM visit contrast is model-conditional is the most
original sentence in the manuscript. Current practice (MMRM by
default, per the DIA working group lineage) makes a
practitioner-facing "when is simple defensible" paper useful.

(c) Software/tools paper around an R implementation of the
summary-measures workflow (an item already listed in Future
research). Not currently viable: the repository contains no
package code (`R/` is empty, the test suite is a placeholder
`expect_true(TRUE)`; inspected).

Recommendation: framing (b). Reasoning: the simulation, even
after correction, cannot claim novelty as a standalone
methods-comparison because Vossoughi and Ard occupy that ground;
but no surveyed paper positions summary measures inside the
estimands framework or gives practitioners decision rules tied
to E9(R1) attributes and to the documented MMRM software
fragmentation. Implications:

- Title: lead with the decision question and the estimand
  angle, for example "When is a summary-measure primary
  analysis defensible? Estimands, efficiency, and robustness
  for longitudinal trials." Retain the current subtitle idea
  only if the target journal tolerates it.
- Abstract and introduction: open with the E9(R1) estimand
  definition problem, not with MMRM's practical burdens; state
  the decision-guidance contribution explicitly and demote
  "novel comparison" language.
- Comparators: MMRM with correctly and incorrectly specified
  covariance, plus the misspecified-summary (nonlinear truth)
  arm, because the decision guidance is only credible if both
  failure directions are shown (this aligns with the M5
  remediation).
- Target journal: Pharmaceutical Statistics or Statistics in
  Medicine (tutorial/practice article); Clinical Trials is an
  alternative if the estimands and operational-risk material is
  expanded. JASA/Biometrics are not realistic for this content.
- Emphasize: the estimands discussion (currently one paragraph
  in the Discussion; promote it to a section), the
  software-fragmentation and regulatory-risk material, and the
  decision checklist in "When are summary statistics
  adequate?" (promote from bullet list to a structured table
  with evidence citations per row).
- De-emphasize or move to supplement: the full ADEMP audit
  narrative, the complete 36-condition tables (headline
  conditions in text, full grid in supplement), and most of
  the Future research list (retain two or three items).

## 6. Assessment

Verdict: reject in current form, with encouragement to
resubmit after major rework. The literature synthesis is
genuinely serviceable and the reporting infrastructure (ADEMP
structure, MCSE machinery, inline computation) is better than
average. But the empirical core is void: the MMRM arm produced
zero results due to a one-token formula bug, two of three
methods are scored against an estimand they do not target, no
null scenario exists, and the rendered manuscript contains
self-refuting sentences and blank table rows that demonstrate
the results were never read after rendering. Because the
required fixes change every number, table, figure, and
conclusion in the Results, this exceeds "major revision" as the
term is normally used; the simulation must be redesigned (M4,
M5) and rerun before the paper's thesis can be evaluated at
all.

## 7. Revision history

- 2026-08-16: Initial review. No prior whitepaper found in
  `docs/`. Major findings: total MMRM fit failure from the
  `mmrm::us` formula bug (M1); estimand mismatch invalidating
  bias/MSE/coverage for SMA-change and MMRM (M2); incoherent
  rendered headline paragraph (M3); no Type I error assessment
  (M4); design does not test the paper's motivating claims
  (M5); stale unprovenanced `sim_09.rds` (M6); false
  "Compliant" ADEMP verdict nested inside References (M7);
  plus ten minor issues.
