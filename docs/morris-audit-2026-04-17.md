# Morris et al. (2019) ADEMP Audit: 09-summary-stats-efficacy
*2026-04-17 09:02 PDT*

## Scope

Files audited:

- `analysis/scripts/sim_study.R`
- `analysis/report/report.Rmd`

## ADEMP scorecard

| Criterion | Status | Evidence |
|---|---|---|
| Aims explicit | Partial | literature-review preamble; no ADEMP header |
| DGMs documented | Partial | DGM described in report text; `sim_study.R` implements it |
| Factors varied factorially | Partial | scenario grid implicit |
| Estimand defined with true value | Met | within-individual effect parameterised |
| Methods justified | Met | summary-stats, paired-t, mixed model compared |
| Performance measures justified | Partial | power and coverage listed in report narrative |
| n_sim stated | Partial | script sets a default, but simulation never runs in the rendered report |
| n_sim justified via MCSE | Not met | no MCSE derivation |
| MCSE reported per metric | Not met | not computed |
| Seed set once | **Not met** | no `set.seed()` call anywhere in `sim_study.R` |
| RNG states stored | Not met | not stored |
| Paired comparisons | Met | same dataset fed to all three analyses per rep |
| Reproducibility | **Not met** | simulation chunks are `eval=FALSE`; report never runs the sim |

## Overall verdict

**Not compliant.**

## Gaps

- Simulation chunks in the report are `eval=FALSE`:
  - `run-sim` chunk at `report.Rmd:405`
  - `power-plot` and `coverage-plot` chunks
- Placeholder text `*[Results to be populated...]*` at `report.Rmd:453`
  indicates results have never been produced.
- `sim_study.R` contains **no** `set.seed()` call — results
  non-reproducible even if the chunks were activated.
- No MCSE machinery.
- `RNGkind()` not pinned.

## Remediation plan

1. Add `set.seed(<fixed integer>)` and `RNGkind("L'Ecuyer-CMRG")` at the
   top of `analysis/scripts/sim_study.R`. Propose `20260310` (YYYYMMDD
   convention used in sibling repos).
2. Activate the `run-sim` chunk in `report.Rmd` (set `eval=TRUE`); cache
   its output to an RDS in `analysis/data/derived_data/`.
3. Activate `power-plot` and `coverage-plot` chunks.
4. Remove the placeholder text and wire the plots to the cached sim RDS.
5. Add `mcse_*` columns (bias, empirical SE, coverage, power) to the
   simulation summary per Morris Table 6.
6. Add an n_sim justification paragraph deriving the chosen n_rep from a
   target MCSE.
7. Add ADEMP-headed Methods section to `report.Rmd` that cites Morris
   Table 6.
8. Store `.Random.seed` per replicate for diagnostic reproducibility of
   any failing rep.

## References

Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate
statistical methods. Stat Med 2019;38:2074-2102. doi:10.1002/sim.8086

---
*Source: ~/prj/res/09-summary-stats-efficacy/summarystats/docs/morris-audit-2026-04-17.md*
