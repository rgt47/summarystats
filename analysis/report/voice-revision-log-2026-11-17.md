# Voice Revision Log: report.Rmd
*2026-11-17 17:17 PST*

Scope: a prose-only revision of `analysis/report/report.Rmd` into the
`/use_rgtvoice` register. The edits are in the working tree only; none
have been committed, staged, or rendered.

**Stricter sentence rule applied.** The revision follows the updated
rule '2. Short, simple sentences' in `use_rgtvoice.md`:

- one idea per sentence;
- a target of 15 to 20 words, with a ceiling of 25;
- subordinate clauses split out, using a plain connective such as
  'This means', 'As a result', or 'However,';
- the common word preferred over the Latinate one.

Content and hedges were kept; long sentences were split, not cut.

## Kinds of changes made (48 Edit operations)

1. **Sentence splitting.** Most edits are of this kind. Long sentences
   with stacked subordinate or participial clauses were divided at
   clause boundaries, with an explicit connective at the head of the
   new sentence. Examples are the Frison and Pocock (1997) paragraph,
   the Senn et al. (2000) paragraph, the MMRM design item, the MAR
   mechanism bullet, the Estimands bullet, the Discussion estimand
   note, and the Appendix verdict.
2. **Plainer wording.** 'demonstrated' became 'showed'; 'requires
   specification of' became 'requires the analyst to specify';
   'requires estimation of' became 'requires estimating'; 'exhibit'
   became 'show'; 'compared ... against' became 'compared ... with';
   'disentangling' became 'separating'; 'subsequent' became 'later'.
3. **Signposting in voice.** 'We begin with' now opens the headline
   results. The research question in the Introduction is posed as a
   question and then answered ('We ask a simple question. Are there
   ...?'). 'That is,' and 'As a result,' carry the logic that the
   removed subordinate clauses used to carry.
4. **Mechanical 'We note that' trimmed.** Two or three of these were
   removed where they introduced a plain fact rather than an
   observation, for example in the Data availability section and in
   the Vossoughi paragraph. They were kept elsewhere.
5. **Em dashes removed.** All nine em dashes in the Appendix
   gap-resolution list became colons (global rule). The manuscript
   now contains none.
6. **Quotation style.** Inside the edited Discussion estimand note,
   double-quoted phrases now use single quotes (global rule).
7. **One paragraph break added** in 'Present study', separating the
   objectives from the scope statement. There is no change to the
   heading structure.

## Verification (structural checks, not a render)

I compared the edited file against `HEAD` with a script:

- Code chunks (all chunk text and chunk options): identical.
- Inline R expressions: identical, in the same order.
- Citation keys: identical as a multiset.
- Prose LaTeX math and display equations: identical as a multiset.
- Numbers in prose: identical as a multiset.
- YAML header, abstract, and every heading: identical.

The document was not rendered (instructed), so LaTeX compilation of
the edited prose is unverified. That said, no LaTeX commands were
added or removed.

## Passages deliberately left alone

- **Abstract.** It is excluded by the skill and by instruction.
- **Code chunks, chunk comments, table and figure content, captions
  inside chunks.** These are excluded by instruction. This includes
  the British 'colour' argument and 'grey50' in the ggplot code,
  which are code identifiers, not prose.
- **Verbatim quotations from Vossoughi et al. (2012).** These keep
  their original wording and British spelling ('analyse'), and keep
  their LaTeX double quotes. The quotations exceed 25 words, but a
  quotation cannot be split.
- **Bullet lists of adequacy conditions and performance metrics.**
  These are already short enumerations.
- **Substantive claims flagged in the referee report.** These were
  left with their meaning unchanged, even where the prose was split.
  Examples are 'slightly conservative ... throughout', 'MMRM estimand
  ... conditional on the covariance model', 'exact p-values',
  'convergence failures are common', the Fitzmaurice 'two-stage
  random-effects lineage' sentence, and 'They' for the
  single-authored Dawson (1994) paper. Fixing them is a separate
  decision for the author.
- **Misattributed citations** (`ardMixedModelRepeated2011`,
  `ardSimulationStudyComparing2015`,
  `liuMixedModelsCovariateInteractions2021`). These were untouched,
  per the instruction not to change cite keys or claims.
- **Remediation-history prose** (the 'whitepaper M1' references and
  the Appendix audit narrative). This was revised for sentence length
  only and not removed. The referee report recommends removing it,
  which is a content decision.
- **The '30 non-null-missingness conditions' sentence** (about 36
  words, mostly math). It was left as is, because splitting it would
  separate the arithmetic from its conclusion.
- **'offers a compelling alternative'** (Introduction). This is
  borderline promotional, but softening it would change the strength
  of the claim, so it was left.
