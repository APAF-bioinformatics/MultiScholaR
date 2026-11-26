# KIT-272 outcome: report text for panels d to f of Figure 1

Branch `kit-272` from origin/v0.9-staging 7f17547f9. Text only: no R code, no figure, no file or tab name changed.
All sentences are DRAFTS for the owner's wording pass (VAL-332); the tests hold them in file-local constants.

## What was run (observed)

| Run | Commit | Purpose | Result |
|---|---|---|---|
| 45726305082 (R1) | 84ee54db | new file only, old templates | 116 passed, 281 failed, 0 errors. Every predicted row failed by value at the predicted count; 12 unpredicted failures in the guards T8 and G1' (variant `unrecorded`): instrument defect 1 |
| 45728034972 (R1b) | 57d83190 | R1 after the guard fix | 104 passed, 269 failed, 0 errors: exactly the prediction. T8, G1', T9 pass; T1 to T7, T10, T11 fail by value |
| 45729807390 (R2) | 085e57e5 vs 57d83190 | fix against tests-only, 19 files, check | all files as predicted; new file 369 passed 4 failed (T7 framing: instrument defect 2); R CMD check identical (0 ERRORs, 3 WARNINGs, 3 NOTEs both sides) |
| 45735465246 (R2b) | 6e2b61ab vs 57d83190 | new file only after the T7 fix | 373 passed, 0 failed, 0 errors on head (base 104/269): as predicted, T7 passes |

Core-hours used: 1.37 + 1.33 + 7.48 + 1.47 = 11.65 (job-minutes x 4 / 60) against expected 12.06 and ceiling 60.33. Instrument defects: 2 of 3 allowed (the guards against the frozen base; the T7 framing), no result re-run.

Standing failures, identical base and head: test-kit220-da-threshold-consumers.R 2, test-omic-summary-module.R 28,
test-omics-summary-core.R 25 (not in the ticket's list of standing failures for this scope; unchanged by this ticket).

## Deviations from the ticket text (what is true)

1. **The declared edits of KIT-227 G1, KIT-236 G1 and KIT-235 T6 are in the fix commit, not the tests-only commit.**
   Ticket A1 asks that G1' "passes on the old templates". An edited G1 substitutes the NEW text on the base side, so
   against the OLD head templates its equality cannot pass; and R2's prediction (`same` on both sides) holds only if
   the edits arrive with the fix. The new file's G1' is a guard instead: the frozen base copies hold each old element
   exactly once (the match assertions of the declared edits), passing on both sides.
2. **Guards T8 and G1' cover the three-stage and two-stage variants only.** The frozen KIT-227 base copies predate the
   "RUV not recorded" wording (NA prints their two-stage text), so the unrecorded variant has no base counterpart in
   them. The unrecorded variant is pinned by T1, T5, T10.
3. **T6 for the limpa template** is exercised on both routes (the route paragraph through `record_present` and
   `dpc_engine` set on the numbers, and the other paragraph): five printed Methods phrases in all.
4. **R2 file count is 19 on staging** (the ticket's 19 included KIT-266's file, which does not exist; the grep also finds
   test-kit264-trx-export-current.R and the new file).
5. **The lipidomics and metabolomics Figure 1 caption** (`lipidomics_report.rmd:796`, `metabolomics_report.rmd:783`,
   "density boxplots of the samples' scores on the first two principal components") also says "density". It is not on
   the ticket's list (section 2.2 and section 5 item 5), so it is untouched; recorded for the coordinator.

## Acceptance criteria

| # | Evidence |
|---|---|
| A1 | R1b: T1 to T7, T10, T11 fail by value at the stated values; T8, T9, G1' pass; deviation 1 for the declared edits |
| A2 | R2 / R2b: every row passes on the fix commit (R2b for T7) |
| A3 | R2: every other file `same` (the three declared-edit files gain only the substitution assertions, 0 failed on both sides); check identical |
| A4 | the sentences are those of ticket section 5 (VAL-332), asserted through the draft constants |
| A5 | `git diff` of the templates shows only the lines of section 2.2, the PCA bullets and the four lipid/metab bullets |
| A6 | NEWS, ledger row `VAL-<merge>`, this file |

## Not checked

- No R was run locally (runner only). A rendered report was not looked at: the tests call `FormatResults`,
  `FormatMethodDetails*` and the placeholder chunk, not a knit.
- Which wording the owner will finally choose; the notch statement rests on the audit's reading of ggplot2 4.0.3.
- Merge with KIT-266 (no branch); its R3 reference and test file.
- Windows. Real data.

## For the change report

**What was wrong.** Since KIT-055 (new to 0.9; `main` is right) the four proteomics report templates said that panels d
to f of Figure 1 (c to d with two stages) were "Density plots: the distribution of protein intensities per sample", and
the Methods paragraph called them "density boxplots". They are notched box plots, per colour group, of the samples'
scores on principal components 1 and 2 (`buildProtNormDensityPlot`, `R/utils_workflow_dia_norm.R:229-249`;
`plotPcaBox`, `R/func_prot_s4_qc_plot_methods.R:311-398`). Every proteomics report of a 0.9 build was affected; no number
was wrong, the description of a figure was. Lipidomics and metabolomics reports said "Density boxplots to assess
distribution consistency" for box plots that come from a separate PCA.

**Code references.** Base 7f17547f9: `inst/reports/proteomics/{DIANN,DIANN_limpa,LFQ,TMT}_report.rmd` Results bullets
(DIANN 703, 718, 734; limpa 595, 612, 628; LFQ 655, 670, 686; TMT 646, 661, 677), Methods phrase (DIANN 677; limpa 506,
568; LFQ 627; TMT 618), placeholder label (DIANN 1562; limpa 1500; LFQ 1515; TMT 1506); lipidomics 762, 812;
metabolomics 749, 799. Fix commit 085e57e5: the same places (lines shift by the PCA bullets' collapse to one line).
`main` carried the right text at `main:inst/reports/proteomics/DIANN_report.rmd:675-677`.

**What changed.** Text only. The bullet now names the scores, the group, the box, whiskers and notch (McGill's rule of
thumb for comparing medians, not readable below about ten samples per group), and says each stage is a separate
analysis; the PCA bullet states the feature subset; the Methods phrase is "box plots of the principal component
scores"; the placeholder label is "principal component scores by group"; lipidomics and metabolomics say their box plots
use every lipid or metabolite with at least two values, separately from the PCA panels. Projects that already generated
a report keep their template copy (VAL-202's stated limit; NEWS says to delete the copy).

**Evidence.** Runs 45726305082, 45728034972, 45729807390, 45735465246 (see above and the records in this folder).
Limits: by reading plus text tests; no render.
