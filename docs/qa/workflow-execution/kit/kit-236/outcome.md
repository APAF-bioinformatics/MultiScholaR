# KIT-236 outcome: the reports print only what the analysis computed, ran or wrote

**State:** done, ready for its implementation audit. Two sentences are owed after KIT-242 merges (listed below).

- Branch `kit-236`, from staging 29065062.
- The gate G0 ran on the scripts-only branch `kit-236-g0`, which is never merged. Its records are copied here, under `g0-45577606908/`.
- Ticket: `tickets/kit/KIT-236.md` r1.3. The owner's OD-1 is VAL-245: "take all recommendations."

## Commits

| Commit | What |
|---|---|
| 1b829e04 | Tests first: `tests/testthat/test-kit236-report-sentences.R` and the base templates `tests/testdata/kit236/base/` (the six templates verbatim at 29065062). No package code. |
| da192fd4 | Tests: instrument fix 1 after R1 (VAL-54). No expected value changed. |
| 6a45d288 | The fix: six templates under `inst/reports/`, the declared test edits (a) and (b), NEWS, the ledger row `VAL-<merge>`, `wording.md`. |
| (records) | This file, `runs.tsv`, and the copied result tables. |

## Runs (`runs.tsv`)

**G0, run 45577606908**, on kit-236-g0 10198a0d = KIT-235 dc2c638f plus the G0 script.
- Every section-4 prediction held. G0a observed 54 lipid and 55 metab sentence rows on the GUI lanes.
- G0b covered V0 to V7 through `RenderReport`; G0c the chunk cases and part B. Details are in `g0-45577606908/outcome.md`.

**R1, run 45606119057**, on tests-only 1b829e04: an instrument defect.
- In the test script, the vector of proteomics template names (`k236Prot`) was overwritten by the runner function of the same name. B1, B3 and W4 therefore stopped with "invalid for() loop sequence".
- G1's log2 check compared head with base minus the zeros sentence, which encodes the change, so it cannot pass on base.
- Fix (da192fd4): the vector is renamed `k236ProtNames`, and G1 compares both sides without the zeros sentence.
- Conditions (i) to (iv) of VAL-54 hold. No expected value changed. This is instrument rerun 1 of 3.

**R1 rerun 1, run 45608054910**, on tests-only da192fd4: the valid failing-first run.
- 63 passed, 174 failed, 0 errors.
- Every non-guard test fails on the base templates: T1, T1c, T2, T3, T3i, T4, T5/T6n/T6r, T7/T7b, T8, T9/T9b, T10/T10b, T11/T11b/T11g, T12/T12g, T13, T14, T15, T17/T17b/T17g/T17n, T18/T18g, T19/T19g, B1, B3, B4, B5, B6.
- The guards pass on base: G1, W1 to W3 and W4.
- Stop rule 1 is not triggered.

**R2, run 45609991452**: fix 6a45d288 against tests-only da192fd4, 34 files, 3 shards per side, with check. Result tables are in `r2-45609991452/`.
- `test-kit236-report-sentences.R` is `fixed`: head 237 passed, 0 failed.
- `test-kit227-report-ruv-sentences.R` is `count-changed`, both sides 0 failed. G1 went from 46 to 42 passed expectations: the declared edit (a) removed the four lipid/metab `results` assertions.
- `test-kit235-lipid-metab-report-parser.R` is `count-changed`, both sides 0 failed. T6 went from 16 to 20: the declared edit (b) compares two paragraphs where it compared one text.
- Every other file (31) is `same`. That includes the three GUI-entry files. Their existing failure lists are identical on both sides: proteomics 4, `test-omics-summary-core.R` 25, `test-omic-summary-module.R` 28, `test-kit220-da-threshold-consumers.R` 2. Those are runner-image failures of stored references, present on base.
- R CMD check: findings identical on both sides (6). Nothing only on head, nothing only on base.
- Stop rules 2 and 3 are not triggered.

**Core-hours:** 312.1 job-minutes in all (G0 10.8, R1 10.5, R1 rerun 11.2, R2 179.6), times 4 / 60 = **20.8 core-hours**. The ticket expected 17.0, with a ceiling of 278.67. No own R3, as the coordinator instructed: R3 is the shared full suite.

## Acceptance criteria

- **G0 observed as predicted:** met (run 45577606908).
- **Section 5's rows fail on base and pass on head; guards pass on both:** met (R1 rerun 45608054910; R2 45609991452).
- **R2 `same` on every file but the new one and the two declared edits; the declared files differ only by those edits; R CMD check unchanged:** met (R2 `compare.tsv`, `check-compare.txt`).
- **Every other sentence of the lipid/metab reports unchanged (G1):** met for the Methods paragraphs G1 guards:
  - the RUV paragraph;
  - the quality-control bullets;
  - the ITSD paragraph;
  - the between-sample normalisation paragraph;
  - the assay sentence;
  - KIT-172's log2 sentence.

  R2 shows this for both lipidomics and metabolomics, with and without a record.

## Owed after KIT-242 merges (left exactly as they are; tests W1 and W3 record today's text)

Correction pass (I12, the coordinator's ruling): D8 is split. Its "none" form ("No between-sample normalisation was applied.") and its "not recorded" form ("The between-sample normalisation method is not recorded for this analysis.") are released now in lipidomics and metabolomics (T6, T6b), and the proteomics "none" form in the four proteomics templates (B2). What still waits: the lipid/metab sentence of a method applied ("Between-sample normalization was performed using the '{method}' method …", W1), the proteomics sentences of a method applied and the proteomics not-recorded form (still "'cyclicloess'": B2 was defined as the "none" sentence, and KIT-227's G1 guards the proteomics default text of exactly this not-recorded case), the internal-standard sentence D19 (W3), and the method named in the Figure 1 caption ("after {method words}"), which shares D8's KIT-242 dependency (I12). The earlier list follows.

- **D8, the between-sample normalisation sentence.** Lipid/metab "Between-sample normalization was performed using the '{norm_method}' method from the 'limma' R package …". This includes its "cyclicloess" default when nothing is recorded (W1, W2). Proteomics B2 is the same: "'none' method from the 'limma' R package" for "none" (W4).
- **D19, the internal-standard sentence.** "Internal standard (ITSD) normalization was applied using the average_centered method with median aggregation …". It still prints when no standard was found (W3).
- Source for D19, from G0 V6: KIT-235's parser returns a written `* ITSD Features Per Assay:` block as its own section, `study_params[["ITSD Features Per Assay"]]`. This ticket's Figure 1 caption already reads that section (R1-N07).
- The Figure 1 caption's `{method words}` (I12): they name the chosen method for an assay whose normalisation silently fell back, as D8 does.
- At the second merge of KIT-236 and KIT-242, T10b's end-to-end case runs: this template's D11 reader against a record KIT-242's code wrote (R1-N09).

## Observed and inferred

**Observed**
- The G0 lanes printed all 103 lipid/metab sentences predicted. Verbatim texts are in `g0-45577606908/g0a/*-sentences.tsv`.
- R2 shows the new texts on the real code path: Main, AnalyzeSamples and the record overlay for lipid/metab; AnalyzeSamples and the Format* functions for proteomics.

**Inferred, not run**
- How a full render looks with the new chunk `eval=` conditions. The GUI-entry lanes rendered the reports on head (R2, `test-gui-entry-lipidomics.R` and `-metabolomics.R` `same`; the docx is recorded as present), but the docx text of head was not read in R2.

## Deviations from the ticket, stated

- **Declared edit (a) is extended to DIA with limpa** (r1.3 section 15, from G0c). In KIT-227 G1, the DIANN_limpa method text is compared with the base text minus the averaging paragraph.
- **Template code outside the listed functions.** These are needed by R1-N05 (the record read in the template) and by the scope addition F6:
  - one line in each AnalyzeSamples (`report_facts`);
  - the helper functions in each Sample Analysis Functions chunk;
  - FormatDESummaryTable and Main in the lipid/metab templates;
  - the five figure chunks of each lipid/metab template.

  No `R/` file is changed. Not touched: the parser, the RUV paragraph, the log2 sentence, the enrichment text and chunks.
- **ASCII hyphen.** Benjamini-Hochberg is written with the ASCII hyphen the templates already use, not the ticket's en dash (`wording.md`).
- **Ledger row cost.** The proteomics overlay's "Parameter not recorded; default shown" note now shows `NA` instead of `0` for the TREAT cutoff, when a record has no DA `lfc` and the file has no cutoff. An `NA` cutoff is how "not recorded" reaches the TREAT sentence (R1-N02). Such records are not written by the application's DA step.

## Not checked

- The rendered docx text of head on the GUI lanes. R2 ran the lanes, but no G0-style text check was run on head.
- The parse shape of a written `ITSD Features Per Assay` block, beyond its being a section with at least one entry (T3i).
- Scripted (non-GUI) lipid/metab routes.
- The proteomics records' design rows: B1 is tested from the design table; the record form is the same code as T11.
- The wording itself. Every new sentence is a draft for the owner's pass (`wording.md`).

## For the change report

**What was wrong and who was affected.** The lipidomics and metabolomics reports printed fixed conclusions for every data set, whatever the data:
- the Figure 1 captions ("batches merged effectively", "correlation … was high and consistent");
- Comments ("All samples passed initial QC criteria");
- Opinions ("data quality is sufficient for differential analysis").

Their Methods described things that were not done or not as done:
- an integration step and a cross-assay check that do not exist;
- zeros "coded as NA" whatever the transformation;
- "each pair of experimental groups" instead of the declared contrasts;
- q-values without the Benjamini-Hochberg fallback;
- "the standard empirical Bayes method … zero" when TREAT ran. The template read the DA settings under a section name the writer never uses;
- placeholders such as "SOP xxx" and "Salesforce ID".

The sample details always said "no technical replicates". The report named a results file that is not written (`DE_lipids_results.xlsx`) and Figures 2 and 3 when they were absent, printing "Warning: … not found" lines. The count sentence was always empty, because the template read a file name Copy does not write.

The proteomics reports shared four kinds: "no technical replicates", the BH fallback, the TREAT reader and the missing block clause. The DIA-with-limpa report said technical replicates were averaged, which the application never does.

Who was affected:
- **Lipid/metab:** users of builds from 2026-01-09 (metabolomics) and 2026-01-18 (lipidomics) to 2026-03-06; users who re-rendered older `•`/`-` files; and every lipid/metab user once KIT-235 merged. The report failed with "Error running Main function" between 2026-03-06 and KIT-235.
- **Proteomics:** users who declared technical replicates, set a fold-change cutoff, hit the BH fallback, or used DIA with limpa without its record.

**Legacy on `main` (436cdbe8) and new to 0.9**
- **Legacy, on `main`:**
  - every lipid/metab sentence of the ticket's section 1 (section 7 lists main's lines);
  - the TREAT reader;
  - the empty count sentence: main's template reads `de_<omic>_num_sig_de_molecules.tab`, `inst/reports/lipidomics/lipidomics_report.rmd:1125`, while main's Copy writes `da_<omic>_num_sig_da_molecules.tab`, `R/func_general_filemgmt.R:2234`;
  - the proteomics "no technical replicates", Storey and TREAT sentences.
- **New to 0.9, not on `main`:** the DIA-with-limpa "Technical replicate averaging" paragraph. `git show main:inst/reports/proteomics/DIANN_limpa_report.rmd` has no such text.

**Code references (base 29065062 → fix 6a45d288)**

Lipidomics (metabolomics is the same functions in `metabolomics_report.rmd`):

| Function or chunk | Base | Fix |
|---|---|---|
| FormatSampleDetails | `:493` | `:496` |
| FormatMethodDetails | `:526` | `:513` |
| FormatResults | `:771` | `:764` |
| FormatComments | `:867` | `:810` |
| FormatOpinions | `:890` | `:834` |
| FormatDEDetails | `:926` | `:842` |
| ReportRecordFacts and helpers | none | `:969` |
| FormatDESummaryTable | `:1121` | `:1281` |
| Main | `:1160` | `:1332` |
| Figure chunks | `:1285-1342` | `:1461-1495` |

Proteomics:
- LFQ: FormatSampleDetails `:365`, FormatDEDetails `:714`, ReportRecordFacts `:972`.
- TMT and DIA-NN: the same functions.
- DIA with limpa: FormatMethodDetails `:352` (`replicate_text` `:363`).

**What changed in results.** No computed value. The changes are only to what the report says: the sentences above are now computed from the analysis record, the design and the files beside the report, or removed.

**Evidence.**
- G0 run 45577606908: each false sentence observed on the GUI lanes, with seven planted variants through the real render.
- R1 runs 45606119057 (instrument defect) and 45608054910 (failing-first).
- R2 run 45609991452 (fixed, `same`, check identical).

**Limits.**
- A project keeps the template copy made at its first report; older projects keep the earlier text unless `scripts/<omic>/<template>.rmd` is removed.
- D8 and D19 are unchanged until KIT-242.
- The lipid/metab FDR sentence states the fallback as a rule until KIT-242 records it per comparison.
- The volcano and count figures missing on the lipid/metab GUI route stay missing (KIT-220 H2, plan section 13); the report no longer refers to them.

## Correction pass (AUDIT-KIT-236-impl-20261125; the coordinator's rulings of 2026-11-26)

Branch `kit-236` merged with staging c4505c7f (8caf2d7d; ledger: staging's rows kept, this ticket's after the last; NEWS: both entries; no template or declared test file changed on staging since 29065062). Tests-only commit d06d424d; correction commit after it.

| Finding | What changed | Test (fails on the tests-only commit, passes on the correction) |
|---|---|---|
| I01 | D10 under the declared RUV design: "on the uncorrected abundance matrix, with each assay's RUV factors as covariates in the model ({assay}: {k} factors; …)" and "the formula '{formula}' plus each assay's RUV factors". Source: the record's `<prefix>_ruv` step (its `ruv_design`, parameter or outcome, and its `best_k` per assay). The record's DA `model` does not carry `tested_matrix` or `w_terms` for lipid/metab (READ: `acceptWorkflowOmicsDaManifest` keeps only the limma fields and the block record, `R/utils_workflow_omics_declarations.R:1765-1768`; the per-contrast model with them is not copied to the record), so the RUV step is the record's statement of the design | T20 |
| I02 | "k of n": n = k + the filter's own removals (the filter step's columns less the columns of the step before it); k = the record's analysed design. Wording unchanged; the redefinition of n is parked row 40 for the owner (marked draft in `wording.md`) | T1d |
| I05 | Unlabelled steps (Skip's `<prefix>_norm_complete`, a reset, `legacy`) are left out of the step list; `<prefix>_ruv` is "RUV correction" | T1e |
| I06 | "p-values were calculated using the 'eBayes' function" only with no fold-change threshold; TREAT: "the variances were moderated with limma's empirical Bayes method ('eBayes'), and p-values were calculated with limma's TREAT test against the fold-change threshold"; not recorded: says so. Lipid, metab, LFQ, TMT, DIA-NN, DIA with limpa (non-limpa branch, from the record's model only) | T21, B7 |
| I07 | At a fold-change threshold the count sentence and the volcano sentence read "(q-value < t and an absolute log2 fold change of at least c)" | T19c |
| I08 | The count sentence names the comparison as the DA paragraph does (the record's friendly name); "summarised" as KIT-220; T19's table in the application's shape ("groupKO-groupWT") | T19 |
| I09 | D6: "Each assay was processed and tested for differential abundance separately; the results of the {n} assays are reported together." ("removed from every assay" left out: false until KIT-255); with a failure the results are said to be reported "for the assays and comparisons where the test could be computed", and D10 names each failure with its error | T8, T22 |
| I10 | NEWS, the template comment and T17b's label: 'Differential Expression UI Parameters' is the only name any writer has used; the template read a name nothing writes | (records) |
| I11 | The ledger row in the five-column form, "KIT-236 implemented VAL-245 …", with what is owed and no overclaim, after VAL-261 | (records) |
| I12 | D8 split: the "none" and "not recorded" forms released (lipid/metab), the proteomics "none" form released (B2); the rest waits (above) | T6, T6b, B2 (W2 and W4 leave) |
| I03 | Recorded: **KIT-249 must merge before or with KIT-236.** On staging the proteomics BH fallback leaves the q-value column empty when `qvalue()` errors (`R/func_prot_da_model_stats.R:686`) while the record says "BH", so the proteomics FDR sentences of this ticket are true only with KIT-249's fix | — |
| I13 | Recorded as a correction: a DIA report (DIA-NN, DIA with limpa) without an analysis record now prints the q-value threshold the user set (the DA UI block the DIA route writes is read) instead of a fixed 0.05 | — |
| I16 | Change report: the missing volcano figure on the lipid/metab GUI route is new to 0.9 (`main`'s Copy copied `Volcano_Plots` for every omic); the sentence text is legacy | — |
| I17 | Routed to the owner, not this ticket: the proteomics volcano sentence "… and no log fold-change threshold", unconditional since 9c8c6830 (2025-04-08), reads against D20 at a fold-change threshold above 0 | — |

### Render predictions (I04; written before the run)

The correction run's scripts job (`render.R` in this directory) runs the six GUI lanes at the correction commit and reads each report's text; for lipidomics and metabolomics it renders the lane's project again under V1 to V7. Predicted (each phrase present or absent; the full table is in `render.R`):
- lipid_canonical: "The analysis included 6 samples, each a separate biological sample (no technical replicates were declared in the design)."; "Data analyses were performed with MultiScholaR 0.9.0 (build …); SOP reference: [to be completed by the analyst]."; D6 without failures; KIT-172's log2 sentence and no zeros sentence; the held ITSD sentence and the held 'cyclicloess' sentence; D9; D13; caption "… on the intensities before log2 transformation (left); after cyclic loess normalisation (normalizeCyclicLoess) (right). RUV-III correction was not applied to this dataset. …"; "Features remaining after each step (total across the 2 assays): import and design 20; intensity filter 17; duplicate resolution 17; quality-control finalisation 17; internal-standard normalisation 17; log2 transformation 17; between-sample normalisation 17; sample-correlation filter 17."; "6 of 6 samples were kept by the sample-correlation filter (threshold 0.5)."; the two analyst fields; D10 "(1 comparison: KO_vs_WT)" with "p-values were calculated using the 'eBayes' function"; the FDR rule; "tested for any difference from zero log2 fold change"; "q-value below 0.05"; the `DA_results_lipidomics.xlsx` sentence; no volcano sentence; the count sentence "summarised in Figure 2 (q-value < 0.05). For the comparison KO_vs_WT in LCMS_Pos, there were …" and the Figure 2 caption; no Figure 3 caption; no "not found" warning of the four chunks this ticket owns.
- metab_lc: the same with 6 features at every step, no count sentence and no Figure 2 caption.
- prot_lfq, prot_tmt, prot_dia, prot_dia_limpa: no "no technical replicates."; no "The design is not available"; "False discovery rates were estimated as q-values (Storey, 2002)"; no "standard empirical Bayes … zero log-fold-change" sentence; no "Technical replicate averaging"; no "'none' method"; LFQ, TMT and DIA-NN "tested for any difference from zero log2 fold change" (the limpa route prints no TREAT sentence).
- V1 (none): "No between-sample normalisation was applied." and the caption's "after the normalisation step, at which no normalisation was applied (right)". V2 (offset 1): the transform and the zeros sentence, no "coded as NA". V3 (record technical replicates): "The analysis included 6 measured runs from 5 biological samples; 1 biological sample was measured more than once …". V4 (images): the volcano sentence and the Figure 3 and Figure 2 captions. V5 (no xlsx): no results-table sentence. V6 (ITSD block): the caption "before internal-standard normalisation and log2 transformation (left)" and the held ITSD sentence. V7 (record lfc 0.5 and `treat_applied` TRUE): the TREAT sentence, the TREAT p-value clause, the significance rule with 0.5, no "'eBayes' function" p-value clause.
