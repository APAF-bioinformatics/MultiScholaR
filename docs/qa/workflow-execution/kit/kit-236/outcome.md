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

## Owed after KIT-242 merges (left exactly as they are; tests W1 to W4 record today's text)

- **D8, the between-sample normalisation sentence.** Lipid/metab "Between-sample normalization was performed using the '{norm_method}' method from the 'limma' R package …". This includes its "cyclicloess" default when nothing is recorded (W1, W2). Proteomics B2 is the same: "'none' method from the 'limma' R package" for "none" (W4).
- **D19, the internal-standard sentence.** "Internal standard (ITSD) normalization was applied using the average_centered method with median aggregation …". It still prints when no standard was found (W3).
- Source for D19, from G0 V6: KIT-235's parser returns a written `* ITSD Features Per Assay:` block as its own section, `study_params[["ITSD Features Per Assay"]]`. This ticket's Figure 1 caption already reads that section (R1-N07).
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
