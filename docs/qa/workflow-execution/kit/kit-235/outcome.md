# KIT-235 outcome: the lipidomics and metabolomics reports read the study-parameters file the application writes

Ticket `tickets/kit/KIT-235.md` r1.2. Branch `kit-235`. All computation ran on GitHub-hosted runners; nothing was run on
the laptop.

## Commits

| Commit | What |
|---|---|
| 28916e2d | merge origin/kit-227 5e7efe5f into staging 67569959 (the ticket's base, "staging after KIT-227's merge"; KIT-220 touches neither lipid/metab template nor the writer) |
| 96e3c813 | records: G0 |
| 80fa4c48 | tests only: `tests/testthat/test-kit235-lipid-metab-report-parser.R` (T1 to T6, branch A), `tests/testdata/kit235/base/{lipidomics,metabolomics}_report.rmd` (the ten chunks at the branch point, verbatim); no package code |
| f53db55a | merge v0.9-staging 3eeaf100 (R2's base) |
| dc2c638f | the fix (F1) and NEWS |

G0 ran from a separate branch `kit-235-g0` (staging 67569959 + `tools/kit235/g0.R` only; not for merging): 90802bdf,
then c1c03299 (instrument fix 1).

## Runs (`runs.tsv`)

| Run | Commit | Purpose | Job minutes | Result |
|---|---|---|---|---|
| 45525539989 | 90802bdf | G0 (a, b, c) | 14 | G0a observed; G0b/G0c instrument defect (below) |
| 45573071526 | c1c03299 | G0 instrument rerun 1 (VAL-54) | 12 | G0b, G0c observed |
| 45575079261 | 80fa4c48 | R1, tests only | 12 | T1 to T5 fail at the stated values |
| 45577137662 | dc2c638f vs f53db55a | R2, 15 files, 3 shards per side, check | 122 | new file fixed, 14 files same, check identical |

Total 160 job minutes = 10.7 core-hours (job minutes x 4 / 60, set-up included), against 13.0 expected and 278.67
ceiling. No correction pass, one instrument rerun.

## Observed

### G0a (run 45525539989; repeated identically in 45573071526): `g0-45525539989/g0a/`

The writer's own files (`createWorkflowArgsFromConfig` on the fixtures of
`test-general-filemgmt-s4-contracts.R`: M metab per-assay RUV, L1 lipid single form, L2 lipid per-assay, S lipid and
metab without RUV) went through the templates of main 436cdbe8, staging 67569959 and kit-227 5e7efe5f. The chunks were
evaluated in `new.env(parent = asNamespace("MultiScholaR"))`, as `RenderReport` does.

- With the writer's `*` bullets, every case and every version gives the Sample Details "Error running Main function:
  missing value where TRUE/FALSE needed", with method_details empty (0 characters).
- The extractor called alone returned 41 fields; AnalyzeSamples returned 46 numbers and swallowed no error; no chunk
  errored. The stop is the template's, not the harness's.
- The same files with `•` bullets print "METHOD DETAILS" on all three versions.
- Fields: as section 1.3 of the ticket predicts. For example, staging M `*` gives num_assays NA,
  itsd_applied FALSE, normalisation_method "cyclicloess", ruv_assay_list `"* Status"`, `"* Grouping Variable"`. kit-227 M `•`
  gives "2", "assay_pos, assay_neg", TRUE, "median", "sum", TRUE, "1", "PQN", "group", with assays k 2/3 and
  controls 1/2.
- Preconditions: pandoc 2.9.2.1 available, locale C.UTF-8, rmarkdown 2.31, knitr 1.51.

### G0b, the application's real report (run 45573071526): `g0-45573071526/g0b/`

The GUI-entry lanes `lipid_canonical` and `metab_lc` ran the real Summary save and Generate Report (20 captures each).
The Word reports' text (`lip-report-text.txt`, `met-report-text.txt`, from `word/document.xml`) holds, in place of
the Sample Details:

> Error running Main function: missing value where TRUE/FALSE needed

and no "METHOD DETAILS". The captured `study_parameters.txt` (`*-c17-src-study_parameters.txt`) has 23 lines starting
`* ` and none starting `•` or `-`.

The first run's capture failed: the trace could not find `kit235Capture` from the helper's frame, so 0 files were
captured and nothing was judged. This was fixed by embedding the function in the tracer call. No expected value, rule
or ceiling changed.

### G0c, the GUI files through the fixed reading: `g0-45573071526/g0c/checks.txt`

Each lane's file was rewritten to `•` and run through the staging and kit-227 chunks with its own record and design.
The values were the same on both templates and both lanes:
- num_assays "2";
- assay_names "LCMS_Pos, LCMS_Neg";
- log_transform_applied TRUE;
- normalisation_method "cyclicloess";
- the RUV paragraph "RUV-III batch correction was not applied".

All as predicted.

- **ITSD, as predicted (decision rule item 1):** neither file has an "ITSD Features Per Assay" block, and the
  Methods say "Internal standard (ITSD) normalization was applied using the average_centered method with median
  aggregation …". This sentence is false: no standard was found on either lane. By the coordinator's routing it goes
  to KIT-236's sentence list, and the recording fault to a separate `R/` ticket. It is not changed here and is not a
  stop.
- **RUV form:** "Status: Not Applied"; no single form, no per-assay form, no Automatic section. G0c did not stop.
  **Branch A** is selected (E-6): F3 and F4 are withdrawn.

### R1 (run 45575079261, tests only): `r1-45575079261/`

29 failures, every one at the ticket's stated base value:
- T1 and T4: "Error running Main" present, no "METHOD DETAILS".
- T2: names `"* Status"`, `"* Grouping Variable"`.
- T3: num_assays, assay_names, itsd_method_type, itsd_aggregation and ruv_grouping_variable are logical NA;
  itsd_applied and itsd_removed are FALSE; offset is 0 (double); normalisation "cyclicloess".
- T3b: aggregation NA.
- T5: no "was not applied".

No T row other than T6 passed on base, so stop rule 1 did not fire.

T6's head-against-base comparison on the `•` and `-` shapes passed on base. Its second clause ("identical to head on
the `*` file") failed on base, 4 times (M and L1 by 2 shapes). It cannot pass there, because base stops on the `*`
file. So the ticket's "T6 passes on base by construction" is true only of its first clause. This is recorded for the
ticket text; the test was not changed.

### R2 (run 45577137662): `r2-45577137662/`

Base f53db55a = staging 3eeaf100 + kit-227 5e7efe5f + this ticket's tests and records. Head dc2c638f = base + F1 + NEWS.

- `test-kit235-lipid-metab-report-parser.R`: base 21 passed / 29 failed, head 50 / 0. Verdict `fixed`; T6 passes on head.
- The other 14 files are `same`:
  - `test-gui-entry-lipidomics.R` 35/0 and `test-gui-entry-metabolomics.R` 105/0 on both sides;
  - `test-kit227-report-ruv-sentences.R` 166/0;
  - `test-general-filemgmt-s4-contracts.R` 30/0 (the writer's pinned `*` shape);
  - `test-omic-summary-module.R` 1675, `test-omic-gui-config.R` 1911, `test-omics-summary-core.R` 676;
  - the remaining files as listed in `compare.tsv`.
- R CMD check: 0 errors, 3 warnings, 3 notes on both sides, findings identical (0 only on head).

## Acceptance

| Criterion | State | Evidence |
|---|---|---|
| G0 observed as predicted | met | G0a, G0b, G0c above; ITSD as predicted and routed; branch A |
| Section 5's rows other than T6 fail on base | met | R1, 29 failures at the stated values |
| They pass on head | met | R2 head 50/0 |
| T6 passes on both | partly: passes on head; on base its head-vs-base clause passes, its `*` clause cannot | R1, R2 |
| R2 `same` on every file but the new one | met | `compare.tsv` (14 same) |
| R CMD check unchanged | met | `check-compare.txt` |
| Stored references and existing expectations unmoved (stop rule 3) | met | R2: every existing file same, no file re-recorded |
| Diff only section 4's lines (stop rule 4) | met | `git diff f53db55a dc2c638f -- inst R`: 8 lines in the two templates (lipid :123, :172, :173, :187; metab :120, :169, :170, :184), no `R/` file |

## Inferred (not observed)

- On the GUI route the fixed report prints the Methods with the values G0c shows. G0c applied the `•` rewrite with the
  templates of staging and kit-227. F1 makes the `*` file parse identically to the `•` file: T6 checks this on head for M and L1.
- Proteomics and transcriptomics templates are untouched, so their reports are unchanged.

## Stated limits

- A project keeps the template copied at its first Generate Report. Projects that generated a lipid/metab report
  before this version keep the broken parser unless the copy is removed (KIT-227 OD-3, VAL-202).
- A report generated after a save that failed (the short file the Summary writes on an error) still shows "Error
  running Main function: …"; that file records no assays.
- (Branch A) A study-parameters file in the older single RUV form states no k and no control count. The application
  does not write that form: by reading on the applied route, and by G0c on the skip route.
- The newly visible sentences (ticket section 8) belong to KIT-236, KIT-219, KIT-227's wording pass and the wording
  pass. That includes G0c's false ITSD sentence on runs where no standard was found.

## Not checked

- The applied-RUV GUI route: every lipid/metab GUI-entry lane skips RUV. That the applied route writes the per-assay
  form is by reading.
- A real render of the fixed template on the GUI route. G0b observed the base only; head was checked through the chunks.
- Whether tags v0.5.0-rc1 or v0.9.0 were distributed.
- Hand-edited files with other bullet characters.
- The R3 full suite (shared, coordinator's E-2).

## For the change report

**What was wrong.** Lipidomics and metabolomics reports showed no Methods; this is a legacy defect. Commit 0281da66
(2026-03-06, "fix: remove UTF-8 BOM and non-ASCII characters from R source files") made the study-parameters writer
write `*` bullets in place of `•`. That was `R/func_general_filemgmt.R` on main; on staging it is
`R/func_general_filemgmt_study_parameter_writer.R`. The two report templates' parser read only `•` and `-`. So
"Number of Assays" was read as missing, the Methods code stopped at `if (num_assays > 1)`, and the report's handler
printed "Error running Main function: missing value where TRUE/FALSE needed" as the Sample Details. Methods, Results,
Comments, Opinions and the differential-abundance text were empty, while the application's status line said "Report
generated [OK]". The assay count, internal-standard settings, log transformation, normalisation method, RUV grouping
and per-assay k and controls were not read.

**Who was affected.** Every lipidomics and metabolomics user of the application since 2026-03-06, and scripted users of
the writer and the templates since then. This holds on main and in tags v0.5.0-rc1 and v0.9.0. Proteomics reports had
the same `*` loss for their RUV values only; that was fixed on staging in f91741f6.

**Code.**
- Base (main 436cdbe8, and staging at f53db55a):
  - `inst/reports/lipidomics/lipidomics_report.rmd:123`, `:172-173`, `:187`;
  - `inst/reports/metabolomics/metabolomics_report.rmd:120`, `:169-170`, `:184`.
- Fix dc2c638f: the same lines now read `*` as a bullet (`^[-•*]`, `^[•*-]`, `^[•*[:space:]-]+`), as the proteomics
  parsers have since f91741f6.

**What changed in results.** No computed value changed. The report's text changed: a lipid/metab report made from the
application's file now prints its Methods and the other sections. Files in the older `•` and `-` shapes print exactly
what they printed before (T6, head against the base templates' chunks).

**Evidence.**
- G0 runs 45525539989 and 45573071526: the defect observed on the writer's files under main, staging and kit-227, and
  on both GUI lanes' real Word reports.
- R1 45575079261 and R2 45577137662: T1 to T6 fail first and then pass; 14 neighbouring test files are unchanged; R CMD check is unchanged.
- Limits as stated above. The GUI evidence covers the RUV-skipped route only, with fixture data that hold no internal
  standard.
