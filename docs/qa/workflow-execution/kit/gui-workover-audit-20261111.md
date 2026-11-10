# GUI workover series: post-hoc audit findings (2026-11-11)

Audit of KIT-101 to KIT-113 on `performance-v2`, requested by the user against `gui-workover-audit-handover.md`. Thirteen independent read-only auditors: one per ticket, one for the series-level items. Audited at `b8300748` to `ec7db20e` (the lead's session committed KIT-113's run 9 during the audit).

**Method and limits.** Audit only: no repairs, no R, no tests, no runs, no KC data. Findings come from reading code, git history, committed TSV/JSON and the run traces. Every test count, lane result, e2e result and `R CMD check` result is taken from the lead's committed records, not re-run. Findings marked "by reading" were not reproduced.

## 1. Summary

**What holds.** No statistical kernel changed during the series (no commits to the normalisation, RUV, DA, NIPALS, imputation or limpa `func_` files after KIT-101). The shared modules pass the same arguments to each kernel as the deleted per-omic modules. Every oracle re-freeze is attributed and changed only what was claimed. The skip-baseline overlay, including KIT-113's rebuild, reproduces exactly. Every reported memory figure re-derives from the traces. `R/mod_*` fell from 56,123 to 22,320 lines; 27 of 28 omic × stage cells are on the shared frame.

**What does not hold, most serious first.**

| # | Finding | Ticket | Severity |
|---|---|---|---|
| 1 | Managed DA publication leaves provenance pending; a later DA re-click overwrites the previous step's real provenance entry (F112-05). Proteomics only. One-line fix; no oracle depends on the wrong value. Deferral to KIT-109 is not sound. | 112 | high |
| 2 | When `qvalue()` falls back to BH, the `fdr_qvalue` column holds BH values and nothing durable in the analysis record or exports says so; on the managed path there is no notice either (F112-01). | 112 | high |
| 3 | A manual internal-standard selection can resolve to different features in process than on the worker (P-19), by reading. | 105 | high |
| 4 | Managed path: a correlation click queued behind a running normalisation runs on the previous normalisation's result and commits it as head, by reading. | 111 | high |
| 5 | P-43's "+0.66 GiB" holds in one of three probe cases; the persistent cost (about 180 MiB per contrast kept by the app) is under-reported; the cheaper design (no plot objects in the results) is not offered. | 113 | high |
| 6 | P-05 describes the readers wrongly: in a managed GUI session all five proteomics readers turn non-numeric quantity text into silent NA. KIT-102 also changed PD-TMT behaviour without declaring it. | 102 | high |
| 7 | P-01's headline numbers (88 false positives; q-value against BH) were measured only on an unnormalised, group-confounded fixture, one seed. After F-18, RUV may refuse on real LFQ data (P-03 understated). | 101 | high |
| 8 | Lipid and metab keep 0 as a value (log2(0+1)) while their import validation counts 0 as missing. Same mechanism as F-18; recorded nowhere. | 101 | medium |
| 9 | No whole-suite run since KIT-105 (seven tickets, 35 `R/`-touching commits); 278 of 393 baseline test files not rerun. Every `check.log` is `--no-tests`. Check-mode blocks error, not skip; the fixture ships, root resolution (`here::here()`) fails. | series | high |
| 10 | 111m's memory floor was called a one-off; it recurs in 112m (2 of 4 members on builds with KIT-111 code, 0 of 9 before). P-34 was required by the ticket's rule and never raised. | 111/112 | medium |
| 11 | Enrichment provenance does not determine the result: recorded correction can differ from the one applied (clusterProfiler branch hard-codes BH), database version missing when there are no hits or the service failed, annotation table and DA files not identified. | 107 | medium |
| 12 | KIT-109's ticket names none of the items queued to it. | series | high |
| 13 | Ticked boxes that are not met: KIT-101 (three), KIT-103 (applet), KIT-104 (skip rule ticked on a table with failures), KIT-106 (two whole-state boxes), KIT-102 (`parity.md` never committed). | several | medium |

## 2. Findings by ticket

Verdict columns: implementation verified / evidence valid / claims supported.

### KIT-101 (design, diagnosis): partly / partly / partly

| Id | Sev | Evidence | Finding | Smallest repair |
|---|---|---|---|---|
| A101-01 | high | `kit-101/diagnosis/f19.txt` §2 | The 88 false positives (c') and q-value 41 against BH 23 were measured only on the unnormalised fixture (null logFC about -0.26). Centred variants give 19/18 and 22/20; c' was never run centred. One seed. | Reword P-01 and KIT-114 (a)/(b); re-measure on a fixed fixture over several seeds. |
| A101-02 | high | `fixes7/corpus-lfq-300-mnar.txt`; `func_prot_norm_optimization_helpers.R:983-993` | After F-18, RUV refuses when no control is complete in every sample. The corpus was built at 8% missing (18.9% complete); real LFQ is recorded at 51.8% missing. | Mark P-03 as blocking the real-data gate; count complete controls (aggregate only) before P-04. |
| A101-03 | medium | `func_lipid_s4_normalization_methods.R:583-586`; `func_lipid_import_detection.R:352`; `func_metab_import_detection.R:317` | Lipid and metab keep 0 as log2(0+1)=0 while import validation counts 0 as missing. Not in findings, parked decisions or KIT-114. | New parked decision (methodology, parity). |
| A101-04 | medium | `mod_prot_qc_protein_normalisation.R:44` | F-23 and KIT-114 (c) say E-level is exact for median-type methods; the GUI offers none, cyclicloess, quantile, scale. No exact method is offered. | Correct the text. |
| A101-05 | medium | `tickets/kit/KIT-101.md:245,255,259` | "No production file changed" ticked against about 20 `R/` commits. F-23, F-24, F-26, F-27, F-28 rest only on quotes in `findings.md`; F-28 has none. | Revision listing the authorities; user confirms the quotes. |
| A101-06 | medium | `member-results.txt`; `reverify.txt` | 101m/l/t recorded as rejected; "accepted and exact" is an offline re-verification re-run after a script edit; l/t judged against same-build oracles. | Say so in the outcome. |
| A101-07..11 | low/info | see auditor notes | Two dropped diagnosis questions; duplicated normaliser switch (`func_prot_s4_norm_methods.R:76-88`, `:108-112`); UniProt failure writes an empty table with no notice seen; mutation check never run on a proteomics lane; a stale comment. | KIT-109 or KIT-114. |

### KIT-102 (import): yes for structure, partly for parity / partly / partly

| Id | Sev | Evidence | Finding | Smallest repair |
|---|---|---|---|---|
| A102-01 | high | `utils_workflow_dia_sessions.R:1127-1134`; `utils_workflow_dia_state.R:1068,1106-1108` | Managed GUI: a character quantity column is converted at the design step with `suppressWarnings(as.numeric())`. P-05's "still fails at the design step" is wrong. | Rewrite P-05; one GUI-path test. |
| A102-02 | high | `func_prot_import_tmt.R:156,174`; `w3/test-first.txt` | PD-TMT, undeclared: text "0" was 0, now NA; a token in some columns now refuses with an opaque pivot error; a token in every column still imports silently. | Decide under P-05; name the column in the message. |
| A102-03 | medium | `func_prot_import_readers.R:82-84` | Spectronaut parity tested only with a mocked character column. | Record as defensive, or test through vroom. |
| A102-04 | medium | `git log --all -- kit-102/parity.md` empty | Ticked acceptance box and P-05 cite a file never committed. | Commit the survey. |
| A102-05 | medium | Q-25 to Q-30 | Old-against-new equivalence was a one-time run; about 1,040 expectations quarantined, none frozen. | Freeze the payloads (KIT-109). |
| A102-06..08 | low/info | | e2e file 63 to 61 expectations unlisted; `clean_acc` in `globalVariables` hides the check NOTE; `as.numeric` accepts hex text. | Notes. |

### KIT-103 (design): yes for the refactor, no for the applet / partly / partly

| Id | Sev | Evidence | Finding | Smallest repair |
|---|---|---|---|---|
| A103-01 | high | `tickets/kit/KIT-103.md:388` | Applet box ticked; outcome says unmet. | Untick, cite P-15. |
| A103-02 | high | `func_general_applets.R:236-241`; `utils_workflow_dia_state.R:1149-1153` | The applet's table has no `Run` column, so by reading it fails when the initial design is built. The test uses a made-up table. | Restate P-15. |
| A103-03 | medium | `test-omics-design-core.R:182-198`, `:484-510` | Two parity tests now compare the shared code with itself. | Freeze literals. |
| A103-04 | medium | `test-omics-design-core.R:787,864` | No live test pins lipid's refusal of an invalid imported design. | One live test. |
| A103-05 | medium | `utils_workflow_omics_design.R:276-287` | Metab import without `group`/`replicates` half-writes state, then fails. | P-11. |
| A103-06 | medium | `func_omics_design_shared.R:439-446`, `:708-712` | Pre-existing, unrecorded: reopening a design whose contrasts lack the `group` prefix and saving writes "A-B-". By reading. | Record; fix separately. |
| A103-07..11 | low/info | | Wrong explanation of a count; missing declared evidence files; P-08 row stale; proteomics `run_order` config never read. | Notes. |

### KIT-104 (QC): yes for P-16 and shell, partly otherwise / partly / partly

| Id | Sev | Evidence | Finding | Smallest repair |
|---|---|---|---|---|
| A104-01 | info | diff of `prot_dia_limpa.json` at `735ad485` | P-16 re-freeze verified: 32 captures in 14 steps, tab status and stepper only. | — |
| A104-04 | medium | `kit-104/skips-check.tsv` (`ca35c462`) | Skip-rule box ticked on a table showing 17 + 1 failures; post-fix proof is KIT-105's. | Note in outcome. |
| A104-05 | medium | same tables | Check-mode lost coverage starts here (`test-prot-qc-core.R`: 3 error blocks of 5). | Fixture root fix. |
| A104-07 | medium | `test-omics-qc-core.R:60-368` | 10 of 15 core blocks quarantined within the ticket; equivalence proof is one-time. | Freeze expected values (KIT-109). |
| A104-06, -08..12 | low | | `test-metab-01v` 7 to 3 expectations with no Q-row; metab parity narrower than stated; `is.numeric` sample-column discovery in new shared code (`mod_omic_qc_tabs.R:600,640-647`) against `R-shiny.md:316`, unrecorded. | Notes. |

### KIT-105 (normalisation): yes / partly / partly

| Id | Sev | Evidence | Finding | Smallest repair |
|---|---|---|---|---|
| A105-01 | high | `utils_workflow_omics_actions.R:410-425` against `:259-268`; `mod_omic_norm.R:275`; `mod_omic_norm_panels.R:444-461` | Manual ITSD: managed sends ids of the shown rows; in process rebuilds the table at Run time and applies the same row numbers. A second Run without Reset, or a QC refilter after the tab opened, can select different features. By reading. | Resolve from candidate ids on both paths. |
| A105-02 | medium | `kit-110..113/skips-*.tsv` | 278 of 393 baseline files not rerun after KIT-105. | Full suite, both modes. |
| A105-03 | medium | `kit-105/skipcmp.py` | Compares net counts per file, not block identity; ignores `file_error`, `static_blocks`; a targeted file not run is filtered out silently. | Assert the targeted set is present. |
| A105-04..09 | low/info | | P-18 row names readers that do not exist; removal of F-20's branch may show a spurious warning on a restored session (unconfirmed); trx `session_export` +7.9%; two functions over the parameter rule, unrecorded. | Notes. |

### KIT-110 (trx normalisation): yes / partly / partly

| Id | Sev | Evidence | Finding | Smallest repair |
|---|---|---|---|---|
| A110-01 | medium | `kit-110/findings.md` row 2; `test-trx-norm-module.R:219-220`; F111-01 | "Live payload" error on a managed retry after a refusal was filed as a harness note; KIT-111 fixed a real defect with that message. No test drives trx refusal then retry in one managed session. | One test; give it an id. |
| A110-02..07 | low/info | | Two existing test files grew unreported; `lane-changes.md` missing; validator accepts `run$handler` without `run$settings`. | Notes. |

### KIT-106 (DA tab): partly / yes for what it covers / partly

| Id | Sev | Evidence | Finding | Smallest repair |
|---|---|---|---|---|
| A106-04 | medium | scoping README item 7; plan `:37`, `:61` | "Interactive volcano off by default" is attributed to the user; no user quote exists, only "built on demand". | User confirms or reverses. |
| A106-01 | medium | `mod_omic_da.R:274-283`, `:457-461` | Two "no whole-state" boxes ticked but not literally true. | Reword. |
| A106-02 | medium | `correction/trx-e2e.md` | Managed trx volcano page never ran in a real browser. | Correct the outcome. |
| A106-03 | medium | `workflow_kit_driver.R:298-307` | F106-01 confirmed: proteomics' DA views rest on `testServer` only. | One e2e lane with assertions. |
| A106-05..11 | low/info | | DIA-specific second implementation inside the shared file; unread proteomics config fields; attribution proofs for 12 re-frozen values not committed; 106a `correlation` +19.9% unrecorded. | KIT-109. |

### KIT-107 (enrichment tab): yes / yes / partly

| Id | Sev | Evidence | Finding | Smallest repair |
|---|---|---|---|---|
| A107-01 | high | `utils_workflow_dia_enrichment.R:389`; `mod_omic_enrich.R:493-499` | Mixed-species import auto-ticks the organism filter, so every enrichment Run fails until it is unticked. | P-07 (a). |
| A107-03 | medium | `func_general_enrichment_orchestration.R:571-594` | clusterProfiler branch hard-codes BH; the record stores the selected method (default gSCS). | Record the effective correction. |
| A107-04 | medium | `utils_workflow_dia_da.R:~511-531`; `engine_helpers.R:182-185` | No hits and exhausted retries both return NULL; either commits "completed" with service `unavailable` and no database version. | Record the version independently; do not commit after a failed call. |
| A107-06 | medium | `runWorkflowEnrichment:846` | Annotation table not identified; the kernel reads DA TSVs from disk, not the recorded entities. | Content hashes in provenance. |
| A107-07 | medium | `replayWorkflowEnrichmentStep:1124-1132` | Replay compares nine parameters only, after the commit. | Compare the service resource. |
| A107-02, -05, -08..12 | low/info | | "Applied" flag set on a no-op; real version recording untested; failed commit swallowed; status text from current inputs; kernel: query strips `:.*` but background does not, organism names matched by substring. | With P-07/P-28; KIT-114 for the kernel. |

### KIT-108 (summary/report): partly / partly / partly

| Id | Sev | Evidence | Finding | Smallest repair |
|---|---|---|---|---|
| A108-01 | high | `src/workflow_memory.c:19-29`; `utils_workflow_saved_entities.R:1427,1435`; `utils_workflow_dia_state.R:2877` | P-30 prices route B as new native code; `reclaimWorkflowMemory()` already exists and is used. | Correct P-30. |
| A108-02 | medium | `DESCRIPTION:494` | Route D rejected because `readxl` "would be a new dependency"; it is imported. | Correct the diagnosis. |
| A108-03 | medium | `m4-replica.md:49` | After-report gain (0.60/0.36 GiB) is not isolated to route A; the replica shows the retained rise unchanged. | Reword. |
| A108-04 | medium | `func_general_filemgmt_results_summary_helpers.R:186,700,713,723` | F108-01: a failed workbook save is silent and the UI shows "[OK]". Covered by the ticket's own fix rule. | Fix with one test. |
| A108-05 | medium | `func_general_filemgmt_study_parameter_support.R:60-90` | F108-04: installed package falls back to `getwd()`, writing an unrelated repository's branch and commit into `study_parameters.txt`. | KIT-109. |
| A108-06 | medium | `skip-compare.md` | `test-omics-summary-core.R` fails 16 expectations in the main checkout. | Mask the version block. |
| A108-08..10 | low/info | | Workbook identity not pinned beyond sheet names; near-duplicate lipid/metab helper families moved. | KIT-109. |

Memory figures re-derived and equal to the lead's. The report step is not the workflow's peak (slice 2.82/2.90 GiB against 4.04/4.10 at normalisation).

### KIT-111 (proteomics normalisation on the worker): partly / partly / partly

| Id | Sev | Evidence | Finding | Smallest repair |
|---|---|---|---|---|
| A111-01 | high | `mod_prot_norm_workflow_helpers.R:621-643`, `:315-333`; `utils_workflow_execution.R:2933` | Source resolved at click time, then queued: a correlation click during a normalisation re-run commits the previous result as head. Managed only; by reading. | Resolve at start, or refuse while the runner is active; one test. |
| A111-02 | medium | stage JSONs; `kit-112/memory-outcome.md:49` | 111m's floor recurs in 112m; "one-off" unsupported; 112m was compared with the outlier (against 111n it is +0.096 GiB). P-34 not raised. | Restate; raise P-34. |
| A111-03 | medium | `utils_workflow_execution.R:748`, `:1717-1722`; `test-prot-norm-core.R:1619-1622` | Managed users see "worker_failed" without the kernel's reason; the test truncates notices at the first colon. Unrecorded. | Pass the error text through. |
| A111-04 | medium | `publishDIAReferenceRuvCovariates` | F111-02 confirmed, pre-existing, and wider: covariates then correction also compounds. | P-35. |
| A111-05 | medium | `test-prot-norm-core.R:397,420-422,473,504,1522` | Identity masks drop `config_list` (the RUV audit trail), most `@args`, dependency sets, record lineage; files by name only. | Add the audit fields to the capture. |
| A111-06 | medium | `helper-workflow-checkpoint-fixtures.R:11-12` | Check-mode: 79 of 1,553 expectations run; managed tests are load_all only. | `testthat::test_path` root (unverified). |
| A111-07..10 | low/info | | Dispatchers duplicate the in-process gate and resolution; 9 to 12 parameters; F111-05 real, display only; re-freezes verified. | KIT-109. |

### KIT-112 (proteomics DA publication on the worker): partly / partly / partly

| Id | Sev | Evidence | Finding | Smallest repair |
|---|---|---|---|---|
| A112-01 | high | `utils_workflow_execution.R:3002`, `:2630`; `utils_workflow_dia_state.R:830-838` | F112-05 confirmed: a later DA click with changed settings replaces the head state's real provenance entry with `da_publication`. Proteomics only. No recorded artefact found wrong (`.rds` oracles not opened). | Clear the pending provenance in `completeDIADaPublication`; flip three assertions; no re-record. |
| A112-02 | medium | `test-workflow-saved-da.R:1311` | `pending` and `provenance` masked in all 8 cases; only two differ. | Drop the masks after the fix. |
| A112-03 | high | `func_prot_da_model_stats.R:337-353`; `utils_workflow_dia_da.R:308` | F112-01: BH fallback recorded only in the store's contrast record; no analysis-record field, no export, no managed notice. Predates KIT-112. | Notice from the manifest; per-contrast method in the record. |
| A112-04 | medium | `utils_workflow_dia_da_validation.R:143` | F112-04 confirmed; fails after the arguments are committed. | Refuse all-NA or one-level terms. |
| A112-05 | medium | `mod_prot_da_handlers_run_analysis.R:68-79`; `utils_workflow_omics_actions.R:63-70` | Guards: cross-tab pending lets DA commit then resume the other tab's retry; retry ignores changed settings silently; retry can refuse forever if the head moved. By reading. | Third guard; notice; clear a superseded pending. |
| A112-06 | medium | `utils_workflow_dia_da.R:716-723`, `:782-796` | New with F112-03's fix: the reused Glimma page's identity excludes annotations. | Add the annotation identity. |
| A112-07 | medium | `skips-check.tsv:50`; `check.log:10` | All 8 identity blocks error in check mode. | Fixture root fix. |
| A112-08, -09 | low/info | | Session-class sniffing incl. `MockShinySession`; observer over 50 lines; 112m ended on battery. | KIT-109. |

Memory: the supportable claim is the DA step's own rise (+0.46 and +0.50 GiB in two in-process runs, +0.02 in one worker run; same-build spread 0.05). The absolute 2.208 to 1.575 includes an unexplained 0.15 to 0.36 GiB offset arising in `design`. DA is not the app's peak step; the app's maximum on this corpus is unchanged. Slow early steps confirmed, cause inconclusive.

### KIT-113 (enrichment on the worker): partly / partly / partly

| Id | Sev | Evidence | Finding | Smallest repair |
|---|---|---|---|---|
| A113-01 | high | `kit-113/w0/enrichment-heap-groups.tsv` | "+0.66 GiB" in one of three cases (others +0.0009, +0.0004); decision case ran first with an order effect. Robust: no benefit in any case. | Reword P-43 and both outcome files. |
| A113-02 | high | same, `live_after` | Managed app keeps +602 MiB (3 contrasts) or +207 MiB (1) after the Run, against +7 to +16 MiB in process. | Add to P-43. |
| A113-03 | high | double rows; ticket D5 | Without plot objects in the results the managed app peaks at 0.08 GiB against 1.08 to 1.16. The commit of plot-carrying results (+0.94 to +1.73 GiB transient) is the real enrichment peak on both paths. | Fourth P-43 option, own ticket. |
| A113-04 | medium | ticket D8 | The lead followed D8's letter; D8 had no "worse" branch and its no-worse rule could not fail. Defensible; a stop at W0 was better. | Future D8s need a "worse" branch. |
| A113-05 | medium | `utils_workflow_execution.R:1161`; `func_general_enrichment_engine_helpers.R:17,117-179` | Two `gost` calls per contrast, no timeout; a deadline kill mid-write leaves mixed files under the previous record. | Reword P-41/P-42. |
| A113-06 | medium | `mod_omic_enrich.R:590-596` | Another tab's pending publication blocks enrichment with "already running". | Distinct notice; with P-39. |
| A113-07 | medium | `test-omic-enrich-module.R:2573-2596` | The applied organism filter completes in no test on either path. | With P-07. |
| A113-08..13 | low/info | | Exclusions asserted exactly; transfer values may accumulate in the store (inconclusive); `runOmicEnrichAnalysis` about 105 lines; skip rebuild verified; check `--no-tests`. | — |

Open from KIT-107's audit and not traced: `recordWorkflowExternalResource` returns early inside a worker (`utils_workflow_dia_state.R:2924`), so the database version may not be recorded on managed sessions.

### Series level

| Id | Sev | Finding |
|---|---|---|
| AX-01 | high | `tickets/kit/KIT-109.md` names none of the queued items; its scoping numbers are stale. Re-author from HEAD. |
| AX-02 | high | Last whole-suite run: KIT-105 run 8 at `98da5939`. Seven tickets since. The handover's "none ran inside this series" is wrong (KIT-102 to KIT-105 each ran all files). |
| AX-03 | medium | All 12 `check.log` files are `--no-tests`. Overlaid check-mode error blocks grew 393 to 411. |
| AX-04 | medium | Overlay order and filtering verified for every ticket. Load_all passes across full runs fell 33,602 to 30,645 through quarantines (all 99 Q-ids are in the manifest). |
| AX-05 | medium | `m` member slice peak 1.957 to 1.974 GiB on 10 runs; 111m 2.096; 112m 2.062. |
| AX-06 | low | Power-state records exist for all 77 members; always low-power; start load of 2 or more on many members. Timings are weak evidence. |
| AX-07 | info | Six lane re-freezes, each attributed, line counts match. |
| AX-09 | medium | F112-01, F112-05, F108-01 should not wait for a removal ticket. |
| AX-10 | medium | KIT-110's refusal-then-retry note has no id and no disposition. |
| AX-11 | low | Lead-decided behaviour changes with no review row: F107-03, the blank-taxon silent stop, F111-03, F113-04/06, KIT-111 W2 error-path writes, F102-01. None numeric. |
| AX-12 | info | 58 finding ids: none lost. P-34, P-36, P-38 were contingent rows never raised. |

Worktrees (139 entries, nothing removed): 58 under `.claude/worktrees/` merged and clean (4.9 GB), 1 with two cherry-picked commits, 61 proof scratch trees merged and clean, 13 with one modified tracked file, 5 others clean. Look before deleting `scratch-081-run` (a modified helper and an untracked test file). 38 `worktree-agent-*` branches hold no work.

## 3. Decisions for the user

"Lead" is the lead's recommendation in `gui-workover-parked-decisions.md`.

| Item | Lead | Audit recommendation | Why |
|---|---|---|---|
| F112-05 | Defer to KIT-109 | **Fix now**, own commit | Wrong provenance in a user's record; one line; nothing to re-record. |
| F112-01 | Notice text, KIT-109 | **Fix now**: notice plus per-contrast method in the record | The method actually applied must be reportable. |
| F108-01 | KIT-109 | **Fix now** | A report can look complete after a failed save. |
| A111-01 | not known | Fix before KIT-109 closes, test first | A stale result can be committed as head. |
| Full suite | at KIT-109 | **Run now**, both modes, before any removal | 278 files unobserved since KIT-105. |
| Lipid/metab zeros | not raised | New parked decision; measure with KIT-114 | Same mechanism as F-18; parity. |
| Volcano default off | "your decision" | Confirm or reverse | No recorded approval. |
| P-01 | Measure first | Measure first, **redo the measurement** | Existing numbers are confounded. |
| P-02 | Keep E-level, default none | Keep default; label or withdraw quantile and scale until measured | No offered method is exact at E level. |
| P-03 | Message now | Message now; **settle before the real-data LFQ gate** | RUV may refuse on real LFQ. |
| P-04 | Re-record at the gate | Agree; count complete controls first | |
| P-05 | Keep each reader | **Refuse at import naming the column, all proteomics readers**; allow-list of vendor tokens to NA with a count | Silent NA is later imputed with no trace. |
| P-06 | Keep | Agree | File and oracle contract. |
| P-07 | Define | Agree, soon | Restores the original rule; mixed-species enrichment is unusable. |
| P-08 | Keep out of provenance | Agree | Derivable; affects no statistic. |
| P-09 | Keep both | Keep now; align the import to the builder at the next re-record | Export and re-import relabels contrasts. |
| P-10 | Keep none | Keep for the series; add afterwards | Parity; no valid design changes. |
| P-11 | Keep | **Add the preflight to metab import** | Only path that half-writes state. |
| P-12 | Keep | Reword at the next trx re-freeze | "Lipidomics" shown in the trx GUI. |
| P-13 | Keep | Agree | Display and order only. |
| P-14 | Keep | Align with lipid at a re-record point | Better provenance; low priority. |
| P-15 | Own store | Decide whether the applet is used: retire it, or own store plus a real table | It cannot start, by reading. |
| P-16 | Fixed | Agree; acknowledge | Status only; verified. |
| P-17 | Keep | Agree; correct the rationale | Alignment comes from per-step realignment. |
| P-18 | Keep per omic | One rule (lipid's) at a re-record point | Metab's rule can miscount; low stakes. |
| P-19 | Keep | **Candidate ids on both paths**, test first | In process can apply ticks to the wrong rows. |
| P-20 | As decided | Agree; superseded by P-32 | |
| P-21 | Keep | Agree | Suspended placeholders. |
| P-22, P-24, P-29 | Next harness revision | Agree; do it right after KIT-109 | KIT-109 makes more labels stale. |
| P-23 | As done | Agree; caveat DA comparisons across KIT-106 | Display work changed by up to 2.4%. |
| P-25 | Keep hidden | Agree; fix the wording before ever rendering | The note is wrong for lipid. |
| P-26 | Your call | **Build the download; hide the button until then** | Standing defect; parity. |
| P-27 | Drop | Agree, GUI and notebook together | |
| P-28 | Record | Agree; settle A107-03 and A107-06 in the same revision | The filter changes query and background. |
| P-30 | Your call (B, C, leave) | **B through the existing `reclaimWorkflowMemory()`**, measure once; retire the peak criterion | Already in the package; report is not the peak step. |
| P-31 | Drop | Agree | No consumer. |
| P-32 | Keep the exception | Agree; name it in KIT-109's closure | A container join shares nothing visible. |
| P-33 | Accept | Agree | Normalisation is not the session peak. |
| P-34 (new) | not raised | Raise; record the floor as intermittent | 2 of 4 members. |
| P-35 | Fix | Agree; cover the mode switch too | Double correction is a defect, not a method. |
| P-37 | Leave | Agree | No measured cost. |
| P-39 | Keep | Agree; add the third guard and notices | |
| P-40 | Accept | Accept as "no benefit shown" | Harness ticket only with the redesign. |
| P-41 | Accept | Accept after rewording; moot on revert | |
| P-42 | Accept | Accept only on revert; else stage | A killed worker can leave mixed files. |
| P-43 | Revert | **Revert; open the no-plot-objects redesign as its own ticket** | No benefit, persistent cost, new failure modes. |
| KIT-109 timing | Wait | Re-author now; hold removals until P-32, P-43, P-35, F112-05 are answered | |
| Worktrees | undecided | After KIT-113 closes; no `--force`; inspect the 13 modified | |
