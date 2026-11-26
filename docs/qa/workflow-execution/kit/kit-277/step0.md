# KIT-277 step 0 (implementer record)

Branch `kit-277` from origin/v0.9-staging **38443af8** (merge train 36 on top of 0bdf7213: KIT-266 touched
`R/utils_workflow_dia_norm.R`, `R/utils_workflow_dia_products.R`, `R/mod_prot_norm_workflow_helpers.R`, the limpa report;
none of the files this ticket edits). The ticket file (r1.2) is not edited; the repairs of AUDIT-KIT-277-plan-20261126
"Re-check 2" are recorded here, and the code and tests are built to the REPAIRED text.

## Bounded work item
Outcome: a formula-only, organism-only or FASTA-only design Save asks and is carried through (ticket sections 4.1 to 4.5,
three stages, gates SG1 (a run) and SG2 (by reading), R2 the comparison). Exit: acceptance criteria A1 to A8 of the ticket,
checked in `outcome.md`. Owner decisions: VAL-345, VAL-353 (ceiling 200.00 core-hours), VAL-377, VAL-366.

## Repairs recorded (exact replacement text of Re-check 2, section 6)
- **R11 (F1).** After the set-back the DA status is "none". F1's head expectation reads: Copy and the Report refused with
  the `workflow_summary_da_stale` condition (the DA files on disk named as made by no DA of this session); the Summary Save
  passes and its record and study parameters carry no DA block (the gate without `files` refuses only a stale DA,
  `R/utils_workflow_summary.R:153-162`). F4's "as F1" follows. Built as: the Copy/Report gate (`files = TRUE`) "refused",
  the Summary gate (`files = FALSE`) "allowed", the analysis record has no `differential_abundance` analysis.
- **R12 (coordinator's ruling, option (b)).** W288a and W288b in `tests/testthat/test-kit247-design-resave.R` are EDITED,
  not deleted (old and new text in `outcome.md`). The file (66.7 min, `merge-train-20261126v-tests.tsv`) is given an explicit
  120-minute per-file timeout in every request that runs it, is dealt to a shard that holds no other heavy file, and no
  shard is planned above 75 % of the 300-minute job limit.
- **R13.** SG1 runs the S1 rows and G2 together with `tests/testthat/test-omics-module-globals.R`: the stage-1 commit edits
  `R/mod_omic_da.R` and `R/mod_omic_design.R` (the scan covers them); the stage-1 filter is `MULTISCHOLAR_KIT277_STAGE=1`.
- **R15.** After option A's cut (P4) the head is not one of the DA's four states: Copy is refused by the DA gate
  (`workflow_summary_da_stale`; the status is "stale", the DA's source left the history); the Report, the Summary Save and
  Export Session State are refused by KIT-216's head rule before anything is written (VAL-377; not exercised: P4 asserts
  Copy, the DA status and the record). After a formula or organism Save the enrichment tab's Run is refused until the DA runs
  again (`omicEnrichManifestCurrent`, `R/mod_omic_enrich.R:623-637`); its text `enrichment_da_stale`
  (`R/utils_workflow_dia_sessions.R:735-737`) names a changed table as the reason: **recorded for the wording pass (VAL-41);
  the text is not changed.** O4 adds: the enrichment record stays in the states map outside the history (KIT-216); the row
  asserts the marker (`snapshot$enrichment_results` NULL) and the export's files, not the absence of the record.
- **R16.** KIT-279 (r1.1) adds no QC commit rule; its QC guard refuses a QC step whose state is in the active history, so
  after option A's cut re-applying the cleanup is a first application. P5 and P6 run the cleanup once per lane (each path on
  its own lane; no second run in one history).
- **R17, R18 (bookkeeping).** Stale statements: VAL-353 (1) ceiling 200.00 and (2) E-3 = (a) are decided; KIT-247 (train 29)
  and KIT-216 (train 35), KIT-268, KIT-278 are merged; the allow line of the globals test is
  `tests/testthat/test-omics-module-globals.R:59` (`shown_idle <<- idle`); section 13's step-0 list uses
  `R/utils_workflow_dia_sessions.R:1726` (`collectDIANativeGuiFiles`) and `:1782`; the Summary gate's call sites are
  `:298` (Summary), `:391` (Copy), `:525` (Report).

## Lines re-verified at 38443af8 (every cited line of Re-check 2's table; none moved since 0bdf7213)
`R/utils_workflow_summary.R` `assertWorkflowSummaryDaCurrent` :153, calls :298 :391 :525, `workflowSummarySaveStatus` :199,
`assertWorkflowSummarySaveCurrent` :220; `R/utils_workflow_dia_sessions.R` `collectDIANativeGuiFiles` :1726,
`enrichment_da_stale` :735; `R/utils_workflow_dia_state.R` marker :894, `selectWorkflowStateReference` :1699 (:1709),
`recordWorkflowInputs` :3074, `registerDIAInputFiles` :3120; `R/mod_omic_enrich.R` `omicEnrichManifestCurrent` :623, taxon
prefill :533; `R/mod_omic_da.R` `omicDaFollowHead` :145, run :362, load session :554; `R/mod_prot_da_handlers_run_analysis.R`
annotation :71, formula :72; `R/utils_workflow_dia_enrichment.R` resolver :425; `R/utils_workflow_omics_declarations.R`
manifest accept :2353; `R/utils_workflow_dia_summary.R` `protSummaryIntegrationWrite` :53 (annotation :59).
KIT-247's files are as merged (a52a808d): W288a :709-728 (loop over prot_lfq, lipid_canonical), W288b :730-753, C05 :755-775.

## Things found at step 0 that the ticket text does not say (no change of what must be built)
- The proteomics Import Existing Design reads its configuration from the PROJECT's `source_dir/config.ini`
  (`loadProtDesignImportedConfigAndTables`), not the folder's: F4's proteomics variant changes the formula in that file.
- `omicsDesignCompared` returns the remembered record or the bus; the remembered record now holds the formula (C03), so a
  session opened from a file (no commit in it) compares the formula with the bus's and cannot see a change made by an
  import. Recorded in "Not checked".
- The shell's restore of the session configuration reads `identical(flowResult, "changed")` (staging :545), not the
  `%in%` form of the ticket's text; it is edited to `"changed"` or `"limited"`.
- Worker path of P5: not drivable under `shiny::testServer` (no worker); the row asserts the project copy the worker reads
  is the applied table (equal to the bus table) and that the in-process run names the applied file's hash.
- Lipidomics/metabolomics session exports with a NULL `da_results`: not exercised (proteomics only).

## Budget (core-hours = job-minutes x 4 / 60; python3 -I)
Basis: newest value per file over every committed `*tests*.tsv` without "base" in its name (344 files); the 132 timed
grep files sum to **318.2** min (`test-kit216-prot-session-after-enrichment.R` now has a committed timing, 583.8 s, in
merge-train-20261127c's table, used instead of the 22.7 estimate); `test-omics-module-globals.R` 0.08; the new file 37.0
(an estimate): **355.3** per side, 134 files (grep 132, the new file, the globals test).
Shards: the job limit is 300 minutes, so four shards per side (not six: six were for an 80-minute limit); packing
`test-kit247-design-resave.R` (66.7, file timeout 120) with small files only, the new file (file timeout 90) likewise;
jobs 89.2/88.8/88.6/88.6 test minutes + 10 install, all far below 75 % of 300.

| Run | Jobs | Expected min | Expected core-h | Worst case core-h |
|---|---|---|---|---|
| R1 (new file; edited KIT-247 file) | 2 | 123.7 | 8.25 | 15.33 |
| SG1 (S1 rows + G2 + globals) | 1 | 40.0 | 2.67 | 6.67 |
| R2 (134 files, 4 shards x 2 sides + 2 check) | 10 | 870.6 | 58.04 | 92.00 |
| Expected total | | 1,034.3 | **68.95** | |

One correction pass (R1, SG1, R2 again) would add 68.96 expected. Running ceiling before each request: spent so far
(from the runs' job durations) + this request's worst case must stay at or below **200.00**; a request that would
cross it is not sent. Spent so far: 0.00.
