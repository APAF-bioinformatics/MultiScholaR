# Merge train 27 of 2026-11-26: KIT-270, KIT-246c, KIT-221

**This train is branched from `origin/v0.9-staging` at `623f7550fc6efcd33ddbb28ebeddfe1c9da29b2c`** (fetched before
building; the two commits after 7f17547f are ledger-only, VAL-344). It holds trains 15 to 18 and 24 (KIT-213, VAL-341),
23 (KIT-259, VAL-342, VAL-343) as landed on staging; trains 19 to 22, 25 and 26 are other agents'.

## Train

Base: `623f7550fc6efcd33ddbb28ebeddfe1c9da29b2c` (staging).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-270 (test only) | `origin/kit-270` cdaceeed | becb3a82 | `VAL-<train27-1>` (`VAL-<merge>` on the branch) | AUDIT-KIT-270-impl-20261126 (ready) |
| 2 | KIT-246c (wording) | `origin/kit-246c` 09ccd396 | a01ad394 | `VAL-<train27-2>` (`VAL-<merge>` on the branch; the register calls the edits CN-1) | AUDIT-KIT-246-impl-20261126, closure C1 (ready) |
| 3 | KIT-221 | `origin/kit-221` 43107d2b | 3c4c10d0 | `VAL-<train27-3>` (`VAL-<merge>` on the branch) | AUDIT-KIT-221-impl-20261126 (ready) |

Merge order as asked, each `git merge --no-ff`. Record-fix commit (documents only): 47a9ad3f. **Head bound by the run:
`47a9ad3f42827ff4aed2e4050b9a5c9d37faa110`** (`R/`, `tests/`, `inst/`, `man/` identical to the third merge commit
3c4c10d0). This file is in the declaration commit after it, under `docs/` only.

### Conflicts met and their resolutions

- **`docs/methods-decision-ledger.md`** in each of the three merges, and only there: staging's rows VAL-330 to VAL-344
  against the member's row. In each, all of staging's rows kept in order; the member's row directly after VAL-344, then the
  previous members' rows, in merge order, as `VAL-<train27-1>` (KIT-270), `VAL-<train27-2>` (KIT-246c), `VAL-<train27-3>`
  (KIT-221). The coordinator numbers them at landing. KIT-246c's row is named `VAL-<merge>` on the branch (found: one row,
  "KIT-246 wording follow-up (VAL-305, after VAL-327)"); its record note CN-1 is applied (below).
- `NEWS.md` auto-merged (KIT-221's one bullet, after staging's last; KIT-246c and KIT-270 touch no NEWS).
- No conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`. `R/utils_integration_gui.R` auto-merged
  (staging's KIT-267 hunks are far from `:39-46`).
- **The merged code differs from the base by exactly the three members' changes** (SHOWN by `git diff --name-only
  623f7550 HEAD`, outside `docs/qa/`): `NEWS.md`, `R/func_prot_da_compat.R`, `R/func_prot_da_model_methods.R`,
  `R/utils_integration_gui.R`, `docs/methods-decision-ledger.md`, `tests/testthat/test-integration-record.R`,
  `tests/testthat/test-kit221-da-several-comparisons.R`, `tests/testthat/test-workflow-figure-paths.R`. Each member's own
  non-doc diff from its merge base: KIT-270 the one test file; KIT-246c `R/utils_integration_gui.R` and
  `test-integration-record.R`; KIT-221 `NEWS.md`, the two `R/func_prot_da_*` files and the new test file. The one
  file that staging also changed since the member's base (`R/utils_integration_gui.R`) carries exactly KIT-246c's hunk (same added and removed
  lines as the branch's own diff). No file outside the three members' lists differs.

### Record fixes (documents only, one commit before this declaration; no test and no `R/` edited)

- **KIT-221** (`kit-221/outcome.md`, exact texts from the register): I01 (dated correction: 46 files declared, 45 run;
  `test-metab-02q-da-analysis-helper-characterization.R` missed), I02 (the "In earlier releases" paragraph in "For the
  change report"), I03 (the two head observations), I08 (citation slips: the test header's `kit218MockPlots :86-91`
  (really `:79-83`) is a **test comment and is left as it is**; the R1 price and the "171" in `r2-predictions.md` are
  noted, not rewritten). I04: the NEWS refusal clause replaced with the register's text in `NEWS.md` **and marked in
  `outcome.md` as a DRAFT for the wording pass (VAL-41)**; the owner sets the final wording.
- **KIT-270**: R1 run id 45707283924 and the instrument re-send added to its ledger row (the register's J03 text).
- **KIT-246c**: CN-1's three edits applied to its ledger row (the merge commit, train and run id stay placeholders
  `<commit>`, `<n>`, `<id>` for the coordinator). CN-2 (three code comments) is an `R/` edit: **not done**, as the brief
  says no R edits.
- Not applied (landing, the coordinator's): the ledger numbers, VAL-188 (3) and VAL-193 status cells (KIT-221 I06),
  `kit-221/outcome.md`'s `VAL-<merge>` (acceptance row A5), plan and ticket bookkeeping.

## Declaration (written before the run is sent)

**Purpose.** One merge check of the three members on the merged tree: KIT-270's test edit, KIT-246c's wording passage,
KIT-221's several-comparison fix, against staging. The full suite is not repeated (KIT-221's A4 and its named
contribution to the DA-reaching files fall to the next checkpoint).

### Scope recount (merged tree against staging 623f7550)

- **KIT-221's grep** (`differentialAbundanceAnalysis|protein_deAnalysisWrapperFunction|deAnalysisWrapperFunction|
  daAnalysisWrapperFunction|getSignificantData|printCountDaGenesTable|runDIAReferenceDaAnalysis|
  func_prot_da_model_methods|func_prot_da_compat`, `helper-*` files dropped): **44 on the base, 45 on the head**; the
  difference is the new file `test-kit221-da-several-comparisons.R`. Plus `test-gui-entry-proteomics.R` and
  `test-workflow-operations.R` (not matched by the grep): **46 and 47**. The 47 include
  `test-metab-02q-da-analysis-helper-characterization.R` (I01) and `test-kit213-study-parameters-labels.R` (landed in
  train 24).
- **KIT-270's** two patterns find only `:99`, `:100`, `:471` of `test-workflow-figure-paths.R`: unchanged, one file, which
  is already in KIT-221's 47.
- **KIT-246c's** grep (`integrationMethodsText|integrationStudyParameterLines|integrationRecordTexts|methods_scaling|
  help_scaling|integrationScalingEvidenceText`): **13 files on base and head** (unchanged); the register's train scope
  (C3) is 11 files, the four-file wording run's files being `test-integration-record.R`, `-results.R`, `-gui.R`,
  `-model.R` (all four inside the 11). None of the 11 is in KIT-221's 47.
- **Union: 47 + 11 = 58 test files**, `check: true`. `test-prot-norm-core.R` and `test-omics-module-globals.R` are not in
  any member's scope.

### Settings

4 shards per side (listed order, heavy files first so they land in different shards), `file_timeout_minutes` 45,
`file_timeouts` `test-gui-entry-proteomics.R` 50, `check` true on both sides, no scripts (KIT-221's G1 needs none: neither
`R/` file changed on staging; the digests were equal on R2), `blas` reference, `ubuntu-22.04`, `fetch_refs` empty. **No
`env` switches** (KIT-246c's wording run set `MULTISCHOLAR_TRX_MANIFEST_NOTE`; none of the 58 files needs it, as the
KIT-278 register's plan says). **`extra_packages`: `["sva"]`, the only member request that had one** (KIT-221's R2); the
KIT-270 and KIT-246c requests had none. Base = `623f7550`, head = `47a9ad3f`; one run.

### Ceiling and expectation

**Ceiling 700 job-minutes = 46.67 core-hours** (core-hours = job-minutes x 4 / 60; 8 shards x 75 + 2 check jobs x 45 +
10). Expected about 420 job-minutes (about 28 core-hours): KIT-221's R2 took 342.35 job-minutes for its 45 files at 3
shards without the script job, plus the 11 integration files (about 25 test-minutes per side) and one more shard per side.
The scope fits the ceiling.

### Predictions

Base counts are the newest committed tables on staging (trains q, o, n, m, k; KIT-221's R2 `r2/compare.tsv`; train h for
the metabolomics helper file); trains 23 and 24 agree with them for every file they ran. **The train's own base run is the
reference where they differ: a base that differs from the figure is not a result; a head that differs from its own base
is.** Counts are passed/failed/errors/skipped. A file listed "failures identical" fails the same expectations on both
sides (stored references recorded on the laptop; KIT-156).

| File | Scope of | Base p/f/e/s | Head p/f/e/s | Verdict |
|---|---|---|---|---|
| `test-dia-limpa-ruv-rules.R` | 221 | 333/1/0/0 | 333/1/0/0 | same (failures identical) |
| `test-enrichment-significance.R` | 221 | 181/0/0/0 | 181/0/0/0 | same |
| `test-gui-entry-proteomics.R` | 221 | 92/101/0/0 | 92/101/0/0 | same (failures identical) |
| `test-integration-enrichment-mofa.R` | 246c | 0/0/0/3 | 0/0/0/3 | same |
| `test-integration-enrichment.R` | 246c | 349/0/0/0 | 349/0/0/0 | same |
| `test-integration-gui.R` | 246c | 405/0/0/0 | 405/0/0/0 | same |
| `test-integration-input-matrix.R` | 246c | 84/0/0/0 | 84/0/0/0 | same |
| `test-integration-model.R` | 246c | 265/0/0/0 | 265/0/0/0 | same |
| `test-integration-outputs.R` | 246c | 183/0/0/0 | 183/0/0/0 | same |
| `test-integration-record.R` | 246c | 200/0/0/0 | 209/0/0/0 | count-changed |
| `test-integration-removed-factors.R` | 246c | 520/0/0/3 | 520/0/0/3 | same |
| `test-integration-results.R` | 246c | 258/0/0/0 | 258/0/0/0 | same |
| `test-integration-workbook-api.R` | 246c | 144/0/0/0 | 144/0/0/0 | same |
| `test-kit128-no-biological-replication-prot.R` | 221 | 6/0/0/0 | 6/0/0/0 | same |
| `test-kit129-run-path.R` | 221 | 687/6/0/1 | 687/6/0/1 | same (failures identical) |
| `test-kit174-ruv-grouping.R` | 221 | 77/0/0/0 | 77/0/0/0 | same |
| `test-kit200-prot-da-no-value-proteins.R` | 221 | 538/0/0/0 | 538/0/0/0 | same |
| `test-kit201-treat-moderation.R` | 221 | 89/0/0/0 | 89/0/0/0 | same |
| `test-kit204-prot-da-not-tested.R` | 221 | 81/0/0/0 | 81/0/0/0 | same |
| `test-kit205-enrichment-background.R` | 221 | 228/0/0/0 | 228/0/0/0 | same |
| `test-kit208-gprofiler-identifiers.R` | 221 | 60/0/0/0 | 60/0/0/0 | same |
| `test-kit212-da-current.R` | 221 | 309/0/0/0 | 309/0/0/0 | same |
| `test-kit213-study-parameters-labels.R` | 221 (grep) | 56/0/0/0 (train q head) | 56/0/0/0 | same |
| `test-kit215-not-tested-lists.R` | 221 | 112/0/0/0 | 112/0/0/0 | same |
| `test-kit218-da-wrapper-tail.R` | 221 | 80/0/0/2 | 80/0/0/2 | same |
| `test-kit220-da-threshold-consumers.R` | 221 | 71/2/0/0 | 71/2/0/0 | same (failures identical) |
| `test-kit221-da-several-comparisons.R` | 221 | absent | 124/0/0/0 | only-on-head |
| `test-kit240-integration-default.R` | 246c | 96/0/0/1 | 96/0/0/1 | same |
| `test-kit243-summary-current-outputs.R` | 221 | 246/0/0/0 | 246/0/0/0 | same |
| `test-kit248-correlation-skip-record.R` | 221 | 206/0/0/0 | 206/0/0/0 | same |
| `test-kit249-prot-da-qvalue-fallback.R` | 221 | 144/0/0/0 | 144/0/0/0 | same |
| `test-lipid-07-s4-da-plot-methods.R` | 221 | 16/0/0/0 | 16/0/0/0 | same |
| `test-lipid-08-s4-da-analysis-methods.R` | 221 | 20/0/0/0 | 20/0/0/0 | same |
| `test-metab-02k-da-numsig-barplot-characterization.R` | 221 | 4/0/0/0 | 4/0/0/0 | same |
| `test-metab-02p-da-analysis-characterization.R` | 221 | 18/0/0/0 | 18/0/0/0 | same |
| `test-metab-02q-da-analysis-helper-characterization.R` | 221 (I01) | 19/0/0/0 (train h) | 19/0/0/0 | same |
| `test-metab-s4-da-methods-direct-shared.R` | 221 | 12/0/0/0 | 12/0/0/0 | same |
| `test-module-ci-prot-da.R` | 221 | 250/0/0/1 | 250/0/0/1 | same |
| `test-omic-da-module.R` | 221 | 1288/81/0/0 | 1288/81/0/0 | same (failures identical) |
| `test-omics-da-list-routing.R` | 221 | 5/0/0/0 | 5/0/0/0 | same |
| `test-prot-07-da-analysis.R` | 221 | 3/0/0/1 | 3/0/0/1 | same |
| `test-prot-07b-da-handlers-compat.R` | 221 | 35/0/0/2 | 35/0/0/2 | same |
| `test-prot-07c-da-results-characterization.R` | 221 | 23/0/0/0 | 23/0/0/0 | same |
| `test-prot-07c-da-results-output-shared.R` | 221 | 38/0/0/0 | 38/0/0/0 | same |
| `test-prot-07d-da-model-stats-shared.R` | 221 | 60/0/0/1 | 60/0/0/1 | same |
| `test-prot-annotation-uniprot-empty.R` | 221 | 327/0/0/0 | 327/0/0/0 | same |
| `test-prot-da-compat-shared.R` | 221 | 16/0/0/0 | 16/0/0/0 | same |
| `test-prot-design-log2.R` | 221 | 226/0/0/0 | 226/0/0/0 | same |
| `test-prot-ruv-factors-in-model.R` | 221 | 312/0/0/0 | 312/0/0/0 | same |
| `test-prot-s4-norm-methods-direct-shared.R` | 221 | 11/0/0/1 | 11/0/0/1 | same |
| `test-protein-da-golden-master.R` | 221 | 15/0/0/0 | 15/0/0/0 | same |
| `test-reference-mapping-enrichment.R` | 221 | 66/0/0/0 | 66/0/0/0 | same |
| `test-workflow-analysis-record.R` | 221 | 351/0/0/0 | 351/0/0/0 | same |
| `test-workflow-figure-paths.R` | 270+221 | 181/0/0/0 | 181/0/0/0 | same |
| `test-workflow-memory-mechanism-member.R` | 221 | 0/0/0/13 | 0/0/0/13 | same |
| `test-workflow-operations.R` | 221 | 75/0/0/0 | 75/0/0/0 | same |
| `test-workflow-saved-da.R` | 221 | 755/196/0/0 | 755/196/0/0 | same (failures identical) |
| `test-workflow-saved-installed.R` | 221 | 0/0/0/1 | 0/0/0/1 | same |

**Specific predictions** (from the registers):
- `test-kit221-da-several-comparisons.R`: absent on the base, **124/0/0/0 on the head**: `only-on-head`.
- `test-integration-record.R`: 200/0/0/0 to **209/0/0/0**, `count-changed` (KIT-246's block 33 to 42 passed).
- `test-workflow-figure-paths.R`: **181/0/0/0 on both sides, `same`.** KIT-270 removes, on the head only, the failure the
  bare "259" caused when a random temporary folder name contained 259 (about 0.5 percent of file runs). **Head side: it
  cannot fail this way. Base side: it can still flake once in about 190 runs** (`:99`, `:100` or `:471`, raw message
  quoting a folder containing 259): that is the known flake (C01), recorded, **not re-run**, and the verdict is then
  `fixed` (a base-only failure), not a regression.
- `test-workflow-saved-da.R` 755/196 same; `test-omic-da-module.R` 1288/81 same; `test-gui-entry-proteomics.R` 92/101
  same; `test-kit129-run-path.R` 687/6/0/1 same; `test-kit220-da-threshold-consumers.R` 71/2 same;
  `test-dia-limpa-ruv-rules.R` 333/1 same; `test-kit218-da-wrapper-tail.R` 80/0/0/2 same.
- Every other file `same`. **R CMD check: identical on both sides** (staging: 0 ERROR, 3 WARNING, 3 NOTE, 6 findings; one
  warning is the undeclared `cli`, AUDIT-KIT-269 I06), nothing only on head.

### Decision and attempts

- Every prediction holds: TRAIN GREEN for the three members.
- Any other outcome is a result: recorded with its files and counts, not repaired, not re-run, handed back.
- One run. A re-send only for an instrument defect with no result judged and only with the coordinator's word. An HTTP
  500/502 at dispatch with no run created: started once by hand with the README's command, recorded.

## What the coordinator resolves at landing

- The ledger: `VAL-<train27-1>`, `VAL-<train27-2>`, `VAL-<train27-3>` take their numbers (rows are in merge order after
  VAL-344). KIT-246c's row: fill `<commit>`, `<n>`, `<id>` (CN-1 edit 1), and INT-31 (its register C1: number INT-31 and
  fill `<commit>` in INT-31). KIT-270's status cell reads "implemented on `kit-270`; not merged".
- Status cells: **VAL-188 (3)** (append "implemented by KIT-221 (VAL-<n>), merged by train <t>", and plan audit E-2's
  qualification after "no computed value is wrong": "(for text input; with a comparisons column stored as a factor a value
  could carry another comparison's name: shown on base by KIT-221 W7)") and **VAL-193** ("(1) implemented by KIT-221,
  merged by train <t>; (2) KIT-230, parked (VAL-209)").
- `kit-221/outcome.md` acceptance row A5 names `VAL-<merge>`; KIT-221's NEWS refusal clause is a draft for the wording pass
  (I04); optional I05 clause (`enrichProteinsPathways()`).
- Records owed (no run): KIT-221 I05 (plan section 13 row beside F3) and I07 (the checkpoint names KIT-221 and KIT-218 for
  the DA-reaching files); KIT-270 ticket R02 sentences; KIT-246c CN-2 (three code comments, optional, an `R/` edit).
