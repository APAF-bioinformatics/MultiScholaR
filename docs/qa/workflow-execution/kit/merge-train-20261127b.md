# Merge train 35 of 2026-11-27: KIT-216 (full suite, checkpoint)

**Not stacked.** Branched from `origin/v0.9-staging` at `a855e1f0251286061e3b65e388fe63fc83b1dc9d` (fetched before building; the newest head: train 34 merged at b790119a,
then the ledger row VAL-376). Branch `merge-train-20261127b`, worktree `/home/doktersmol/Projects/APAF/MultiScholaR-train35`.

## Train

Base: `a855e1f0251286061e3b65e388fe63fc83b1dc9d` (staging).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-216 | `origin/kit-216` 8278ca44 (it already holds staging a855e1f0) | 761b4301 | `VAL-<train35-1>` (the member's row, was `VAL-<merge>`; the coordinator numbers it at landing, next free VAL-377) | AUDIT-KIT-216-impl-20261126 (ready after one bounded correction, run R2e 45759247038, held; VAL-372) |

Head (bound by the run): `761b43017` (full hash in the request) = the merge commit, the last code-bearing commit. The commits after it are under `docs/` only
(the ledger placeholder `VAL-<merge>` renamed `VAL-<train35-1>`, this record and its predictions table).

### Conflicts

None. `git merge --no-ff origin/kit-216` completed without a conflict (the member branch already carried staging a855e1f0, so NEWS and the ledger merged by themselves,
unlike the implementer's note that they conflicted at an earlier staging). The merge commit's tree is identical to `origin/kit-216`'s
(`git diff origin/kit-216 HEAD` is empty). The only edit after the merge: the ledger row's id.

### The merged tree against the base (`git diff a855e1f0 761b43017`)

Exactly KIT-216's change: 60 files, 23125 insertions, 121 deletions. Outside `docs/qa/workflow-execution/kit/kit-216/` (33 record files) and `tests/testdata/kit216/` (two archives):
- `R/` (11 files): `mod_omic_enrich.R`, `mod_omic_norm.R`, `mod_prot_da_handlers_load_session.R`, `mod_prot_norm_observer_helpers.R`, `utils_workflow_analysis_record.R`, `utils_workflow_dia_products.R`,
  `utils_workflow_dia_sessions.R`, `utils_workflow_dia_state.R`, `utils_workflow_dia_summary.R`, `utils_workflow_enrichment.R`, `utils_workflow_summary.R`.
- `NEWS.md`, `docs/methods-decision-ledger.md` (+1 row).
- Stored reference, edited by the member under its ticket's declared re-freeze (VAL-179 SD-3; not by this train): `tests/testdata/gui-entry/prot_dia.json`.
- Tests (14): `test-enrichment-significance.R`, `test-go-propagation.R`, `test-kit205-enrichment-background.R`, `test-kit206-prot-norm-export-guard.R`, `test-kit213-study-parameters-labels.R`,
  `test-omic-enrich-module.R`, `test-omic-summary-module.R`, `test-omics-enrich-core.R`, `test-workflow-analysis-record.R`, `test-workflow-integration-input.R`, and the new
  `test-kit216-prot-session-after-enrichment.R`; two test-data archives `tests/testdata/kit216/base-enrichment-1-run.msdia.tar` and `-2-runs.msdia.tar`.
- Nothing under `man/`, `NAMESPACE`, `DESCRIPTION` or `inst/`.

## Declaration (written before the run is sent)

**Purpose.** The R3 of KIT-216 (A3 of its ticket; AUDIT-KIT-216-impl F19 (i)): the full suite on the merge, base against head, with R CMD check on both sides.
**This run is also the checkpoint over the merge trains of 2026-11-26 and train 34** (trains a to z of 2026-11-26 and 27a; the last full-suite checkpoint is `checkpoint-20261126`, run 45644467847).
The predictions are therefore against the newest recorded count of every file (the table `merge-train-20261127b-predictions.tsv`, 535 rows: the checkpoint's table overlaid by the later trains' head counts, in train order).
**A base-side count that differs from the file's last recorded count is REPORTED as a finding, not explained away** (it is a result of the checkpoint, whichever cause the records suggest).

**Scope: the full suite**, `tests: "all"` (535 files on the head; 534 on the base, the new file is absent there). Confirmed from the tree now: `tests/testthat/test-kit213-study-parameters-labels.R` and
`test-kit247-design-resave.R`, `test-kit247-import-config.R`, `test-kit247-resave-no-session.R` exist and are included by "all"; the run's plan output (`summary.md`, the test file list) is checked after the run and
any of the four missing from either side is reported.

**Settings.** The form of `kit220-r3-full-1c2d8fdd.json` (the newest earlier full-suite request with a check and the variable's sibling settings; also `kit215-r3-full-94d57634.json`): `tests` "all"; **20 shards per side**; `file_timeout_minutes` 45;
`file_timeouts`: `test-gui-entry-proteomics.R`, `test-gui-entry-lipidomics.R`, `test-gui-entry-metabolomics.R`, `test-gui-entry-transcriptomics.R` and `test-prot-norm-core.R` 120 (as those requests),
`test-workflow-saved-da.R` 90 (as kit220's), and added here `test-kit247-design-resave.R` 100 (its recorded time is up to 67.2 minutes); `extra_packages` none; `env` `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`
(transcriptomics files are in a full suite); `check` true; `blas` reference; `os` ubuntu-22.04; no scripts. **No per-shard timeout field exists** (a test job has the fixed limit of 300 minutes; only `file_timeouts`
and `file_timeout_minutes` are per file). The 20 shards are the runner's round-robin of the sorted file list; computed from the head's list, the six heavy files fall in different shards
(`test-prot-norm-core.R` 5, `test-workflow-saved-da.R` 2, `test-omic-da-module.R` 9, `test-omic-enrich-module.R` 11, `test-omic-norm-module.R` 13, `test-omic-summary-module.R` 16, counting from 0) and the heaviest
shard (the one with `test-kit247-design-resave.R` and `test-workflow-saved-da.R`) is about 82 minutes of recorded file time, inside the 300-minute job limit. Base `a855e1f0`, head `761b43017`; one run.

**Ceiling.** 2,800 job-minutes = **186.67 core-hours** (core-hours = job-minutes x 4 / 60; the ticket's R3 row). Expected about 1,625 job-minutes = **108.3 core-hours** (audit F19: KIT-213's R3 1,423.87 plus about 100 minutes per side
for files added since); the recorded file times of this tree plus 10 minutes of setup per job suggest up to about 1,750.

### Predictions

Per file in `merge-train-20261127b-predictions.tsv` (file, recorded base and its source, predicted head, predicted verdict). Summary: **531 files `same`; 4 not `same`**; R CMD check identical.

| File | Base (recorded) | Head | Verdict |
|---|---|---|---|
| `test-kit216-prot-session-after-enrichment.R` | absent | 152/0/0/1 (the skip is the R1 archive maker) | `only-on-head` |
| `test-omics-enrich-core.R` | 221/20/0/1 | 222/20/0/1 | `message-changed`: only the W2 `kernel_clusterprofiler` reader's run `states` (base digest on base, `"<unchanged>"` on head); F-16 test 4 from 14 to 15 passes (R2e) |
| `test-omic-enrich-module.R` | 1011/59/0/0 | 1011/59/0/0 | `message-changed`: only at W1 `kernel_clusterprofiler` run `states` (F16: base digest `d8f57c99312ad90c` against the stored `1167a2b0dbbd17a5`; head `"<unchanged>"`); see the flake below |
| `test-workflow-integration-input.R` | 241/0/0/0 | 240/0/0/0 | `count-changed` (the SD-1 rename of one test, as in R2 and R2c) |
| R CMD check | 0 ERROR, 3 WARNING, 3 NOTE | identical | nothing only on head |

Named `same` (counts as recorded): `test-kit213-study-parameters-labels.R` 56/0 (E4's printed lines identical on both sides; the head's R2e result); `test-gui-entry-proteomics.R` 92/101 (the re-frozen leaves equal the
head; the KIT-217 digests and artifact lines are identical on both sides); `test-omic-summary-module.R` 1675/28; `test-omics-summary-core.R` 676/25/0/2; `test-kit206-prot-norm-export-guard.R` 116/0; `test-workflow-analysis-record.R` 351/0;
`test-kit205-enrichment-background.R` 228/0; `test-enrichment-significance.R`; `test-go-propagation.R`; the KIT-247 files 188/0, 2/0, 7/0; the transcriptomics files with the variable set (`test-module-ci-trx-summary.R` 87/1,
`test-workflow-trx-declarations.R` 313/31, `test-kit129-run-path.R`, `test-kit224-ruv-feature-levels.R`, `test-integration-input-transcriptomics.R`).

**Standing failures stay as they are** (38 files with a recorded failure or error; each `same` on both sides with its failures identical): among them `test-prot-norm-core.R` 1547/307, `test-workflow-saved-da.R` 755/196,
`test-omic-da-module.R` 1288/81, `test-gui-entry-proteomics.R` 92/101 (and the lipidomics, metabolomics and transcriptomics GUI-entry files), `test-omic-enrich-module.R` 1011/59, `test-omic-summary-module.R` 1675/28,
`test-trx-kernels.R` 238/144, `test-trx-publication.R` 92/82. A file whose failure list differs between the sides is a result.

**Known differences between sources, stated in advance.** `test-omic-norm-module.R` is recorded 639/0 in the newest tables (train v of 2026-11-26; train r's 631/8 was its head, train v's base and head both 639/0);
the coordinator's note and the audit's appendix say 631/8. The run's own base side decides, and a base count of 631/8 would be reported as a finding. `test-module-ci-trx-qc.R` is recorded 34/1/1/0 (the audit's appendix 38/0/1/0): same rule.
Eight test files added since `checkpoint-20261126` have their recorded count from the later trains' tables (see the table's source column).

**The `stringdb_plot` flake (`test-omic-enrich-module.R`).** The known one-sided digest flake (AUDIT-KIT-216-impl I06; train y): a `stringdb_plot` digest expectation in "(b) in process: refusals" that fails on one
side only, so the file reads 1011/59 on one side and 1011/60 or 1010/60 on the other, with the extra failure record naming `stringdb_plot` and no other record differing. It is recognised by exactly that: the only record present on one side
and absent on the other mentions `stringdb_plot`, in that block; the count differs by one on that side only. Recorded as the flake, no re-run. Any other one-sided record, or any change at another digest, is a result.

**Decision and attempts.** Every prediction holds: TRAIN GREEN and A3 met for KIT-216. A miss, a new failure, a message change outside the four rows, or a base count that differs from its recorded one is a result, recorded and not re-run.
No test or expectation is removed, weakened or edited. A shard timeout or an instrument failure with no result: hand back with the evidence; no re-send on my own. An HTTP 500/502 at dispatch with no run created is started once by hand.
