# KIT-100 outcome: on-disk (partitioned) defaults wherever not worse, 2026-11-07

**The four steps the user picked now default to `partitioned`, in the worker and in process, and the record names the variant that ran. Duplicate resolution stays `whole` on both omics. The three screen journeys are accepted and exact on `kit_100`, each recording `partitioned` on the changed steps. They ran under protocol v1.9 (the user, 2026-11-07: "1.9"), which widens ticket keys to three digits. Tests, dependents, `R CMD check` and the install are done. KIT-100 is closed.**

User decision (ticket revision 2.0, 2026-11-07): "yes", widened to parity with proteomics.

## The change (`53e66b4a`)

| Step | Default before | Default now | Measured ratio (partitioned / whole, large) |
|---|---|---|---|
| `lipid_log2` | whole | **partitioned** | 0.925 (KIT-084) |
| `metab_log2` | whole | **partitioned** | 0.954 (KIT-084) |
| `trx_expression_filter` | whole | **partitioned** | 0.941 (KIT-099) |
| `trx_tmm_logcpm` | whole | **partitioned** | 0.979 (KIT-099) |
| `lipid_duplicate_resolution`, `metab_duplicate_resolution` | whole | whole (parity; metab's is 1.042) | |
| `lipid_intensity_filter`, `metab_intensity_filter` | partitioned (KIT-089) | unchanged | |

- **Registry.** Each binding's `default_variants` field (KIT-094) in `R/utils_workflow_lipid_declarations.R`, `R/utils_workflow_metab_declarations.R` and `R/utils_workflow_trx_declarations.R`.
- **Worker.** An unqualified launch runs the registry default (`workflowQualificationVariant`, KIT-008), so no worker change was needed.
- **In process.**
  - Lipid and metab log2 and the trx filter go through `runWorkflowOmicsReferenceStep`, which already passes a default that isn't `whole` (KIT-089).
  - The trx logCPM didn't go through it. `runTrxReferenceRuvs` ran it through `trxRunProduct` → `runSavedOperation` with `list()`, so it would have run and recorded `whole`. `trxRunProduct` now applies KIT-089's rule: with no `variant` given, a default other than `whole` is passed. Its provenance is the RUVs state's `preceded_by` entry for `trx_tmm_logcpm`.
  - The managed path (`omicsManagedProduct("trx_tmm_logcpm", …)`) is the worker's.
- **Identity.** Unchanged: `variant` is not a declared parameter.
- **`whole` stays selectable** (`variant = "whole"`, or a qualification).

## Checks

- **Test first ([test-first.txt](test-first.txt)).** Two new tests:
  - `test-workflow-partitions.R` "KIT-100: log2 runs partitioned by default…": failed 4 of 14 before the change, passes 14 of 14 after.
  - `test-workflow-trx-declarations.R` "KIT-100: the expression filter and TMM/logCPM run partitioned by default…": failed 4 of 16 before, passes 16 of 16 after.
  - Together they check, per step:
    - the worker default (`workflowQualificationVariant(list(), op)`);
    - the variant recorded in process;
    - the result identical to `whole`: the scrubbed state for log2; the trx inputs for the filter; the logCPM samples and matrix, also equal to KIT-040's `log_cpm_tmm` fingerprint;
    - duplicate resolution, and the trx grid, RUVs, DA and RUVIII steps, staying `whole`.
- **Pins that move.**
  - `test-workflow-lipid-declarations.R`, `test-workflow-metab-declarations.R`, `test-workflow-omics-declarations.R`, `test-workflow-trx-declarations.R`: the expected `default_variant`.
  - `test-workflow-partitions.R`: KIT-079's default pin. In KIT-089's test, the "stays whole" example moved from log2 to duplicate resolution.
- **KIT-094 characterisation.**
  - `kit095WithoutTranscriptomics` now also names the KIT-100 log2 changes, step by step. For each omic, the log2 registry default and the log2 step's recorded variant in the fixture journey must be `partitioned`; each is then reset to the snapshot's `whole`.
  - Before the exception, `char-diff.R` found exactly these four differences, and nothing else.
  - The snapshot fixture is unchanged. Transcriptomics is outside the snapshot; its defaults are pinned in `test-workflow-trx-declarations.R`.
- **Module tests faking kernels.** None needed a `whole` pin: no dependent failed. `kit089WholeIntensity()` was not extended to log2.
- **Dependents ([dependents.txt](dependents.txt), [dep.sh](dep.sh)).** 90 files: every test naming the defaults, variants or the changed wrappers, and every lipid, metab and trx module and declaration test.
  - 17,215 pass, 0 fail, 0 error, 142 skip.
  - Each file's skip count equals KIT-099's checkpoint.
- **`R CMD check` ([check.log](check.log)).** `--no-tests --no-manual --no-build-vignettes`, `kit_100` library: 3 WARNINGs, 3 NOTEs, the same sections as KIT-099's baseline.
- **Install.** `kit_100` = `$P/core/browser-smoke.2jjdqX`, from a worktree pinned at `53e66b4a` (`scripts-100/install.sh`, base library `kit_098`'s).
- **Scripts.** [scripts/](scripts/) (install, check, bind), with `dep.sh`, `one.R` and `char-diff.R` here.
- **Bindings (`12190be1`, [bindings-diff.txt](bindings-diff.txt)).** `kit_100` (source `53e66b4a`) and the keys `100a` and `100b` (`094a`/`094b` on `kit_100`, against `lipid86s`/`metab86s`) and `100c` (`098a` on `kit_100` without the notebook reference, against `trx98a`). Nothing else changed except the generation timestamps; the protocol and harness hashes are unchanged. An earlier draft of this file, committed in `db70c3a2` before the bind finished, said the keys were absent. The keys are inert: the runner refuses their stage names (below).

## Protocol v1.9 and the screen journeys

**The selector question.** Protocol v1.8 accepted only `0NN` ticket keys, so the runner refused `kit-journey-100a` at launch: no member ran and no attempt was used. The user chose amendment v1.9 on 2026-11-07 ("1.9"; `protocol.md` section 11).
- `5bd87081`: `0[0-9]{2}` becomes `[0-9]{3}` at the ten sites (`run_workflow_kit.sh`, `workflow_kit_coordinator.R`, `workflow_kit.R`).
- The harness test `KIT-100 (protocol v1.9)` failed before the change (`kitSelectorStage` refused `kit:journey-100a:candidate:1`). After the change the harness file passes 733/0 and `test-workflow-analysis-record.R` 239/0. The trx, lipid and metab declaration files that source the harness pass 324/0, 259/0 and 259/0 on `kit_100`.
- While writing it, one test assertion had the gui-mode selector's id wrong (`member`, not `candidate-1`). That was existing behaviour, so the test was corrected and the code was not changed.
- `14388cc7`: the rebind changed only `protocol_sha256` and the four harness file hashes.

**Journeys** (one attempt each, 1,800 s ceiling, `scripts/all.sh`; on battery in the low-power profile, not changed):

| Key | Corpus | Oracle | Decision | Exact | Elapsed | Member anon peak |
|---|---|---|---|---|---|---|
| `100a` | lipid_screen_missing | lipid86s | accepted | yes, every step | 269 s | 1.84 GiB |
| `100b` | metab_screen_missing | metab86s | accepted | yes, every step | 262 s | 1.87 GiB |
| `100c` | trx_screen, 098a's settings | trx98a | accepted | yes, every step | 181 s | 1.92 GiB |

**Recorded variants ([variants.txt](variants.txt), from each member's analysis record):**
- `lipid_log2`, `metab_log2`, `trx_expression_filter` and `trx_tmm_logcpm` ran `partitioned`;
- both intensity filters ran `partitioned`, unchanged since KIT-089;
- both duplicate resolutions ran `whole`, and so did every other step.

The stage summaries are `kit100-kit-journey-100[abc]-stage.json`.

## Rollback

`git revert 53e66b4a` and the binding commit `12190be1`. `whole` stays selectable throughout.
