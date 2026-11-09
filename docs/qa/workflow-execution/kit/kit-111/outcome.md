# KIT-111 outcome: proteomics' normalisation tab on the worker, measured (2026-11-10)

**Done.** On a managed session, proteomics' four whole-matrix normalisation-tab declarations now run on the worker's `omics` role: `normalization-ruv`, `pearson-pairs`, `pearson-selection` and `protein_ruv_covariates`. Without a managed application, everything runs in process exactly as before.

**Results and records unchanged.** Every state, product, provenance entry, record, file and entry lane is identical. The only re-frozen values are the 86 attributed to F106-02 and F106-03.

**The memory claim is not made.** On the 300-sample DIA corpus the kernel left the application: 0.13 to 0.18 GiB lower while it runs. The step's peak (2.62 GiB) is unchanged, because it comes from the publication half that stays in process. That gives D8's declared "no worse, not better" outcome: the move is kept, and P-33 is raised ([memory-outcome.md](memory-outcome.md)).

- **How it is built:**
  - each in-process runner is split, keeping its name and signature, into a compute half (worker-safe), a commit description and a publication half;
  - the DIA binding adaptor gains `norm` entries, and the tab dispatches through the existing omics runner;
  - `completeWorkflowOmicsAction` commits the worker value's `commit` description;
  - one optional registry flag, `named_input` (F111-03), lets these four take a named saved input in place of the head.
- **Defects fixed, none changing a method or a result:**
  - **F111-01:** an omics cache hit on the worker failed ("live payload"). This affects every omic's repeated action.
  - **F111-03:** the worker refused a non-head input.
  - **F106-02:** the in-process proteomics DA left its working notice.
  - **F106-03:** the no-clusters summary printed a trailing `NULL`.
- **Build:** `kit_111` = `core/browser-smoke.MlQroA` (`c568dfa5`; R as `341dab60`). `R/` +599/−221 lines in 10 files ([lines.tsv](lines.tsv)).

## Work items

| Item | Commit(s) | What |
|---|---|---|
| W0 | `4971fae7`, `7c7c599e`, `a3b69c8d`, `e13c953f` | The large in-process baseline 111p on `kit_108` (M1): valid and complete; the verifier's 18 mismatches are KIT-065's set ([memory/baseline-large.md](memory/baseline-large.md)) |
| W1 | `280de978`, `ce660a68`, `8dd32418` | The in-process captures frozen on the start commit; the F106 tests failing first; probes. F111-01 reproduced; F111-02 found (P-35) |
| W2b | `0be353c1`, `704a9838`, `4a36f847` | F111-01 fixed, test first (`test-workflow-omics-worker.R` 283/0 → 289/0) |
| W2 | `85b5e1a1`, `7ba1cb8a`, `51461ac4` | The split, the handlers, the registry, the adaptor and the completion extension. F111-03 fixed (registry test 1015/4+1E → 1030/0). F111-04 retargeted |
| W3 | `476fd7c8`, `20329fe3`, `df9c601f` | The tab dispatches through the runner on a managed session (`test-prot-norm-core.R` 1541/0). F111-05 queued to KIT-109 |
| W5 | `5f7b6771`, `029d58e4` | F106-02 and F106-03 fixed. The one attributed re-freeze: 62 lane values, 24 module-capture values ([lane-changes.md](lane-changes.md)) |
| W4 | `341dab60`, `bcd3396b` | Merge checks; the bindings block; M0; `compare-stores.R`, checked on a known-identical pair |
| W6 | `5797d191` onward | The proof and the measurement, below |

## Proof

| Check | Result |
|---|---|
| Merge checks ([merge-checks.txt](merge-checks.txt)) | `Collate` 342; NAMESPACE identical; `man/` and `inst/` unchanged; back-references 2 |
| Simulation Setup ([proof/chain-a.md](proof/chain-a.md)) | `prot_dia` and `prot_dia_limpa` identical to their oracles; the traced runner, compute, commit and publish calls once each; the worker handlers 0 (in-process session) |
| Run 1, targeted ([targeted.txt](targeted.txt)) | 49 files, 15,753 passed, 0 failed, 0 errors |
| Run 3, entry suite | All 10 lanes identical, 0 step diffs: trx 29, lipid 35, metab 105, proteomics 193 |
| Run 4, DIA session and memory set ([dia-memory-set.txt](dia-memory-set.txt)) | 733/2, 191/1, 38/1, 174/4, 80/11, each equal |
| Run 5, `R CMD check` ([check.log](check.log)) | 3 WARNINGs, 3 NOTEs; the R-code NOTE is byte-identical to KIT-108's |
| Run 7, journeys ([proof/chain-b.md](proof/chain-b.md)) | 111l, m, t, a, b, c accepted and exact. 111l/m/t add exactly the moved operations, 1/1 each; 111a/b/c dispatch as 108a/b/c |
| M2, lifecycle | Completion as above. 111f fails at `normalisation_ruv` (`normalization-ruv` 1/0, `worker_failed`, head `protein_replicate_filtered`) |
| M3, large pair ([memory-outcome.md](memory-outcome.md)) | Stores identical: 16 states, 98 comparisons, 0 different. Same 18 verifier mismatches. "Better" not met; every "no worse" rule met |
| Run 8, e2e ([e2e-final.txt](e2e-final.txt)) | lipid 88, metab 190, trx 17, DIA 42, DIA_LIMPA 80, enrichment 98, TMT/LFQ 153: 0 failures, 0 skips |
| Run 9, skips ([skip-compare.md](skip-compare.md)) | Two changes, each accounted for: `test-workflow-execution.R` (the unset `NOT_CRAN`, rerun once, equal to the baseline) and `test-prot-norm-core.R`'s check-mode errors (all the known source-tree read) |

## Budget use

- **GUI proof:** one instrumentation rerun (`test-workflow-execution.R` with `NOT_CRAN=true`). The correction pass was not used.
- **Memory:** the correction pass was used once, as 111n, to rerun 111m, whose application floor (+0.10 to +0.13 GiB from `rollup` on) did not reproduce.
- **Host:**
  - every member under `systemd-inhibit` with no other `multischolar*` unit live ([power-state.txt](power-state.txt));
  - two read-only ticket agents (KIT-112) worked during 111q's window, recorded in `memory-outcome.md`;
  - W1's agent R runs overlapped 111p's post-run verification, not its measured window (recorded at W0).

## Hand-offs

- **To KIT-109:**
  - **F111-05:** the RUV covariates QC image paths in the staging copy.
  - **The diagnostic's exit boundary** (`workflow_checkpoint_gui_diagnostic.R`) fires at dispatch on a managed session (the W3 note).
- **To KIT-112:** reuse the split and completion pattern, the `named_input` flag and `compare-stores.R`. Check whether the controller peak is reachable before changing code (P-33's lesson).
- **Parked for the user:**
  - **P-32:** proteomics' normalisation tab and the shared frame: the ticket recommends a container-join ticket; the lead prefers the recorded exception.
  - **P-33:** the publication half to the worker, or accept.
  - **P-35:** a repeated DIA_LIMPA RUV covariates click re-applies RUV.
