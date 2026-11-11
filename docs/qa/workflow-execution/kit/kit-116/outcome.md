# KIT-116: real-path verification and test repairs after the GUI workover audit (2026-11-11)

**Outcome:** done. W1 and W2 complete; W3's e2e lanes all pass and P-34 closes (0 of 6). One declared run (the real-worker probe) was not run, with the reason below. F116-01 was found and then fixed on the user's word; the deviation is reported. Every count below is `load_all`, one process per file, with `R_LIBS` at kit_920's library, unless a row says otherwise.

## W1, product code

| Item | Result | Commits | Evidence |
|---|---|---|---|
| P-02 interim (user, 2026-11-11) | Quantile and scale are labelled "approximate on limpa" on the DIA_LIMPA normalisation tab, with one line of help saying why (their standard errors are not rescaled after DPC). The choice values are unchanged. | `08e9b25e`; oracle `2ab897a0` | `test-module-ci-prot-dia-limpa.R` 38/0/2 (new block failed first). One entry oracle re-frozen, `prot_dia_limpa`: the only change is the QC tabs' UI hash at each step, plus the oracle's SHA and date. The lane is 38/0 against it. |
| AX-11, blank taxon (user, 2026-11-11) | An enrichment Run with a blank taxon now says "Enter the organism's taxonomy ID before running the enrichment analysis." and writes nothing, as before. A taxon of blanks only is now stopped there too; before, it went on and failed later. | `722b6126` | `test-omic-enrich-module.R` 1070/0 (new block failed first). |

## W2, tests and tooling

| Item | Result | Commits | Evidence |
|---|---|---|---|
| A102-05, A104-07 (Q-30, Q-64, Q-75) | **20 quarantined comparison blocks frozen and live again.** Their expected values were computed from the deleted old functions at the commits where each comparison last passed (`d469cea3`, `7b99e7e8`, `827f2328`), by scripts kept in `freeze/`. Each script refuses to run off its green commit and checks every old value against the core of that commit. The current core equals every frozen value; nothing is left skipped. Three expectations that tested only the old functions were dropped, each named in a comment. | `223a4b6b`, `8324547d`, `4eab0760`, register `b4016d62`, `459b9fea` | `test-workflow-omics-import.R` 24/0/8 to **290/0/0**; `test-omics-qc-core.R` 32/0/10 to **115/0/0**. |
| A103-03 | The two design parity tests that compared the shared code with itself now compare it with the old lipid and metab functions' values, frozen at `616e1147` (the parent of the commit that made them wrappers). | `f7bd77a4` | `test-omics-design-core.R` 61/0/10 (the 10 skips are Q-35/Q-37, outside this item). |
| A104-03, A105-06 | **Confirmed at the function level, then fixed (F116-01).** On DIA_LIMPA, `runProtNormTabEntryWorkflow` warns "Please complete the Quality Control filtering steps" at `peptide_count_filtered` and `ruv_covariates` whenever the QC tab status is not `complete`. That includes the status the DA tab's load-session leaves (`normalization = complete`, QC untouched). Whether a real restore lands on those two states was not shown: the DA load-session's head is normally `correlation_filtered`, where no warning is given. **Fixed on the user's word (2026-11-11, "fix any outstanding issues"):** on DIA_LIMPA those two states need no status, because they are reached only once the QC has finished; `ruv_corrected` keeps its pre-QC generation. | `5b3a93c9` | the probe's table: `peptide_count_filtered` and `ruv_covariates` warned with either status; `correlation_filtered` generated the pre-QC and did not warn. After the fix: `test-prot-05b-norm-module-contracts.R` 509/0/9 (six new cases failed first); `test-gui-entry-proteomics.R` 193/0; `test-prot-norm-core.R` 1690/0. |
| A101-10 | **The DIA_LIMPA entry lane detects a changed scientific default, but not a changed DA heatmap default.** See [mutation.md](mutation.md). | `1d5206be`, `70d95073` | attempt 1 (heatmap `top_n`) 38/0, not detected; attempt 2 (the normalisation default) 15/23, first data difference at `qc_protein_normalisation`. |
| A111-05 | The managed-against-unmanaged capture now also compares the bus's `config_list` (the RUV audit trail) and each committed state's whole `@args`, with a guard that neither is empty. | `d32082f8` | `test-prot-norm-core.R` 1674/0 (from 1642; the guard block 389/0). Managed and unmanaged are equal on both new fields in every case. |
| AX-10 (A110-01) | A managed trx refusal followed by a run in the same session completes as a first run does, with the frozen refusal and the frozen run's history, states, outcomes, tab status and files. There is no "live payload" error: KIT-111's F111-01 fix covers it. A coverage test: it did not fail first. | `3c2f9866` | `test-trx-norm-module.R` 20/0 for the new block. |
| A105-03 | `kit-105/skipcmp.py` takes an optional list of targeted files and stops when one has no row. It also reports a new file error and fewer static blocks. | `fad473cf` | checked on two small tables (a missing targeted file stops with exit 1; a parse error and a lost block are reported) and on the suite's own table against itself (no rows). |

## W3, runs

**Build:** `kit_921`, install `browser-smoke.uWhPRZ` at `aa729025`, from kit_920's library. Bindings at `ab4ee4cc` add the build and members 921a to 921f, nothing else.

**Not run, against the declaration: the real-worker probe for A111-01.** The question A111-01 answers is which input reference a queued click uses, and the controller decides that before dispatch. Two tests already prove it: R2's `2db1c3e5` (failed first) and KIT-111 W3 (R5), "a click while the tab's action runs waits for it". A separate worker process adds only the resolution of a named non-head input. Every journey member does that at correlation, the six 921 members included. A probe would need a new real-application test harness for no new evidence.

**e2e lanes** ([e2e.txt](e2e.txt), `6cdfbbc8`): the seven lanes once each on `kit_921`, the app on the install's library, the tests from a worktree at `aa729025`.

| Lane | Pass/fail | Wall | Against KIT-112 |
|---|---|---|---|
| lipidomics canonical | 88/0 | 130 s | 88 |
| metabolomics LC-GC | 190/0 | 281 s | 190 |
| transcriptomics screen | 16/0 | 89 s | 17: R5 removed the "no page before the tick" expectation (the volcano is on by default). **A106-02:** the R5-edited volcano block ran in a browser for the first time, and the managed page arrived. |
| proteomics DIA | 42/0 | 285 s | 42 |
| proteomics DIA_LIMPA | 86/0 | 216 s | 80: **A106-03,** the six new DA-view assertions (the volcano page, the heatmap image, the results table's rows, the CSV download with rows). |
| proteomics enrichment | 98/0 | 585 s | 98 |
| proteomics TMT/LFQ | 153/0 | 295 s | 153 |

**P-34** ([p34-cut.md](p34-cut.md) fixed the measure and the 0.80 GiB cut before any run; evidence in `p34/`, commits `461c074c` to `16c0a734`): six `m` members on `kit_921`, in sequence, nothing else live, I2 attached. All six were accepted and exact.

| Member | `rollup` min (GiB) | max | Step time | Verdict |
|---|---|---|---|---|
| 921a | 0.732 | 0.822 | 25 s | released |
| 921b | 0.722 | 0.813 | 25 s | released |
| 921c | 0.732 | 0.822 | 27 s | released |
| 921d | 0.730 | 0.864 | 28 s | released |
| 921e | 0.734 | 0.824 | 29 s | released |
| 921f | 0.734 | 0.824 | 28 s | released |

**Decision (as declared): 0 of 6 at or above the cut, so P-34 is closed as intermittent and rare.**
- **All members:** the floor was not released in 2 of 19 `m` members so far (111m, 112m).
- **Builds with KIT-111's code:** 2 of 10.
- **These six:** the minimums lie inside the released group (0.709 to 0.745), and the sampler traces are not needed.
- **Step time:** `rollup` took 25 to 29 s against 13 to 17 s before. The sampler ran beside it and the host load was 2 to 5. This is recorded, not judged; the memory measure is the minimum, not a time.

## Checks

- Each W1/W2 item's targeted file and dependents, as listed per row.
- `R CMD check --no-tests` on a build of `16c0a734`: 3 WARNINGs, 3 NOTEs, the series baseline.
- Install: `kit_921` (`browser-smoke.uWhPRZ`).
- Not run: the full suite. It is KIT-109's entry run, per the plan.

## For the user

- **The dropped probe** (above).
