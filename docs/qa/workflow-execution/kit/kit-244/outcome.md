# KIT-244 outcome: DIA with limpa, "Run Normalization & RUV" pass-through

Ticket r1.2 (`tickets/kit/KIT-244.md`). Part A implements the owner's decision VAL-220 (1). Part B is a guard that
changes no number. The body-digest re-pin is the owner's decision VAL-246. Branch `kit-244`, cut from `v0.9-staging`
16112917 (merge train 3, with KIT-225), then merged with staging 7327607f.

## Commits

| Commit | What |
|---|---|
| 5cc6b032 | G0 records (script and tables from the scripts-only branch `kit-244-g0`) |
| ffb90f3a | Tests only: `tests/testthat/test-kit244-limpa-pass-through.R` (nothing under `R/`) |
| e4a37f2e | Fix: `R/utils_workflow_dia_norm.R` (two conditions and the comment), `R/utils_workflow_dia_products.R` (the guard) |
| 5600e8c8 | R1/R2 records, NEWS entry, ledger row `VAL-<train8-1>` |
| 8fa690d0 | Re-pin of one body digest, `tests/testthat/test-prot-norm-core.R:45` (VAL-246) |
| 7daada75 | Merge of staging 7327607f (ledger conflict only: staging's VAL-230 to VAL-254 kept, this row after them) |

## Step 0 (by reading, at the branch point 16112917)

- **Identity with G0's binding.** These files are byte-identical between 600dcc6e and 16112917 (`git diff --stat` empty):
  - `R/utils_workflow_dia_norm.R`
  - `R/utils_workflow_dia_products.R`
  - `R/mod_prot_norm_session_helpers.R`
  - `R/utils_workflow_dia_state.R`
  - `tests/testthat/test-prot-norm-core.R`
  - `tests/testthat/helper-prot-norm-core.R`
- `R/mod_prot_norm_workflow_helpers.R` differs only by KIT-174 (a)'s two insertions (`omicNormWarnRuvGrouping`). Neither is on a path G0 ran:
  - The first is in `runProtNormRuvCovariatesObserver`, after the `skip` mode's `return`. G0 ran that observer in mode `skip` only.
  - The second is in the classic branch of `runProtNormNormalizationWorkflow`, after the `DIA_LIMPA` branch has returned. G0 pressed only DIA_LIMPA owners.
  - So G0 was not repeated.
- **Cited lines re-read.** Each cited line is where the ticket says:
  - pass-through `:597-619`, condition `:604`
  - `diaLimpaPassThroughSource` `:627-631`
  - guard `R/utils_workflow_dia_products.R:939-947`
  - body digests `test-prot-norm-core.R:43-45`, block `:156-167`, expectation `:165`
- **Tests that press the pass-through, and the head at each press.**
  - `test-prot-norm-core.R` (W1, W2, P-35): at L7, and again at the RUV states.
  - `test-prot-norm-input-current.R` `:481`, `:596`, `:599`, `:702`, `:761`, `:890`, `:904`: all at `peptide_count_filtered`.
  - `test-kit206-prot-norm-export-guard.R` `:321`, `:436`, `:615`: L7.
  - `test-workflow-saved-da.R:731`: L7.
  - `test-kit225-prot-norm-reset-source.R`: L7.
  - None presses at `correlation_filtered` or at a peptide head.
- **`_snaps` and `tests/testdata`:** no file holds the pass-through's texts.
- **Regression list.** The ticket's grep gives 41 files at 16112917: the 39 at 600dcc6e plus `test-kit174-ruv-grouping.R` and `test-kit200-prot-da-no-value-proteins.R`, which came with the train. Their union with KIT-225's R2 list (30 files) is 58 files; with the new file, 59.

## G0 (run 45583972146; binding 1a9ffe4c = `origin/kit-225` 600dcc6e plus the script, coordinator's ruling)

Tables: `g0-base-values.tsv`, `g0-best-partner.tsv`, `g0-pairs.tsv`. Every action ran in process and through the worker; the two paths agree on every row. There were no instrument errors.

**Observed.**
- **t_hi.** The pairs product on L7 has 30 within-group pairs, and every one of the 12 samples is in one. b(s) runs from 0.96258 to 0.97206.
  - m = 0.96257597326065381 (run_0001); m2 = 0.96464941373759827 (run_0004).
  - **t_hi = 0.96361269349912604.** It removes run_0001 only.
  - 0.5 removes nothing, so C1 (L7 filtered at 0.5) has L7's samples and the same value.
- **(A) at 0.5.**
  - Press 2 at `correlation_filtered` records C1-state, draws its figures from it, and `diaLimpaPassThroughSource` gives the same.
  - Correlation 0.5 again then reads C1-state, commits nothing, and gives the same value.
  - A Skip after press 2 commits `correlation_filtered` with `skipped = TRUE`. Reset after that Skip is refused (`prot_norm_reset_source_unknown`).
  - Reset after press 2 changes no state; its notice says "(the peptide quality-control steps)".
  - RUV set to Skip commits nothing and press 2 then records C1-state.
- **(A) at t_hi.**
  - A-hi: correlation 0.5 after press 2 reads Chi-state and keeps 11 samples (run_0001 absent).
  - A-skip: the Skip commits Chi's 11 samples; Reset is refused.
  - A-C2: Reset is refused.
- **B1.** At all five peptide heads the table is `PeptideQuantitativeData`.
  - In process, the press completes: it records the peptide state and draws figures from it.
  - Managed, the press queues one action.
- **B2** (states made by the tab's own observers; stop rule 5 not needed):
  - Apply: `workflow_operation_input`, "pearson-pairs expects ProteinQuantitativeData, the source is PeptideQuantitativeData".
  - Skip: commits the peptide table, then `saved_entity_unsupported`.
  - RUV covariates: `saved_entity_outputs_pending`, "correlation outputs are pending; retry that action".

**Stop rules.**
- Rule 1: not triggered. Press 2 recorded C1, and no peptide press was refused on base.
- Rule 5: not needed.
- Rule 6: there is a t_hi.

**Inferred.**
- E2b passes on base: there was no commit and the value was identical. It is kept as a G row, for that reason.

## R1 (run 45588177269, tests-only commit ffb90f3a, new file only)

- 124 passed and 80 failed, all in the E rows (E1, E1b, E2, E3/E3b, E4, E4b, E5/E5b). Each failure shows G0's base value: "C1-state", "Chi-state", 11 samples, `prot_norm_reset_source_unknown`, the "(the peptide quality-control steps)" notice, no refusal, and one queued action.
- No G row failed: G1, G2b, G3, G4, G5, G7 and E2b-as-G all passed.
- Tables: `r1/`.

## R2 (run 45591340671; base ffb90f3a, head e4a37f2e without the re-pin)

- **Settings.** 59 files, 4 shards a side, check. Env `MULTISCHOLAR_KIT105_W3_DUMP=1` and `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` on both sides.
- **Results.**
  - 57 files `same`.
  - The new file is `fixed`: 204 passed, 0 failed.
  - `test-prot-norm-core.R` has exactly one new failure, the body-digest block (`:165`), with one element differing: `runDIAReferenceLimpaNormalisationPassThrough` went from `fe852d7a…` to `e991f0d2aa83d45b4ee1df5484793a1ea60118e7ae59123401f5548d79293f79`. Its 306 base failures are identical in message (the file has one entry in `compare-failures.txt`).
  - The head's dump: the other 19 pinned values are unchanged.
- **Check.** Identical on both sides: 0 ERRORs, 3 WARNINGs, 3 NOTEs; nothing only on head.
- Tables: `r2/`.

## Re-pin (VAL-246) and R3 (run 45605539062; base staging 7327607f, head merge 7daada75)

- **The re-pin.** `test-prot-norm-core.R:45`: old `fe852d7a3086ab62b9bb8fadf6214da74668265e0d7bccf60d001d1eae01fbbf`, new `e991f0d2aa83d45b4ee1df5484793a1ea60118e7ae59123401f5548d79293f79`. The new value is taken from R2's head dump (run 45591340671), not computed locally. It is the only value edited, and the comment at `:43-45` names KIT-244 and VAL-246.
- **R3 settings.** 8 files, 2 shards a side, check:
  - `test-prot-norm-core.R`
  - `test-gui-entry-proteomics.R`
  - the new file
  - `test-kit225-prot-norm-reset-source.R`
  - `test-prot-norm-input-current.R`
  - `test-kit206-prot-norm-export-guard.R`
  - `test-workflow-saved-da.R`
  - `test-kit174-ruv-grouping.R` (staging changed `omicNormRuvHeadTechnicalReplicates`)
- **R3 results.**
  - All 7 shared files are `same`. `test-prot-norm-core.R` is 1548/306 on both sides; the body-digest block passes on both, and the head's dump equals the pinned values for 20 of 20 functions.
  - The new file is only on head: 204 passed, 0 failed.
  - The check is identical on both sides.
- **Effect on R2's failure.** The one expectation that failed by design in R2 now passes.
- Tables: `r3/`.

## Acceptance criteria

| # | Verdict | Evidence |
|---|---|---|
| A1 | met | G0 run 45583972146; `g0-base-values.tsv`; B1 cites KIT-225 G0-2 and was observed again here at all five heads; binding as ruled |
| A2 | met | R1: every E row fails at G0's value, every G row as base. R2 and R3: every E row passes on the fix. 05-core passes on R2's base side (the tests-only tree) |
| A3 | met | R2: every file `same` except the new file (`fixed`) and the body-digest expectation at `:165` (one element). The 306 base failures are identical. Nothing new only on head in the check. KIT-224 was not straddled: neither side contains it (it merged later, in merge train 6, c4505c7f). The merge train's two sides will both contain it, so its two `test-prot-norm-core.R` messages (VAL-213) are equal on both |
| A4 | pending | the merge train's run. With the re-pin committed, the block passes (shown in R3) |
| A5 | met | NEWS entry under the unreleased heading; ledger row `VAL-<train8-1>` after the last VAL row; `man/`, `NAMESPACE`, `DESCRIPTION` unchanged |

## Core-hours

| Run | Job-minutes | Core-hours |
|---|---|---|
| G0 | 17.1 | 1.14 |
| R1 | 18.6 | 1.24 |
| R2 | 506.5 | 33.77 |
| R3 | 208.7 | 13.91 |
| **Total** | **750.9** | **50.06** |

- Expected 33.8, ceiling 107.0.
- R2 ran above its expected 457 job-minutes, within its 660 ceiling. Its longest shard took 85.6 minutes, above the 75 per shard assumed.
- R3 was not in the ticket's plan (section 11: no own R3 unless the merge train's attribution fails). It was declared
  before it was sent (purpose, base 7327607f, head 7daada75, 8 files, ceiling 300 job-minutes), cost 208.7
  job-minutes, and keeps the total at 750.9 of the 1,605 job-minute ceiling. It is not a correction pass and does not
  replace A4. Coordinator's ruling (2026-11-26): accepted as such; it changed no deciding comparison; recorded as told
  to the owner, with the budget figure 50.06 core-hours against 33.8 expected.
- No VAL-54 instrument rerun was needed.
- Module CI runs started by this ticket's pushes (not requested; audit I11): 727.6 job-minutes over six runs, listed in
  `runs.tsv`. Their only failures are staging's transcriptomics list (MCI-030.1, -031.1, -032.1); the merge push's
  `e2e impacted (proteomics-dia-limpa)` lane passed (one press per journey).

## Wording (VAL-41; the owner sets it)

The refusal (class `prot_norm_limpa_peptide_input`) is shown after "Error in normalization: ". It is the ticket's draft, verbatim: "on DIA with limpa, \"Run Normalization & RUV\" needs the protein table that DPC Quantification makes, and the current table is still a peptide table. In the Quality Control tab apply the remaining peptide steps and DPC Quantification (normally followed by Between-Sample Normalisation and the Peptide Count Filter), then press \"Run Normalization & RUV\" again."

## Not checked

- No R was run on the laptop. Everything above ran on GitHub runners with R 4.6.1, reference BLAS and ubuntu-22.04.
- **Fixture.** On the "small" fixture, 0.5 removes no sample. The sample-level difference between before and after is therefore shown at t_hi only (one sample, run_0001). The 0.5 rows tell before from after by the recorded table and by Reset.
- **Running application.** A browser was not driven; the press at `correlation_filtered` in the running application is by reading only.
- **Other omics and routes.** Lipidomics, metabolomics and transcriptomics were not traced.
- **Other press sequences.** RUV-mode runs after press 2 were not run (their input is unchanged by reading).
- **Not this ticket.** Windows and real data; I08 (F2) and U1 (handed off).
- **U1.** After an Apply and then a Skip, the record keeps the Apply's threshold. For a second press followed by a Skip, the live table now differs from what replay recomputes until U1 is repaired.
- **F1, F3 to F5.** As in the ticket, section 12.
- **Merge train.** The train's full run (A4).
- **Replay after the lower-threshold sequence.** Replay after press 1, Apply at t_hi, press 2, Apply at 0.5 was not
  run. By reading, replay filters L7 at 0.5, the live table; KIT-248's G0b (45597013061) observed replay's source on
  this route. Two of the three links are observed (the live table, R2 E2; replay's source, KIT-248 G0b); the third,
  the record naming 0.5, is read from the code.
- **F6 (sessions saved before this change).** A session exported before this change after a second press keeps the
  correlation result as its recorded start; in it, re-applied filters still act at the first threshold, and Reset may
  be refused; the way out is Revert in Quality Control to Peptide Count Filter, then press. Development sessions only.

## For the change report

- **Where.** On the development branch only. The DIA with limpa route does not exist on `main` (ticket 2.5), so nothing released was affected.
- **What was wrong.**
  - (A) On DIA with limpa, after the correlation filter (or after RUV was set to Skip, which saves nothing), pressing "Run Normalization & RUV" again made the filtered table the new starting table. A filter re-applied at a lower threshold could not bring back the samples the first filter had removed, although the analysis record named the lower threshold. Skipping the filter kept those removals while recording "skipped". "Reset to Pre-Normalization" either changed nothing or was refused with a false reason.
  - (B) Pressing the button before DPC Quantification was accepted on a peptide table. The correlation filter and RUV covariates then failed with internal messages, and the RUV covariates said "retry", which cannot work.
- **Who was affected.** Only users of DIA with limpa (an opt-in route), and only for these sequences.
  - (A) changed results only when the first filter removed a sample and the second action was a lower threshold or a Skip. At the same or a higher threshold the samples kept were the same. The reason: a sample the first filter keeps keeps its best within-group partner, and the correlation between two kept samples does not change when other samples are removed (each pair uses only its own complete observations, and the filter removes samples, not proteins). So filtering the filtered table again at a second threshold keeps exactly the samples the unfiltered table keeps at the higher of the two thresholds. No sample was ever removed twice at the same threshold.
  - No result recorded in this programme went through the sequence.
- **Code.**
  - Base 16112917: `R/utils_workflow_dia_norm.R:604` and `:628` (the kept-input condition); `R/utils_workflow_dia_products.R:939-947` (the guard).
  - Fix e4a37f2e: `R/utils_workflow_dia_norm.R:606` and `:630` add `correlation_filtered`; `R/utils_workflow_dia_products.R:942-951` refuses a peptide table first.
- **What changed in results.**
  - (A) After the second press, a re-applied filter acts at the chosen threshold on the table the first press designated, and a Skip commits that table. This is what the record states and what replay recomputes (by reading for replay: two of its three links observed, see Not checked; a Skip after an Apply excepted: U1, repaired by KIT-248).
  - (B) No number changes.
- **Evidence.** G0 45583972146 (base values; t_hi = 0.96361 removes run_0001; at base a re-apply at 0.5 after t_hi kept 11 samples, after the fix 12); R1 45588177269; R2 45591340671; R3 45605539062. The fixture is the "small" synthetic recipe (12 samples). It shows the mechanism, not effect sizes on real data.
