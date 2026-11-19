# Harness v1.11 change set: runner evidence (2026-11-19)

## Run 44461036735: full suite and check

- Base: 039ce71d (`kit-129-integration`). Head: 4a7c211f (`kit-harness-v1.11`, commits b1 to b4, a1, a2).
- `run-tests` 8fe3e11b; 20 shards; `sva`; `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`; reference BLAS; R 4.6.1.
- **427 files**: 410 `same`, 9 `new-failure`, 6 `message-changed`, 1 `count-changed`, 1 `only-on-head`.
- `R CMD check`: 0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides, the same six findings.
- Failing-first (`kit-129/final-method/failing-first.R` on base 039ce71d):
  - the five tests of `test-workflow-folder-names.R` fail on base;
  - the new `test-workflow-kit-harness.R` block (protocol v1.11 item 1) is absent from base.
  - On head: folder names 60 passed, 0 failed; harness 773 passed, 1 failed. That one failure is the same on base (KIT-099 corpus digest, laptop-recorded).

**How each difference was classified.** Every non-identical failure message was compared with base's after mapping:

| New name | Mapped back to |
|---|---|
| `.n-<ID>` | `.dia-native-<ID>` |
| `.we` | `.workflow-entities` |
| `.wx` | `.workflow-export` |

A difference that disappears under that mapping is attributed to the folder names. The only lines left were:
- the step labels ("at step ..."), because the folder keys now differ from the first step on;
- context names that `testthat` prints around the changed keys;
- the key order (sorting `.n-` and `.dia-native-` puts the keys in a different order).

| File | Verdict | Class | Cause |
|---|---|---|---|
| `test-gui-entry-proteomics.R` (prot_dia, prot_tmt, prot_lfq, prot_lfq_fragpipe; prot_dia_limpa was already failing) | new-failure / message-changed | expected, re-record (R1) | `<store entries>` keys `.dia-native-<ID>/.workflow-entities` become `.n-<ID>/.we` (H6) |
| `test-gui-entry-lipidomics.R` (lipid_canonical) | new-failure | expected, re-record (R1) | H6 keys only |
| `test-gui-entry-metabolomics.R` (three lanes; already failing on runners from the cyclic-loess fixtures, KIT-134 defect 4) | message-changed | expected, re-record (R1) | H6 keys; the state differences are base's own |
| `test-gui-entry-transcriptomics.R` (trx_screen; already failing, KIT-128) | message-changed | expected, re-record (R1) | H6 keys; the rest is base's own |
| `test-omic-da-module.R` (new: trx_refusals, prot_refusals, prot_load, prot_main, prot_invalid_preflight, prot_no_contrasts; the lipid/metab/trx cases were already failing, KIT-128) | new-failure / message-changed | expected, re-record (R2) | H6 frozen `<store entries>` keys and their order |
| `test-omic-enrich-module.R` (8 cases; one error, `files[[state_entries]]` subscript out of bounds at the `..workflow-entities` literal, :3688) | new-failure | expected, re-record (R3) | H6 |
| `test-omics-enrich-core.R` (4 cases: the core compared with the module's W1 (b) captures) | new-failure | expected, re-record (R3, the same captures) | H6 |
| `test-omic-summary-module.R` (prot_owner, prot_owner_variants, prot_native new; trx cases already failing) | new-failure / message-changed | expected, re-record (R4) | H6 |
| `test-omics-summary-core.R` (trx_managed new; trx_main already failing) | new-failure / message-changed | expected, re-record (R4/R6) | H6 |
| `test-workflow-figure-paths.R` (the Summary session archive's limit: 64/125 becomes 77/138) | new-failure | expected, re-record (R5) | H6: the staging root `.wx` is 13 characters shorter |
| `test-module-ci-trx-da.R`, `-norm.R`, `-summary.R`; `test-workflow-trx-declarations.R` | message-changed | expected | H4: the same actual and expected digests on both sides; the message text changed (the expectation is now `expectTrxDeclDigest`; the manifest note is the rule) |
| `test-workflow-kit-harness.R` | count-changed | expected | H2: the new block, 39 expectations pass |
| `test-workflow-folder-names.R` | only-on-head | expected | H6: new file, 60 passed |
| **`test-workflow-saved-projection.R`** ("installed adapter members preserve exact consumers and close on late failure": "reader input integrity failure") | new-failure | **unexpected: a harness defect of a1** | The checkpoint driver's installed-adapter reader cases named their payload `.workflow-entities/...`, but the package wrote `.we`. Repaired in 04bc186c (the one correction pass); confirmed in the run below. |

Nothing unexpected is left beyond the saved-projection row, which was repaired.

## Run 44467185692: the correction

- Base 039ce71d; head 04bc186c.
- Files: saved-projection, folder-names, kit-harness, checkpoint-harness, execution-adapters, rollup-qualification, bulk-relational, memory-mechanism-member.
- Result:
  - saved-projection: 406 passed, 0 failed, `same` (the defect is repaired);
  - folder names: 60 passed, 0 failed;
  - kit-harness: 773 passed, 1 failed, the same failure as base;
  - checkpoint-harness and bulk-relational: `same`;
  - execution-adapters, rollup-qualification and memory-mechanism-member: skipped on runners on both sides (0 expectations; they need a systemd user slice or the installed vehicle). They are not exercised off the laptop.
