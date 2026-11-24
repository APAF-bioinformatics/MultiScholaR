# KIT-199: re-freeze of the GUI-entry reference `prot_dia_limpa` on run-tests, with a control record

Ticket `tickets/kit/KIT-199.md` revision 1.2 (the plan re-check's Q01 to Q06 applied at step 0). Owner decisions:
methods ledger VAL-109 (5) and VAL-120. Branch `kit-199-limpa-reference-refreeze`.

## Commits

| Name | Commit | Content |
|---|---|---|
| B | 19997323dd9a0efac5830696c7d3ac617560dae6 | `v0.9-staging` at the branch point |
| A | 12d68ca78442a50b46f52e33525dd85cf86df5e7 | the anchor (G1D's head, run 44953187101's base) |
| F1 | 84ec739515b1a44cc4dc8d1edc78d199f95016a7 | `g1a-reference-dia-limpa.tsv`, `declared-keys.tsv`, `control.R` (S1), `compare-da.R` (S3); `git diff --name-only B F1` lists these four |
| F2 | bf2b04d788c651fd4d875994129d552586511bc2 | parent F1; `tests/testdata/gui-entry/prot_dia_limpa.json` only |
| R | ef38bf2b557f303a72110e2742e377e12a663f36 | records of R1 and R2 (holds `r1/da-runner.tsv`); `R/`, `src/`, `inst/` equal B's |

## Pre-request checks (6.2; reading, no compute), at B

| Check | Observed | Holds |
|---|---|---|
| `git log --format=%H 12d68ca7..B -- tests/testdata/gui-entry/prot_dia_limpa.json` | empty | yes |
| `git diff --quiet 12d68ca7 B -- tests/testthat/helper-gui-entry.R tests/testdata/e2e/prot_dia_limpa tests/testdata/gui-entry/uniprot` | exit 0 | yes |
| `prot_dia_limpa` entry of `tests/testdata/e2e/manifest.json` at 12d68ca7 and B (`git show`, parsed, compared) | equal | yes |
| Harness commit against which the coordinator reconciled ruling C1 (a) | `kit-159-validation-harness` 9d1e87a7313fb2356f8fdbd6a063ab93adef813e (2026-11-24 13:37 +11), the commit the plan re-check read | named |
| `failures.txt` SHA-256 | 1b918062f8c94a3c10859a965f5f1120ea3c692ae70181039846c66a43730bde; byte-identical to `failures-base.txt` | yes |
| Reference rows | ten rows, six steps, one per (step, class). Blocks at lines 1669-1683, 1685-1699, 1701-1753, 1755-1807, 1809-1861, 1863-1915 (the ticket's "1670-1914" counts from the first block's `test:` line to the last block's last trace line; same blocks). `da` and `summary_refusal` hold only the DA key; each summary step holds the DA key and seven integration keys (four `.parquet` live only, three `.tsv` stored only); no `states`, `outputs`, history, tab status, notification or `<store entries>` line | yes |

`v0.9-staging` moved from 7d90a4c8 (the audit's reading) to B by three ledger-only commits.

## Observed

### R1 and R2: the control (runs 45315987153 and 45315990378, F1, reference BLAS)

- Both scripts exit 0 (652 s and 742 s). Every check of `r1/verdict.tsv` and `r2/verdict.tsv` is TRUE (25 checks each):
  the three sides ran the lane inside `test_that()` and passed; each raw DA copy's masked digest equals its record's at
  `da`; `steps-anchor.json`, `steps-base.json`, `steps-head.json` and `steps-frozen.json` byte-identical and parsed
  identical (sha256 ed5f2b87…); the raw DA copies of anchor, base and head byte-identical (sha256 007dea03…).
- Link (i), `frozen_vs_base`: 34 rows, exactly the reference's keys (six DA rows classed `declared_change (KIT-199)`,
  base value 8a4c2dd5…, stored value 80cfedaa…; 28 integration rows classed `integration_artifacts`, four steps of
  seven keys with the absent side as the reference says); no `outside` row.
- `anchor_vs_base`, `base_vs_head`, `head_vs_frozen_new`: 0 rows each.
- Negative control: the comparator reports exactly the two planted rows (`da` DA key, `norm` state `initial`), both
  `outside` (`r1/negative-control.tsv`, `r2/negative-control.tsv`).
- No `mechanism` stage; no `withheld.txt` in either run.
- CPU models: R1 AMD EPYC 9V74 (4 cores), R2 AMD EPYC 7763 (4 cores). Both: image ubuntu22 20261116.365.1, glibc
  2.35, libblas3 and liblapack3 3.10.0-2ubuntu1, R 4.6.1, Rust toolchain 1.100.0, internal pins fastduckdb e724f90b…
  and fastlimpadpc 0ad8c527…, `OPENBLAS_NUM_THREADS` and `OMP_NUM_THREADS` unset (`r1/environment.tsv`,
  `r2/environment.tsv`; they differ only in `cpu_model`).
- Agreement 4, compared on this machine with `cmp` and `sha256sum` on the fetched artifacts: `steps-anchor.json`,
  `steps-base.json`, `steps-head.json`, `steps-frozen.json`, `gui-entry-control.tsv` (4bda7326…), `da-runner.tsv`
  (2b2def97…) and `negative-control.tsv` identical between R1 and R2; the raw DA copies' sha256 identical
  (007dea03… on all six sides); `oracles/prot_dia_limpa.json` identical (`diff` empty, same `frozen` date).
  S6a and S6b do not apply; there is no `r4-s6b/`.

### F2: the freeze

R1's `oracles/prot_dia_limpa.json` committed as returned (its `source_sha` is F1, `frozen` 2026-11-24). `git diff
--stat F1 F2` lists that one file (24 insertions, 20 deletions). Inside it, the diff changes only: the DA key's value
80cfedaa… → 8a4c2dd5… at six steps; at the four summary steps the three `integration/proteomics/*.tsv` keys (content
digests) removed and four `*.parquet` keys (`<present>`) added, `corrected_proteomics.parquet` among them;
`source_sha` 5d381272 → 84ec7395; `frozen` 2026-11-17 → 2026-11-24. Its `steps` member written alone is
`steps-frozen.json`, byte-equal to `steps-head.json` in the job (R1 verdict).

### R3: the freeze commit under test (run 45319239164, F2 against B, check on both)

`test-gui-entry-proteomics.R`: base 183 passed / 10 failed (`prot_dia` 4, `prot_dia_limpa` 6), head 189 / 4
(`prot_dia` only; its four failure blocks identical on both sides), verdict `fixed` for `prot_dia_limpa`.
`test-e2e-proteomics-dia-limpa.R`: skipped on both (`same`). Check: ERRORs 0, WARNINGs 3, NOTEs 3 on both; findings
identical (6 on both, none only on one side). Records: `r3/compare.tsv`, `r3/check-compare.txt`, `r3/tests*.tsv`, `r3/failures.txt`, `r3/failures-base.txt`, `r3/compare-failures.txt`. Shards: base on AMD EPYC 7763, head on Intel Xeon Platinum 8573C; the head reproduces every capture of the re-frozen file on a third CPU model, of another vendor (masked captures only; the raw DA table was not compared there). `r3/failures-base.txt`'s six `prot_dia_limpa` blocks equal run 44953187101's text for text.

### R3b: the merge commit under test (run 45325909248)

Merge 6b3e615a60de619296d178490df939215c6ee892 (this branch at 89b72c9c merged with `v0.9-staging` 9453e79b; the only
conflict was the ledger, resolved by keeping staging's VAL-122 to VAL-130 and this ticket's row after them), pushed
to `run-k199-r3b`; base 9453e79b; `test-gui-entry-proteomics.R`, one shard, no check, reference BLAS. Base 183 passed /
10 failed, head 189 / 4 (`prot_dia` only; its four blocks identical on both sides); `prot_dia_limpa` `fixed`, no
failure at any step. Both shards on AMD EPYC 7763. Records: `r3b/`.

### R4: the route under OpenBLAS (run 45319240524, records commit R; D6, reported, decides nothing)

Measurement mode of F1's `control.R`, then `compare-da.R` (`r4/da-runner.tsv` = a, `r1/da-runner.tsv` = b). Job: AMD
EPYC 7763, OpenBLAS 0.3.20 (pthread), `OPENBLAS_NUM_THREADS` and `OMP_NUM_THREADS` unset; reference job R1: AMD EPYC
9V74, reference BLAS 3.10.0 (R1 and R2 identical, so R2 on EPYC 7763 is the same-CPU reference).

| Quantity | Result |
|---|---|
| Keys, order, header | 393 rows each, no key on one side only, order identical, header identical |
| Missing values | identical in every column |
| Calls (q < 0.05) | 14 up, 5 down on both; 0 changed |
| `log2FC` | 343 of 393 values differ in text; largest absolute 2.9e-13, relative 5.2e-11; 0 outside SWD-083 |
| `raw_pvalue` | 391 differ; largest absolute 2.5e-13, relative 5.4e-12; 0 outside |
| `fdr_qvalue` | 393 differ; largest absolute 1.7e-13, relative 4.3e-12; 0 outside |
| eight `log2norm.*`, eight `raw.*`, `annotation_score` | identical text (17 columns) |
| SWD-083 verdict | `agree` (overall; each moving column `agree`, the rest `identical`) |
| Masked DA digest under OpenBLAS | 80cfedaa… (the old stored value), not 8a4c2dd5… |

- The only difference the mask leaves visible is the notation of one value, the raw p-value of P00038 (0.00019627370306785395 under the reference BLAS, 1.9627370306790967e-4 under OpenBLAS, 2.8e-13 relative); every other difference lies in masked fractional digits. The masked digest therefore separates the two libraries on this table through that one value, not because the mask shows differences of this size in general (computed from the committed `r1/da-runner.tsv` and `r4/da-runner.tsv` with the lane's mask; shown for the 22 columns, the 15 text columns inferred unchanged).

### The Linux bundle's lane summary (no compute)

Latest `distributable.yml` run with the GUI-entry lanes: 45151293709 (2026-11-23, commit a4f60900, branch
`run-k172-bundle`; installed application, OpenBLAS 0.3.20). `prot_dia_limpa` on ubuntu-22.04 (job 137613133479) and
ubuntu-24.04 (job 137613133536): 34 passed, 4 failed, 24 steps, 4 steps differing, 28 discrete differences (the seven
integration keys at four steps), 0 digest differences: the DA key equalled the stored 80cfedaa…. Later bundle runs
(45297583402, 45273632957 and others) skipped the GUI-entry lanes.

## Inferred (not shown)

- After this merge the Linux bundle's lane (OpenBLAS) will differ from the reference at the six DA steps; no gate reads that step.
- The laptop's stored digest equals the OpenBLAS digest, so the laptop's masked DA text equals OpenBLAS's (shown, barring a hash collision), and its visible difference from the runner is the notation of P00038's raw p-value. That the laptop's numeric differences under the mask are of the same size as OpenBLAS's (about 1e-11 relative) is inferred only; the laptop measurement of 6.3 decides it.
- Note for KIT-159's G13 (audit I05): rounding to 10 significant digits does not absorb the OpenBLAS difference: 2 of the 1,179 values in the three moving columns change at the 10th significant digit, so a 10-digit fingerprint would report `differs` where SWD-083 agrees.
- D7's exact mode: R1 and R2 agreed byte for byte on EPYC 9V74 and EPYC 7763, so under the ticket's D7 text (as
  repaired by the implementation audit's I04) G1(a) compares exactly when its job lands on either model and its
  recorded environment (image, C library, BLAS and LAPACK packages and libraries, R, Rust toolchain, internal pins,
  unchanged lockfile) equals `environment.tsv`'s; otherwise under SWD-083.

## Acceptance criteria

| # | Met | Evidence |
|---|---|---|
| A1 | yes | pre-request table above; `g1a-reference-dia-limpa.tsv` and `declared-keys.tsv` in F1, before R1 |
| A2 | yes | `r1/negative-control.tsv`, `r2/negative-control.tsv` (two planted rows, both `outside`) |
| A3 | yes | `r1/verdict.tsv` all TRUE; `r1/gui-entry-control.tsv`; `r1/sha256.tsv`; `r1/environment.tsv` complete |
| A4 | yes | `r2/` all TRUE; agreement 4 masked and raw (above); CPU models of both jobs recorded |
| A5 | yes | `git diff --stat F1 F2`: one file; diff confined to the rows of 6.3; `steps-frozen.json` = `steps-head.json` |
| A6 | yes | `r3/compare.tsv` (`fixed`; `prot_dia` 4 on both; e2e `same`), `r3/check-compare.txt` (identical) |
| A7 | yes | this file, `runs.tsv`, `r1/`, `r2/`, `r3/`, `r4/`; longest path under this directory 69 bytes |
| A8 | yes (handed over) | ledger row `VAL-<merge>` after VAL-121; plan section 11 row and section 5 sentence: ticket 6.3, for the coordinator; no version change; no NEWS entry |
| A9 | yes | `runs.tsv`: 204.0 job-minutes, 13.60 core-hours against 17.67 expected and 72.0 ceiling |
| A10 | yes | R3b (run 45325909248) at merge 6b3e615a against 9453e79b: `prot_dia_limpa` `fixed`, `prot_dia` unchanged; `r3b/compare.tsv` |
| A11 | yes | `r4/compare-da.tsv`, `r4/measure.tsv`, `r4/environment.tsv`; table above; bundle summary read |

## Not checked

- The laptop's DA table and the cause of the laptop-against-runner difference (6.3, end-of-project laptop session).
- OpenBLAS with a set thread count, or on a CPU model other than EPYC 7763; the Windows bundle.
- Cross-model byte identity beyond EPYC 9V74 and EPYC 7763.
- Image drift after 20261116.365.1.
- Every item of ticket section 13 not named above.

## Merge record (coordinator)

Merged into `v0.9-staging` on 2026-11-24; the ledger row written `VAL-<merge>` above is VAL-136. R3b tested the merge with staging 9453e79b; staging had since gained KIT-175 (d2179e3b), whose own full suite (run 45315053703) showed `test-gui-entry-proteomics.R` failing with identical text on base and head, i.e. it does not move this lane; a confirmation run of that file on the merged staging head is ordered after this merge.
