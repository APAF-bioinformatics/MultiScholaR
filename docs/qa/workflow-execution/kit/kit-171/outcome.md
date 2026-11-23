# KIT-171 outcome: MaxQuant reader, the legacy `Contaminant` column

Ticket revision 1.2. Branch `kit-171-maxquant-contaminant` from `v0.9-staging` 9a53d791. Every run on `run-tests`
(R 4.6.1, reference BLAS, ubuntu-22.04); run ids, commits and job minutes in `runs.tsv`.

| Commit | Content |
|---|---|
| c3aa6f46 | tests only: T1 to T4 in `tests/testthat/test-prot-import-maxquant-direct-shared.R`; script S1 `real-file-counts.R`; `R/` equal to base |
| 770155ba | "filter": the filter block filters `Potential.contaminant`, `Contaminant`, `Reverse`, logs rows removed per column and the "no contaminant column" case |
| 575067a0 | "typing": `Contaminant` added to the typed flag columns; roxygen, `man/importMaxQuantData.Rd` (by hand), NEWS |

## Failing first (A1), R1 on c3aa6f46, and the fix, R2 on 575067a0

| Test | Stated base value | Observed at c3aa6f46 | At 575067a0 |
|---|---|---|---|
| T1 | 4,999 (only row 4321 absent) | 4,999; absent only `REV__P04321`; no log lines (4 failures, 0 errors) | pass: 4,997, rows 1234, 3333, 4321 absent, rows 2222 and 4000 kept; log 2 `Contaminant`, 1 `Reverse` |
| T2 | 299, no `Contaminant` log line | 299; no `Contaminant` and no `Reverse` line (3 failures) | pass: 298, row 11 present, runs H1 and L1, one line each |
| T3 | 298 | 298; no log lines (3 failures) | pass: 297; log 2 `Potential.contaminant`, 1 `Contaminant` |
| T4 | 49 ids, no such line | 49 ids (pass); the line absent (1 failure) | pass |

All failures at the tests-only commit are expectation failures; no test errored. Logs:
`r1-tests-only.log`, `r2-fix.log`.

## Middle commit (A2), R1b on 770155ba

T4 green. T2 and T3 green. T1 red at 4,999 with rows 1234 and 3333 present, and no log line reports 2 rows removed in
`Contaminant` (by the code, the line reports 0; the test captures T1's log lines, so no console shows them): the
column was guessed logical by vroom and its marks lost (vroom's parsing-issue warning in T1). The typing change is
shown red on its own by T1. Log: `r1b-filter-only.log`.

## Real file (A5), script S1

`real-file-counts.tsv` (both runs side by side; the raw outputs `real-file-counts-r1.tsv`, `-r2.tsv`). Zip
38,428,469 B, SHA-256 557a66cc...d7133c in both runs; member 27,782,803 B, SHA-256 d245d4dc...1411ea.

| Count | c3aa6f46 (R1) | 575067a0 (R2) | Ticket |
|---|---|---|---|
| distinct protein groups, `filter_contaminants = TRUE` | 6,610 | 6,566 | 6,610 / 6,566 |
| distinct protein groups, `FALSE` | 6,694 | 6,694 | 6,694 |
| logged rows removed, `Contaminant` / `Reverse` | absent (base logs nothing) | 44 / 84 | 44 / 84 at the fix |
| runs | H1,H2,H3,L1,L2,L3 | same | H1..H3, L1..L3 |
| contaminant-marked groups with `_ECOLI` in `Fasta headers` | 0 | 0 | 0 |

No contaminant-marked group has an `_ECOLI` entry in `Fasta headers` (no positive control was counted).

Parsing checks (text read): 0 double-quote characters; header 112 fields, 0 data lines with another count; 6,694
data rows, 0 duplicated and 0 empty `Majority protein IDs`; 6,695 CR bytes (CR LF endings). Flags: `Reverse` 84,
`Contaminant` 44, both 0, `Only identified by site` 80. The base count equals the field-count parser's 6,610, so
the stop rule of section 8 did not apply.

Prefix counts (2.3 inferred about 7 and 71): groups with a `CON__` member and no `Contaminant` mark: 0 in `Majority
protein IDs`, 7 in `Protein IDs` only. Groups with a `REV__` member and no `Reverse` mark: 1 in `Majority protein
IDs`, 70 in `Protein IDs` only. Marked contaminant groups: 44 with a `CON__` member in `Majority protein IDs`, 22
leading with it, 22 made of `CON__` members only (the 22/22 split of 2.3). Marked reverse groups: 84, all leading
with and made only of `REV__` members. No unmarked group leads with a tagged member.

## Regressions (A3, A4, A6, A8)

- R2 (base 9a53d791 against 575067a0, the nine files of section 9): 8 files `same`, the reader's file
  `count-changed` (36 to 52 passed, 0 failed). `test-gui-entry-proteomics.R` (22) and `test-workflow-kit-harness.R`
  (1) fail identically on base and head (stored references recorded on the laptop; `r2-compare-failures.txt`: no
  non-identical failure); `test-e2e-proteomics-tmt-lfq.R` skipped its 3 blocks on both. `R CMD check`: base and head
  0 ERROR, 3 WARNING, 3 NOTE; findings only on head 0 (`check-compare.txt`).
- R3 (full suite, 448 files): 447 `same`, 1 `count-changed` (the reader's file, the four new tests), no
  `new-failure` or `message-changed`; 36 files fail identically on both sides (`r3-full-suite-compare.tsv`,
  `r3-full-suite-compare-failures.txt`). No stored reference of any omic moves: the `prot_lfq` entry lane's
  failures are the same on base and head.
- One vroom parsing-issue warning at the fix commit in T2: the untyped `Only identified by site` column (out of
  scope, 2.4), not a failure.

The 36 files failing on both sides of R3 are base state, with equal counts and byte-identical failure text on both
sides: 27 are stored or frozen references not re-recorded, 9 are runner-environment failures (missing `devtools` or
`lobstr`, `systemd-inhibit`, a bash diagnostic, an installed-library path, the optional-package load check).

## Not checked

- A positive control for the `_ECOLI` count: the rows whose `Fasta headers` contain `_ECOLI` at all were not counted,
  so the zero cannot be told apart from a column without that token.
- The rule that a line is logged also when 0 rows are removed: implemented, but no test asserts it.
- `R CMD check` at the records commit: the check ran at 575067a0, before the evidence files existed.
- Windows: no lane (by design, section 9).
- The vroom warning on the untyped `Only identified by site` column remains (stated above; out of scope).

## Budget

Core-hours (job minutes x 4 / 60, set-up included): R1 1.21, R1b 0.65, R2 10.54, R3 61.02; total 73.4 against 90.1
expected and 230.7 ceiling without correction. No correction pass used. Longest R3 job 45.0 min (ceiling 70).
