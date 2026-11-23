# KIT-186 outcome: transcriptomics differential abundance without RUV

Ticket `tickets/kit/KIT-186.md` revision 1.3 (VAL-59; owner decisions VAL-66 D1 to D4 as recommended; VAL-74).
Branch `kit-186-trx-da-without-ruv` from `v0.9-staging` ffdfdfb7 (the branch point, "base" in R2).

**State: done, ready for its implementation audit.** R2 showed one verdict that revision 1.2 had not predicted:
`test-omic-da-module.R` `new-failure` at the frozen trx panel's whole-page digest (`:293`), which the help text W2 moves.
The owner decided (VAL-74) that the value is not edited here and waits for the end-of-project re-record (plan section 11
row, "Not here", beside `:312`), and revision 1.3 restates 6.2, A3 and A5 to admit exactly that one new failure; no
correction pass. R3 (full suite, with R CMD check at the coordinator's request) then gave exactly 6.2's verdicts.

## Commits

| Commit | What |
|---|---|
| 5d67e3dc | Tests only: `test-trx-da-without-ruv.R` (T1 to T6), T7 in `test-module-ci-trx-summary.R`, the four updated expectations (6.2), S1 script |
| 2060238d | Tests only: R1 instrument defect fixed (the journey helper reads the grid summary output only where the tab defines one); S1 judges A2 on A2's fields |
| e9e1df74 | The fix: reader and model text (7.1), W1, W2, comment in `mod_omic_da.R`, report W3 and W4, NEWS, ledger row `VAL-<merge>` |

## Runs (observed)

See `runs.tsv`. R 4.6.1, KIT-134 library, `blas: reference`, `ubuntu-22.04`, `extra_packages: ["sva"]`,
`MULTISCHOLAR_TRX_MANIFEST_NOTE=1`.

### R1 (run 45220978512, tests-only commit 5d67e3dc; `r1-tests.tsv`, `r1-failures.txt`)

- Every failure is an assertion; no test errored.
- T1: `out$da` NULL; one error notification "Run the normalisation first: the DA runs on the normalised (RUVs) state.";
  model text not "Formula: ~ 0 + group"; the record's `trx_ruvs` outcomes (`none`, `ruv_applied` FALSE, k 0,
  `edger_qlf_none`), head `trx_ruvs` and no `W_` column all passed. So the tab completes normalisation with `none` and
  refuses the DA, as 2.6 expected. T2 and T3: refusal notification, no DA operations, `out$da` NULL. T3's outcomes
  (`groups_all_features`, `num_sv_be`, k 0, `ruv_applied` FALSE) passed: the declared design with k forced to 0 completes.
- T4 (d) failed at the old text; T4 (a) to (c), T5 and T6 passed.
- T7: the base renders a report (`expect_length` passed) without a DA; the three required statements are absent and
  "RUVs factors W_1" and "inferred empirical negative control" are present (5 failures): the first stated base failure.
- Updated expectations fail at the old text: `test-omic-gui-config.R` (`trx$source$no_data`), `test-omic-da-module.R`
  `trx_refusals`/`run_no_data` notes, the refusal-order test, the frozen help (`(a)` and `(W3)`, help part).
- **Instrument defect.** `out$errors` held, beside the expected failure, an entry `norm`: "The test referenced an output
  that hasn't been defined yet: output$norm-ruv_summary_counts". The journey helper read the grid summary output
  unconditionally; the tab defines it only for today's grid (`R/utils_workflow_trx_norm.R:137-138`,
  `R/mod_omic_norm_panels.R:1193`), so every gate-on declared run records it, on base and head alike, and T1 to T4 could
  not pass on head. It is raised by test code, not the package (ticket section 9: an instrument defect under VAL-54).
  Fix (2060238d): the helper reads that output inside a handler that turns only that one condition into no summary. No
  expected value changed. One further instrument run (R1b) of the item.

### R1b (run 45225338786, 2060238d; `r1b-tests.tsv`, `r1b-failures.txt`)

- `test-trx-da-without-ruv.R`: 36 passed, 21 failed, 0 errors. Each failure at the value of R1 above; the `norm` entry is
  gone. The remaining `out$errors` entry `da` (empty message) is the helper's read of the DA summary output, which does
  not exist after a refusal; it is part of the stated base failure (no DA).

### R2 (run 45227884996, head e9e1df74 against base ffdfdfb7; `r2-compare.tsv`, `r2-compare-failures.txt`, `r2-check-compare.txt`, `r2-tests.tsv`, `s1/`)

Scope search re-run at the branch point (ticket section 9's `grep`): it finds the 23 listed test files it found at
5cbd4005 and no other; the request listed the section 9 scope, the three large files first.

| Verdict | Files |
|---|---|
| `only-on-head` | `test-trx-da-without-ruv.R` (57 passed, 0 failed) |
| `message-changed` | `test-gui-entry-transcriptomics.R` (13 failures on both sides) |
| `count-changed` | `test-module-ci-trx-summary.R` (T7: 6 passed; the 2 failures are MCI-032.1's laptop hashes, on both sides) |
| `new-failure` | `test-omic-da-module.R`: `(a)` and `(W3)`, the trx panel's `html` digest (below) |
| `same` | the other 26 files |

- GUI-entry pair: `diff` of the base and head blocks of `compare-failures.txt` differs only by the `da_refusal`
  notification lines (head: actual W1 against expected old text).
- **`test-omic-da-module.R` new-failure.** The frozen trx UI `omicDaW1FrozenUi$trx$html`
  (`tests/testthat/test-omic-da-module.R:293`, "001abe86...dd56c7ac") is a SHA-256 of the rendered panel, which holds
  the help text W2 changes. Head: "56d6c2da...0b10df04". The help part itself (`:312`, updated) passes on head, and the
  other 22 head failures of the file are the base's 22 (runner-against-laptop captures of lipid, metab and trx runs),
  and the updated notification captures (`trx_refusals`, the refusal order) pass on head.
  6.2 did not list `:293`; R1 could not show it, because the tests-only commit renders the old help.
- **A2 (S1, script ok, 261 s).** All six routes have k at least 1 on both sides (route 5's estimated count gave 3), so
  all six are A2 routes: 36 of 36 A2 values equal (DA table SHA-256, `formula_used`, k, model text, masked
  notifications, route 1's masked report-text SHA-256). The `errors` field differs on routes 4 to 6 by the base helper's
  `norm` entry only (recorded, not an A2 field, as written in the script before the run).
- **A4.** Check base and head: 0 errors, 3 warnings, 3 notes each; identical findings, none only on head.

## Acceptance criteria

| # | Status | Evidence |
|---|---|---|
| A1 | Met, with the R1 instrument defect recorded | R1, R1b (failing first at stated values); R2 head: the new file passes, T7 passes, the four updated expectations pass |
| A2 | Met | `s1/ruv-routes-compare.tsv`: 36/36 |
| A3 | Met under r1.3 (VAL-74) | R2: every verdict as 6.2 predicts; the one `new-failure` is the admitted `:293` digest |
| A4 | Met | `r2-check-compare.txt`; R3's check also identical (`r3-check-compare.txt`) |
| A5 | Met under r1.3 | R3 below |
| A6 | Met in the diff | NEWS under the unreleased heading; ledger row `VAL-<merge>`; no `DESCRIPTION`, `packaging/`, `man/` or `NAMESPACE` change |
| A7 | Met (read) | below |

## A7: the other omics' DA readers, re-read at ffdfdfb7

- Proteomics DIA, label-free, TMT: the DA tab accepts by state name (`R/mod_prot_da_handlers_init.R:19`, `:151`). RUV
  `skip` commits `normalized`; a declared design with k = 0 keeps `ruv_corrected` with `ruv_applied = FALSE`
  (`R/utils_workflow_dia_products.R:789-799`). DIA with limpa without RUV commits `ruv_covariates` (`:1469`, `:1483`).
- Lipidomics and metabolomics: `readLipidWorkflowMetadata` / `readMetabWorkflowMetadata`
  (`R/utils_workflow_lipid_declarations.R:125`, `R/utils_workflow_metab_declarations.R:119`) behind the head fallback
  that needs the normalisation tab complete (`R/mod_omic_da.R:286-289`; `completeOmicNormTab`,
  `R/mod_omic_norm.R:414-420`), whatever k.
- Outside transcriptomics, `"W_` occurs in `R/` only in the normalisation help (`R/mod_omic_norm_panels.R:113`) and the
  trx DA help (`R/utils_omic_gui_config.R:2907`). No other omic has the rule; nothing changed.

## Inferred, not observed

- That the trx panel's `html` digest changes only because of W2: the help part is the only other trx UI part that
  changed between base and head, and all other trx UI parts pass on head. The rendered HTML was not compared text to text.
- That the correction would be one value, `:293` set to the head digest "56d6c2da95d26603a47eb607893f694b1e70098a6262af04766fbe660b10df04".

### R3 (run 45241570517, head af7f0738 against base ffdfdfb7; `r3-full-suite-compare.tsv`, `r3-full-suite-compare-failures.txt`, `r3-check-compare.txt`)

Full suite, 20 shards per side, one process per file, with R CMD check (the coordinator's request; the ticket had
`check: false`). Head af7f0738 is e9e1df74's code plus these records. No file timed out or crashed on either side.

| Verdict | Files | Class |
|---|---|---|
| `same` | 445 | unchanged |
| `only-on-head` | `test-trx-da-without-ruv.R` (57 passed, 0 failed) | this ticket's new tests |
| `message-changed` | `test-gui-entry-transcriptomics.R` | W1: the base and head blocks differ only by the `da_refusal` notification (head: W1 actual, old text expected; `trx_screen.json:907`, section 11 row) |
| `count-changed` | `test-module-ci-trx-summary.R` (71 to 77 passed; the 2 failures, MCI-032.1's laptop hashes, on both sides) | T7 added (D3) |
| `new-failure` | `test-omic-da-module.R`: tests `(a)` and `(W3)`, both at the trx panel `html` digest only (actual `56d6c2da...0b10df04`, stored `001abe86...dd56c7ac`) | W2 (VAL-74), `:293`, section 11 row |

R CMD check: base and head 0 errors, 3 warnings, 3 notes; identical findings, none only on head.

### R2b (run 45264765395, merge 4433f86c against v0.9-staging d12f04ae; `r2b-compare.tsv`, `r2b-compare-failures.txt`, `r2b-check-compare.txt`, `s1-r2b/`)

Section 6.4: `v0.9-staging` gained code after R2 (0.9.5: KIT-172, KIT-180, KIT-173). Merged with a true merge commit
(4433f86c); the only hand-resolved conflicts were `NEWS.md` (the unreleased section kept above 0.9.5) and the ledger
(`VAL-<merge>` after VAL-78); no hand-resolved hunk under `R/`, `tests/`, `inst/`, `DESCRIPTION` or `NAMESPACE`. The
incoming code changes `R/utils_omic_gui_config.R`, a file of 7.3, so S1 was repeated. Same scope, order, shards, check.

| Verdict | Files | Class |
|---|---|---|
| `same` | 26 | unchanged |
| `only-on-head` | `test-trx-da-without-ruv.R` (57 passed, 0 failed) | new tests |
| `message-changed` | `test-gui-entry-transcriptomics.R` | W1: the blocks differ only by the `da_refusal` notification |
| `count-changed` | `test-module-ci-trx-summary.R` (71 to 77 passed; the same 2 failures on both sides) | T7 |
| `new-failure` | `test-omic-da-module.R`: `(a)` and `(W3)` at the trx panel `html` digest only | W2, VAL-74 (`:293`) |

`test-omic-da-module.R` has 51 failures on the new base (22 at ffdfdfb7); the head's are those 51 plus the two admitted
digest failures, every other one identical. Check: identical findings, none only on head. S1 (script ok, 254 s): all six
routes k at least 1 on both sides (k 2, 2, 2, 2, 3, 1), 36 of 36 A2 values equal.

Under 6.4, R3 is not repeated (R3b): the only hand-resolved conflicts were NEWS and the ledger, and R2b shows 6.2's
verdicts.

## Core-hours

R1 1.3, R1b 0.7, R2 16.8, R3 67.2: 86.0 for R1 to R3, against 85.0 expected and a ceiling of 456.4; R2b 15.5 (expected 19.5, ceiling 36.0). Total 101.5. R3 included the check jobs
the ticket's R3 did not (`check: true` at the coordinator's request).

## Not checked

- The report text of a declared design with k = 0 (W3's second sentence): rendered by no test; T3 checks the DA only.
- The RUVIII export after a k = 0 run (F3), the worker's normalisation log at k = 0 (F1).
- Windows; real data.
