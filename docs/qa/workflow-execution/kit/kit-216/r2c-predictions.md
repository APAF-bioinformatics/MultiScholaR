# KIT-216 R2c: predictions, committed before the run

The run compares two commits, with R CMD check:
- **Head:** 7865ab47, which is the correction pass 45b3026c merged with staging 4d588ba3.
- **Base:** 687bccf9 on branch `kit-216-r2c-base`, which is tests-only 4320e668 merged with the same staging commit.

The scope is the ticket's grep on 4d588ba3 (88 files) plus the extras, 91 files in all (`r2c-scope.txt`). That
matches the R2 audit, F12.

The R2 audit's prediction table (section 6) is adopted unchanged. Line numbers are those of the merge tree, and may
shift by a line or two after the correction's edits.

| File | Verdict | Where |
|---|---|---|
| `test-kit216-prot-session-after-enrichment.R` | `fixed` | Every E row fails on base as in R2 and passes on head. G21 and E9's restore check exist on head only, so they pass on head. 1 skip (the R1 archive maker) |
| `test-kit207-summary-head-guard.R` | `same` | G7 passes on both |
| `test-workflow-dia-products.R` | `same` | KIT-024 |
| `test-enrichment-significance.R` | `same` | S10 |
| `test-omic-enrich-module.R` | `new-failure` / `message-changed`, **only** at stored leaves | See the list below this table |
| `test-omics-enrich-core.R` | `new-failure` / `message-changed` | W2 run `history`/`states` (4 cases; readers `:103`, `:135`). F-16: 14 → 15 expectations |
| `test-omic-summary-module.R` | `new-failure` | `prot_native` `save`: `values$session_summary`, `notes`, `files`, `integration` (about `:4294`, `:4296`, `:4306`, `:4310`) |
| `test-omics-summary-core.R` | `new-failure` | the same `prot_native` `save` leaves (readers `:223`, `:230`) |
| `test-gui-entry-proteomics.R` | `message-changed` | `prot_dia`: see the list below this table. The pre-existing parquet/tsv artifact lines are the same on both sides |
| `test-workflow-integration-input.R` | `count-changed` | the test renamed (SD-1) |
| the 11 new staging files and every other file | `same` | — |
| check | identical findings; nothing only on head | — |

**`test-omic-enrich-module.R` leaves:**
- W1 run-step `history` → `"<unchanged>"` at `:1053`, `:1485`, `:1865`, `:2119`, `:2366`, `:2640`, `:3469`.
- W1 run-step `states` → `"<unchanged>"` at the following lines.
- W113 `history` → `c("initial","correlation_filtered")` at `:3921`, `:3957`, `:3992`, `:4028`, `:4063`, `:4101`, `:4163`.
- The readers P-28 `:4327`, A107-04 `:4385-4388` and A107-08 `:4428-4432` fail at those same leaves.
- These are `same`: W3 ×2, P-07, A107-04 (no terms), and A107-08's failed-run part.
- The clusterProfiler known failures are unchanged.

**`prot_dia` leaves** (`test-gui-entry-proteomics.R`). For each of the six steps from `enrichment` to `summary_export`:
- `head` → `correlation_filtered`;
- `history` loses its last element;
- `states` loses `enrichment_completed`;
- `…da_status` → "[OK] Analysis Complete …".

Any other difference is stop rule 2 or 3. There is no second correction.

**Expected cost:** about 955 job-minutes (63.7 core-hours), inside R2's 1,080-minute ceiling.
