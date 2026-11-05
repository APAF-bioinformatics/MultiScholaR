# KIT-092 plan: content fingerprints for tables inside products (protocol v1.7)

Authority: the user, 2026-11-06, "go with 1" on KIT-089's question, then "fix all recommendations" on the transcriptomics re-cut audit (KIT098-F02). Ticket revision 1.1.

**Bounded outcome.** The verifier judges a table nested in a product by its content, as it already judges a state's own tables. A partitioned journey that is content-identical to a whole one is accepted. A step with no bound oracle, or an undeclared in-memory member, can no longer pass. Nothing else in the protocol changes.

## File by file

1. **`tools/profiling/workflow_kit.R`**
   - **`kitCanonical`**: when a store and an api are given, a nested `parquet_duckdb` reference becomes `list(.ref = "table", rows, columns, content = kitFingerprintTable(...)$sha256)`, memoised per entity id. This is the form `kit-089/content-compare.R` used, plus the row count and column names, so a zero-row table still carries its schema. Without a store, as for the args summary and the `.rds` files, the payload hash stays. That keeps the summary hashes the KIT-089 comparison found identical.
   - **The verifier's overlay** (`verifyKitMember`, the oracle branch): it recomputes the products the driver attached to each step from the store. A new helper, `kitStepProducts`, is the driver's and recorder's attachment rule (first-appearance step, ordered by hash, named `<type>:<sha12>`). Before, only the state record's own products were recomputed.
   - **Exactness gaps (KIT098-F02)**, in `kitCases`:
     - A single journey member is refused at binding when an oracle is missing for any step it runs that leaves a state (`initial` excepted, which carries no data). The test is the one `kitOracleJourneyPopulations` uses.
     - It is also refused when it has no oracle source and does not declare `in_memory`. An oracle source is bound `oracles`, or the KIT-005 oracles for DIA on the full corpus.
     - "In memory" is `isTRUE(case$in_memory)`, taken from the ticket binding or `lifecycle_omics`, in `kitSingleJourneyDecision`, `kitLifecycleDecision` and the verifier's failure branch. It is no longer inferred from `is.null(oracles)` on a route other than DIA. `in_memory` enters the binding digest only when TRUE, so every other digest is unchanged apart from the protocol hash.
2. **`tests/testthat/test-workflow-kit-harness.R`**, written first and failing before the change:
   - Two stores hold the same table content written in different bytes (one table, then the same rows as two row groups). Products nesting each fingerprint identically. A changed value inside the nested table mismatches.
   - A single journey with an unbound step oracle is refused. A non-DIA member without oracles or an `in_memory` declaration is refused. An undeclared omics lifecycle member is judged on populations and the head rule.
   - The existing tests that bound single journeys without full oracles, or relied on the inference, are given synthetic oracles or an explicit `in_memory`.
3. **`docs/qa/workflow-execution/kit/protocol.md`**: a v1.7 section that states both changes and the user's decision, and a version line.
4. **`docs/qa/workflow-execution/kit/bindings-generator.R`**: `054l` and `054m` (the KIT-054 legacy lane, no owner) declare `in_memory = TRUE`. The four `kit-086` oracle files are bound as before, now re-fingerprinted.
5. **`docs/qa/workflow-execution/kit/bindings.json`**: regenerated from the build list, with the environment flags of the last bind (`scripts-086/bind.sh`). The expected diff is the protocol hash, the harness hashes, the `lipid86s`/`metab86s`/`metab86`/`lipid86` oracle products, the `054l`/`054m` declaration and the date.
6. **`docs/qa/workflow-execution/kit/kit-092/`**:
   - `refingerprint.R`: recomputes each oracle from its store under the new form. It checks the head and that class, rows, components, support and summary are unchanged, rewrites the `kit-086/oracles-*.rds` file (the old one stays in git history), and records which products changed in `refingerprint.txt`.
   - `reverify.R`: the six re-verifications, offline, on reflinked copies of the member directories (the runtime directories are only read). Each member is verified with its own build's library and judged by `kitSingleJourneyDecision` on its recorded receipt. Results go to `reverify.txt`.
   - `outcome.md`.

## Oracles re-fingerprinted

- `lipid86s`, `metab86s`, `metab86` and `lipid86`: the `kit_086` baselines that KIT-088 and the KIT-089 journeys compare against.
- KIT-088's other baselines don't exist yet. They are recorded under v1.7.
- Every other bound oracle keeps its v1.6 fingerprints. Its recorded decisions stand. A new comparison against it needs a re-fingerprint first (protocol v1.7).

## Checks and budget

- The harness test file, one process. Its dependents by grep (`test-workflow-lipid-declarations.R`, `test-workflow-metab-declarations.R` and `test-workflow-analysis-record.R` source the harness).
- The bindings regenerate.
- The six re-verifications, one bounded correction pass.
- No journeys and no members.

**Backward compatibility.** Decisions already recorded stand under their own protocol hash. Rollback is `git revert`, and the oracle files keep their v1.6 fingerprints in history.
