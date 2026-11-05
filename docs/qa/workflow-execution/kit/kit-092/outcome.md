# KIT-092 outcome: tables inside products fingerprinted by content (protocol v1.7)

**The verifier now judges a table nested in a product by its content, as it already judged a state's own tables. It also no longer passes a step without an oracle, or an in-memory member nobody declared. Re-verified offline under v1.7, the two partitioned journeys `089u` and `089v` are accepted and exact. The four KIT-086 members are still accepted and exact. No journey ran.**

Authority: the user, 2026-11-06. "go with 1" answered KIT-089's question. "fix all recommendations" was given on the transcriptomics re-cut audit, finding KIT098-F02. Plan: `plan.md`. Protocol: `protocol.md` section 9.

## The change

| Where | Before (v1.6) | Now (v1.7) |
|---|---|---|
| `kitCanonical` | A nested `parquet_duckdb` reference was its file hash (`payload$sha256`). | With a store, it is the table's pages (`kitFingerprintTable`), plus its row count and column names. Without a store (the args summary, `.rds` files) it is unchanged. |
| Journey verifier | Recomputed the state's own products. The products the driver attached to each step stayed as recorded in-run. | Recomputes both from the store (`kitStepProducts`, the driver's own attachment rule). |
| `kitCases`, single journey | A step without a bound oracle was skipped at verification (`kitReproduceMismatches`). | The member is refused at binding if any step it runs that leaves a state lacks an oracle. A member with no oracle source must declare `in_memory`. |
| In memory | Inferred: a route other than DIA, bound without oracles. Exactness `NA`, populations and head rules waived. | Only where the binding declares `in_memory = TRUE`. The KIT-054 legacy members `054l` and `054m` declare it. |

## Re-verification (`reverify.txt`)

| Member | Build | Recorded decision (v1.6) | v1.7 decision | Oracle mismatches |
|---|---|---|---|---|
| `089u` lipid screen, missing values, `partitioned` | `kit_089` | invalid (`verifier_accepted`) | **accepted, exact** | 0 |
| `089v` metab screen, missing values, `partitioned` | `kit_089` | invalid (`verifier_accepted`) | **accepted, exact** | 0 |
| `086u` lipid screen, missing values | `kit_087` | accepted, exact | accepted, exact | 0 |
| `086v` metab screen, missing values | `kit_087` | accepted, exact | accepted, exact | 0 |
| `086p` metab large, missing values | `kit_087` | accepted, exact | accepted, exact | 0 |
| `086q` lipid large, missing values | `kit_087` | accepted, exact | accepted, exact | 0 |

**How each member was re-verified.**
- Its stage directory was reflinked into `kit/scratch-092-reverify`. The runtime directory was only read.
- Its recorded binding was checked twice, and it equals its case both times:
  - under the bindings it was launched with;
  - under the current bindings, with the two inputs v1.7 itself changed taken from the launch: the protocol hash and the KIT-005 canonical-form record.
- `verifyKitMember` ran on the member's own build, with the current harness and oracles.
- `kitSingleJourneyDecision` re-decided the member on its recorded receipt.

The decisions recorded under v1.6 stand as recorded. These are re-verifications under v1.7, not new members.

## Oracles (`refingerprint.txt`, `kit005-check.txt`)

- **`lipid86s`, `metab86s`, `metab86`, `lipid86`** (`kit-086/oracles-*.rds`): re-fingerprinted from their stores. On each:
  - the store's head is the recorded head;
  - the v1.6 recomputation equals the record exactly;
  - under v1.7, class, rows, components, support and summary are unchanged;
  - 47 product entries changed name or hash. These are the `_science` identity products and the `<type>:<sha12>` result products, renamed by their new hash.

  The v1.6 files stay in git history.
- **KIT-005 journey oracles:** the bindings generator has recomputed them from control-1's store at every bind since KIT-012, so they now follow v1.7.
  - Two entries changed: `journey_design`'s `peptide_data` component, a record that nests tables, and the identity product `normalization_science`.
  - The products the driver attached reproduce unchanged under v1.7 (`kit005-check.txt`).
- **Every other bound oracle** keeps its v1.6 fingerprints and needs a re-fingerprint before any new comparison (protocol 9.3). KIT-088's baselines don't exist yet and are recorded under v1.7.
- **Unavailable oracles:** none. All four stores exist at their recorded heads.

**The bindings** were regenerated with the build list and flags of the last bind (`scripts-086/bind.sh`). The diff is limited to:
- the protocol hash and the two harness file hashes;
- the four `kit-086` oracle prefixes, their products only;
- the two KIT-005 entries above and the KIT-005 canonical-form record;
- `054l`/`054m` `in_memory`.

## Checks

- **Test first (`test-first.txt`):** three new tests fail before the change, 11 failures and 1 error. After: `test-workflow-kit-harness.R` 551 pass, 0 fail.
  - The tests: the same content in different bytes matches, and a changed value inside a nested table mismatches; attached products are rebuilt; an unbound step oracle, a member without an oracle source and an undeclared in-memory member are refused, and an undeclared lifecycle member meets the populations and head rules.
  - The existing single-journey tests now bind synthetic oracles. The KIT-054 legacy lane and the KIT-068 lifecycle test declare `in_memory`.
- **Dependents** (the tests that source the harness): `test-workflow-analysis-record.R` 239 pass, `test-workflow-lipid-declarations.R` 258, `test-workflow-metab-declarations.R` 258, 0 fail.
- **Attempts:**
  - The first bind failed at once: `perf-107/synthetic-corpus.json` is git-ignored and was absent from the worktree, so it was copied in. The second bind completed.
  - The first re-verification of `089u` stopped on the binding check. The case digest carries the KIT-005 canonical-form record, which v1.7 changed. `reverify.R` now takes that record from the launch bindings, as it does the protocol hash. This was the one correction pass.
- **Third bind:** the KIT-005 finding was added to protocol 9.3, so the bindings were regenerated once more. They differ from the second bind in the protocol hash only. The re-verification ran on the second bind, and the protocol hash is the one input it takes from the launch anyway.
- **Not run:** `R CMD check` and the full suite. No package code changed: only `tools/profiling`, a test file and `docs/`.

## Recorded, not changed

- A DIA single journey on a corpus other than `full` and without oracles (`028`, the replay origin) would now be refused at binding. It has no oracle source, and the verifier never judged it against one. Its recorded decision stands.
- Opening a store to read it appends to its `checkpoint-io.log`, so the four oracle stores' logs grew during the re-fingerprint. Every bind already does the same to the KIT-005 control-1 store. Heads and entities are unchanged. The six member directories were not touched: the verifier read reflinked copies.
- `kitGuiMemberDecision` still recomputes only the state's own products of both members, so the attached products stay as recorded. These DIA `gui` members are outside this ticket.
