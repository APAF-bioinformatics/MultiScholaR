# Janitor Refactor Closeout Plan

## Objective

Close the `janitor` refactor against two independent fidelity baselines:

- `main` is the completeness baseline. Every original export, class, generic,
  method, and effective function must remain accounted for.
- The frozen post-peptide-QC `janitor` checkpoint is the behavioral baseline.
  Structural work must preserve its intentional scientific corrections.

The closeout is complete only when the package has no oversized `.R` files, no
duplicate entity keys, no stale extraction scaffolds, a documented filename
convention, and passing package-loaded isolated contracts.

## Starting Baseline

- `janitor` and `origin/janitor`: `f14e5f8`
- `main` and `origin/main`: `436cdbe`
- `.R` files: `323`
- files over 1,000 LOC: `18`
- files over 2,000 LOC: `6`
- duplicate entity keys: `38` (`21` exact, `17` variant)
- redundant duplicate occurrences: `42`
- stale `TODO: Extract` headers: `18`
- `R/func_general_helpers.R`: `1,789` LOC
- `R/func_peptide_qc_imputation.R`: `1,090` LOC
- `tests/testthat/test-*.R` files: `274`

These figures must be regenerated after every production wave; they are not
permanent expectations.

## Non-Negotiable Rules

1. Move top-level code by exact-source manifest extraction wherever possible.
2. Stage and review every structural wave before applying it to live `R/`.
3. Keep public wrappers and public symbol names stable while files move.
4. Treat `DESCRIPTION` `Collate:` order as executable behavior.
5. Compare effective definitions with the frozen post-peptide checkpoint after
   every wave; compare surface completeness with `main` at closeout.
6. Do not combine structural movement, behavior fixes, and broad naming changes
   in one wave.
7. Do not accumulate extracted code into another generic `*_helpers.R` monolith.
8. Leave unrelated worktree content out of refactor commits.

## Naming Contract

The repository's established `func_` prefix is its `golem::add_fct()`
equivalent. A repository-wide `func_` to `fct_` rename is out of scope because
it adds coupling risk without improving ownership.

Use these runtime filename forms:

- `app_ui.R`, `app_server.R`, `run_app.R`, `launch_app.R`, and `zzz.R` for
  Golem package entrypoints
- `mod_<domain>.R` for a top-level omics workflow shell
- `mod_<domain>_<feature>.R` for a cohesive complete module
- `mod_<domain>_<feature>_ui.R`
- `mod_<domain>_<feature>_server.R`
- `mod_<domain>_<feature>_<responsibility>_helpers.R`
- `func_<domain>_<feature>_<responsibility>.R`
- `func_<domain>_s4_<responsibility>_class.R`
- `func_<domain>_s4_<responsibility>_methods.R`
- `func_<domain>_s4_<responsibility>_helpers.R`
- `utils_<crosscutting_purpose>.R`

All names use lower snake case. Shared S4 generics load first from
`func_general_s4_generics.R`; this is an ownership rule, not a filename
exception.

Within peptide code, `func_pept_s4_*` owns peptide S4 behavior and
`func_prot_qc_peptide_*` owns the proteomics peptide-QC workflow.

## Execution Order

### 1. Freeze The Behavioral Baseline

- Review and commit existing parity remediation separately from audit tooling.
- Exclude unrelated untracked content.
- Run package load, focused remediation tests, recursive function inventory,
  surface inventory, and refactor-tool tests.
- Record the resulting commit as the immutable structural comparison point.

### 2. Assimilate Peptide QC

- Retire the exact duplicate
  `peptideMissingValueImputationLimpa(PeptideQuantitativeData)` definition.
- Keep `R/func_pept_s4_missingness.R` as the S4 method owner.
- Remove the stale extraction inventory from
  `R/func_peptide_qc_imputation.R`.
- Remove confirmed dead audit helpers and correct misleading audit labels in a
  behavior-fix checkpoint before structural movement.
- Split imputation by peptide helpers, protein helpers, and S4 ownership.
- Split replicate filtering by technical-replicate filtering,
  identification-evidence filtering, and confidence/q-value filtering.
- Keep the audit and exclusion implementations intact unless measured
  duplication or responsibility analysis justifies a split.

### 3. Canonicalize Duplicate Entities

- Resolve all `17` variant duplicate keys first.
- Preserve the frozen checkpoint's effective winner unless a separately tested
  behavior fix explicitly supersedes it.
- Remove all `21` exact duplicate keys by ownership family.
- Re-run entity multiplicity and effective-winner checks after each family.
- Exit with zero duplicate entity keys.

### 4. Remove Extraction Scaffolding

- Remove all `18` stale extraction headers.
- Do not implement historical aspiration names that never existed.
- Replace only useful historical text with short ownership notes.
- Remove empty migration breadcrumbs after updating filename references.

### 5. Split Files Over 2,000 LOC

Split by stable responsibility, in this order:

1. Proteomics enrichment module helpers: state, observers, execution, rendering,
   and downloads.
2. General enrichment: preparation, engine adapters, result shaping, plotting,
   and persistence.
3. Multiomics enrichment: preparation, orchestration, pathway engines, and
   result integration.
4. Metabolomics normalization module: state, workflow, QC, rendering, and
   observer registration.
5. Metabolomics S4 normalization: transformation, normalization, RUV, and
   design cleanup methods.
6. File-management config: reading, mutation, formatting, serialization, and
   workflow construction.

### 6. Split Remaining Files Over 1,000 LOC

- Work by tested domain: proteomics, lipidomics, metabolomics, then shared code.
- Reduce `R/func_general_helpers.R` below 500 LOC or retire it as an umbrella.
- Keep normal targets in the 150-500 LOC band.
- Permit 501-800 LOC only for cohesive orchestrators or S4 method families.
- Do not leave any runtime `.R` file over 1,000 LOC.

### 7. Filename Closeout

- Apply the naming contract after function ownership is stable.
- Update `DESCRIPTION`, roxygen source references, tests, dev scripts, reports,
  workbooks, tickets, and refactor documentation.
- Remove empty nonconforming runtime files instead of renaming breadcrumbs.
- Do not rename exported function symbols merely to match filenames.

### 8. Final Fidelity Gate

- Zero parse failures and zero package-load-order warnings.
- Every `main` export, class, generic, and method remains present.
- No unexplained effective-definition drift from the frozen checkpoint.
- Zero duplicate entity keys.
- Zero stale extraction headers.
- Zero runtime `.R` files over 1,000 LOC.
- Focused tests pass after each wave.
- All package-loaded isolated contract files pass.
- Package check introduces no new error, warning, or note.
- Regenerate the file-size, filename-coupling, recursive-function, surface, and
  final parity reports.

## Wave Gate

Every production wave follows the same sequence:

1. Audit filename coupling and identify the focused test surface.
2. Add or confirm characterization coverage.
3. Author and verify an exact-source manifest.
4. Generate and inspect the staged targets.
5. Apply transactionally and update `Collate:`.
6. Run post-apply selector, parse, multiplicity, and effective-winner checks.
7. Run focused package-loaded tests.
8. Commit at the responsibility or module boundary.
9. Update this plan and the current audit counts.

## Progress

- [x] Freeze and checkpoint the post-peptide behavioral baseline at `02d596c`.
- [x] Assimilate peptide-QC code into the refactored ownership model.
- [x] Reduce duplicate entity keys from 38 to 0.
- [x] Remove all 18 stale extraction headers.
- [ ] Split all six files over 2,000 LOC.
- [ ] Split all remaining files over 1,000 LOC.
- [ ] Complete filename normalization.
- [ ] Pass the final dual-baseline fidelity gate.

Current checkpoint after extraction-scaffold cleanup:

- runtime `.R` files: `329`
- files over 1,000 LOC: `15`
- files over 2,000 LOC: `5`
- duplicate entity keys: `0`
- redundant duplicate occurrences: `0`
- stale `TODO: Extract` headers: `0`
- `R/func_general_helpers.R`: `1,489` LOC
- imputation ownership wave: `dc5d6ec`
- replicate/evidence/confidence ownership wave: `fa5ed55`
- duplicate canonicalization wave 1: `67c44b9`
- duplicate canonicalization wave 2: `a1fb983`
- duplicate canonicalization wave 3: `58a3bef`
- extraction-scaffold cleanup: `66001da`
- proteomics enrichment responsibility split: this checkpoint

Wave 3 preserved all `3,339` recursive function occurrences and all `2,875`
unique function ASTs from `a1fb983` exactly, with zero unmatched expressions or
AST drift. It replaces the remaining load-order-dependent cross-omics methods
with one tested list router per generic and retains every prior domain method as
an internal helper. The package inventory now contains one `Collate` entry for
every runtime file and zero duplicate entity keys.

The extraction-scaffold cleanup removed the remaining `17` stale inventory
headers and retired five comment-only runtime breadcrumbs. All `3,345` function
occurrences from `58a3bef` remain AST-identical, the runtime and `Collate`
inventories both contain `322` files, and `75` generated Rd provenance lines
now identify their actual roxygen-owning source files.

The proteomics enrichment module split replaces the `4,224`-line accumulated
helper file with eight exact-source owners between `258` and `754` LOC. Its
focused package-loaded gate passes `1,412` expectations, and all `3,345`
recursive function ASTs from `66001da` remain exact.

The isolated `test-prot-s4-missingness-direct-shared.R` file currently has two
errors in untouched protein missingness paths. Both errors reproduce against a
clean archive of pre-structure commit `7e74163`; they are inherited defects, not
peptide extraction regressions, and remain mandatory remediation before the
final all-contract gate.

The isolated `test-prot-s4-qc-direct-shared.R` file also has one stale
error-message expectation in the protein intensity-filtering path. Its other 71
expectations pass, and the same assertion fails identically against a clean
archive of pre-wave commit `fa5ed55`; it is likewise inherited final-gate debt.
