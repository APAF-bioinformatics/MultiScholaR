# KIT-174 feasibility (F2 for part (a); F1 and F2's part (b) items when part (b) starts)

Read at `v0.9-staging` 0cf0e440 (branch point of `kit-174`). Nothing was run for these answers; "read" means read at
the place cited.

## F2, part (a) items (answered before part (a)'s code, A1)

| Question | Answer (read) |
|---|---|
| Tests or stored references that pin the proteomics RUV grouping help | None. `grep -rn "Experimental factor to preserve\|select 'batch'" tests/ inst/` finds only `tests/testthat/test-omic-gui-config.R:1958`, a synthetic configuration group ("Experimental factor to preserve.") used to test the validator, not the shipped string. `docs/methods-decision-ledger.md` quotes the text (DAS-003, UVR-021); not a test |
| Stored references that pin the proteomics normalisation panel | `test-omic-norm-module.R`'s frozen UI digests cover lipid, metab and trx only (`omicNormFrozenUi`, `:107-182`; the test at `:184`); the proteomics frame tests (`:1896-2070`) assert ids and outputs, not the help text. The gate-off script hashes `app_ui(NULL)`, which does not render the omic tabs (ticket section 7) |
| Does any stored lane choose a grouping column other than `group` (so the warning would fire)? | No. The GUI-entry lanes never set `ruv_grouping_variable` (`tests/testthat/helper-gui-entry.R:520`, `:615`, `:714` set only the mode); the five proteomics lanes record "RUV grouping variable: group" (`tests/testdata/gui-entry/prot_*.json`); the lipid and metabolomics lanes record no grouping. Every test that calls the three run actions passes `group` (`test-omic-norm-ruv-design.R:899`, `:916`; `test-module-ci-prot-dia-limpa.R:161`; `test-prot-design-log2.R:532`; `test-workflow-saved-qc.R:401`; `test-prot-05b-norm-module-contracts.R:1567`; `test-workflow-gui-memory-lifetimes.R:684`; the shared frame through `test-omics-norm-state-guards.R`'s and `test-kit207-summary-head-guard.R`'s inputs, `ruv_grouping_variable = "group"`). Tests that use `batch`, `tech_rep_group` or other columns call the kernels directly, not a run action |
| Does any file under `tests/testdata/gui-entry/` or any value G1(a) compares change? | No (inferred from the two rows above: the help text is not in any lane's captures, and no lane's run reaches the warning). The stop rule does not apply |

Part (b) items of F2 (import notifications, `initial` state value, report paragraphs, DIA-NN filter sentences, the TMT
`Master.Protein.Accessions` question) and F1 are answered when part (b) starts, at its own branch point.
