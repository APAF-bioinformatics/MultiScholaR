# KIT-108 run 9: the skip comparison

The set is run 1's 31 files, in both load modes. The baseline is KIT-105's full-suite rows overlaid with KIT-110's, KIT-106's and KIT-107's, compared with `kit-105/skipcmp.py`. Files outside the set are filtered out. The full comparison is in [skip-compare.tsv](skip-compare.tsv).

- **`load_all`** ([skips-load_all.tsv](skips-load_all.tsv)) is run 1, under strict skips: 9,034 passed, 0 failed, 0 error blocks.
- **Check mode** ([skips-check.tsv](skips-check.tsv), script [skips-check.sh](skips-check.sh)) ran on the `kit_108` check build (`X6BzNV`, `7b52f535`): 7,226 passed, 0 failed, 49 error blocks. Against the baseline, failures fall from 149 to 0 and error blocks from 391 to 49.

**Every change against the baseline is one of these:**
- **A recorded quarantine:**
  - whole files with a line-1 skip: Q-170 (`test-lipid-14a`), Q-171 (`test-lipid-15`), Q-172 (`test-metab-05`), Q-180 (`test-prot-12`), Q-181 (`test-prot-12a`);
  - blocks: Q-173 (MCI-023, lipid summary), Q-174 (MCI-017, metab summary), Q-182 (MCI-011, prot summary);
  - Q-175 and Q-183: the comparison sides in the new core test.
- **A new file:** `test-omic-summary-module.R` (1,703/0 under `load_all`; 1,454/0 and 4 skips under check) and `test-omics-summary-core.R` (665/0, 2 skips under `load_all`).
- **Check mode only:**
  - **The "KIT-040 references are not shipped" skips** in the two new files (4 in the module test, 2 in the core test): the trx entries read the KIT-040 reference data, which the check layout does not ship. The existing trx files in this set skip the same way (`test-module-ci-trx-summary.R` 4, `test-workflow-trx-declarations.R` 6, and others), and all of them run under `load_all`.
  - **`test-omics-summary-core.R`: two error blocks at `7b52f535`.**
    1. "(W4, F104-01) the native export's metadata …": "cannot open the connection". The fixture helper `readWorkflowCheckpointRecipe` reads `tests/testdata/workflow-checkpoints/scenarios.json` through `here::here()`, and the check layout does not ship it. This is the known source-tree read (as KIT-107's). It passes under `load_all`.
    2. "(W5, D8) the copy's DA workbook is built one sheet at a time …": **a defect in KIT-108's own test.** The callr child chose `load_all` when the namespace path had `DESCRIPTION` and `R/`. An installed package has both, so under check the child called `load_all` on the installed copy and failed on its DLL (`getDLLRegisteredRoutines.DLLInfo`). **Fixed** in the test only: the source tree is the one without `Meta/`. This was the run's one correction pass, and only this file was rerun ([w6/d8-test-fix/](w6/d8-test-fix/)):
       - **check, on the `kit_108` check build** with the fixed file: 537 passed, 0 failed, 1 error block (item 1), 4 skips;
       - **`load_all`, in run 1's worktree and environment:** 665 passed, 0 failed, 0 errors, 2 skips, the same as run 1.

**Two instrumentation slips in the rerun, recorded:**
- **The first `load_all` rerun** used the default library, which has no RUVSeq: 575/13, 4 skips ("{RUVSeq} is not installed"; `load_all-default-lib-rows.tsv`).
- **The second** used run 1's library but the main checkout rather than run 1's worktree: 649/16. The 16 failures are all one field, `study_parameters`, one line short. `resolveWorkflowGitInfo()` looks for a `.git` directory in the working directory. The main checkout has one; a worktree has a `.git` file. So the "Version Information" block differs by one line. W1's captures were frozen in a worktree. This is **F108-04**, handed to KIT-109 (its full suite runs these files).

**Decision:** no file has a new file-level skip, more hidden blocks, fewer tests run or a new skip category without a Q-entry. The one new check-mode error that was not the known source-tree read was KIT-108's own test defect, fixed and rerun.
