# KIT-176 feasibility (A1): F3, the package paths, and the two meanings of `package_commit`

Read 2026-11-23 at `v0.9-staging` b052b824, by reading only; nothing was run. Written down while R1 and R2 were
running, after the fix commit 5700c8a1 had implemented D4 (the ticket's section 8 asks for the answer before
implementation). The answer below did not change with the runs, and its run-time leg (A5: no difference family in
one run only, none involving the software block) held in R2.

## F3: does any reference the bundle lane compares hold the record's `source_commit`?

**No.** D4 (VAL-31) is implemented; A5's family comparison in R2 is the run-time check (F3's extended stop rule).

The bundle lane compares two kinds of reference: the GUI-entry lanes (`packaging/verify/run-lanes.R` against
`tests/testdata/gui-entry/<lane>.json`) and the browser e2e lanes (`packaging/verify/run-e2e.R`, the
`tests/testthat/test-e2e-*.R` files as test counts). The Windows figure-path cases run named tests.

| Reader of `source_commit` | What it does | Compared by the bundle lane? |
|---|---|---|
| `R/utils_workflow_analysis_record.R:457`, `:472` | backfilled steps get `NA`; each step's `software$sha` | only inside `analysis_record.json/.xlsx/.docx`, which the GUI-entry recorder lists **by path only** (`tests/testthat/helper-gui-entry.R:18-21`: "sha256 is taken for .tsv/.tab/.csv/.fasta artifacts only (records, logs, ... and reports ... are recorded by path)") |
| `:564` | the record's top-level `software$sha` | the same |
| `:1039` | the record-to-report field map (`source_commit` -> `software`, `sha`) | reports (`inst/reports/proteomics/*_report.rmd`, "Source commit") are recorded by path |
| `R/utils_workflow_dia_state.R:808`, `:819` | the provenance entry and the software environment of a commit | provenance is not a state's content; state digests are taken over `getState` values |
| `tests/testthat/test-workflow-dia-products.R:580` | `expect_true(is.character(...) && length(...) == 1L)` | run in the Windows figure cases (`dia_products_old_layout`); a 40-hex value satisfies it as `NA_character_` did; it is not a stored reference |
| `tests/testthat/test-prot-norm-core.R:432` (also `test-omic-enrich-module.R:3885`, `test-workflow-saved-da.R:627`) | exclude `source_commit` from their comparisons | not run by the bundle lane |
| `tests/testthat/test-workflow-analysis-record.R:452`, `:589` | an overlay set to `"abc123"` in the test | not run by the bundle lane |

Searches: `grep -rln "source_commit\|SourceCommit" R inst tests` (the files above);
`grep -l "source_commit\|Source commit\|SourceCommit" tests/testdata/gui-entry/*` (no hit); no `.tsv`/`.csv` artifact
of the analysis record exists (`R/utils_workflow_analysis_record.R` writes none); the e2e test files do not name it.
Output text in the GUI-entry records is masked for hexadecimal runs of 10 or more characters (`helper-gui-entry.R:18`),
so a commit shown in an output would appear masked; no oracle holds "Source commit" text. `build` (sha256 of
`R/MultiScholaR.rdb`) is not changed by writing a DESCRIPTION field.

## The package paths (`PKG_PATHS`)

`DESCRIPTION NAMESPACE R inst src man data configure configure.win configure.ucrt demo exec .Rinstignore`: what
`R CMD INSTALL` of the checkout (`packaging/r/install-library.R`, `rInstall(checkout)`) reads for code and data, from
R 4.6.1's installer as cited in the ticket's section 2.2 (`src/library/tools/R/install.R`: `configure` `:1305-1312`;
Windows `configure.ucrt`/`configure.win` `:1276-1298`; `:1343-1348`; `demo/`, `exec/` `:1705-1722`; `inst/` with
`.Rinstignore` `:1739-1745`; `tests/` only with `--install-tests` `:1781-1787`, which the bundle does not use).
`data`, `configure.win`, `configure.ucrt`, `demo`, `exec` and `.Rinstignore` are not tracked today; each is listed so
that it counts if it ever is (a pathspec that matches nothing is not an error to git).

Left out on purpose: `NEWS.md`, `LICENSE` and `README.md`. R copies them into the installed package, but they are text
only and change no computation; a NEWS-only commit (as this ticket's own) therefore keeps the staging commit as the
package commit. `tests/` is not installed. `packaging/` is not part of the package; what it changes in the installer is
named by the build commit, the lockfile hash and the snapshot, which About keeps.

## Two meanings of `package_commit`

- Here (the bundle lane, `VERSION.json`, About): the commit whose package source the installer holds, by the rule of
  `packaging/package-commit.sh` (the merge base with the release line when the package source there equals the
  build's, otherwise the build commit with `package_on_release_line: false`).
- In the `run-tests` lane, `package_commit` is a request field meaning "the package under test" (any commit of any
  branch; `run-tests.yml:11-13`, `:97-98` on `origin/run-tests`). The two are different things and are not compared.
