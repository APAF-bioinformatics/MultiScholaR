# KIT-241 outcome: the enrichment results archive without an external zip program

Ticket r1.1 (`tickets/kit/KIT-241.md`, after AUDIT-KIT-241-plan-20261125). Branch `kit-241` from
`v0.9-staging` 3eeaf100. Commits: 2838e59e (probe script, docs only), 89c3f246 (tests-only), e4b0f98e (fix,
NEWS, ledger row). Run branches (never merged): `run-k241-g0` (693c5d86 = 2838e59e + the probe stage of
`distributable.yml`), `run-k241-w1` (9619cb16 = e4b0f98e + the same stage + the Explorer mechanics repair).
Runs and minutes: `runs.tsv`. No R ran on the laptop.

## Observed

### G0 (run 45577719573, commit under test 3eeaf100; `g0-base-values.tsv`)

The commit check printed `val_product_tree_ok` in both jobs; the installed MultiScholaR names SourceCommit
3eeaf100; closure 385 packages with zip in the plan on both platforms; zip 3.0.1 in every installed library.

| Host | R_ZIPCMD | Sys.which("zip") | Case A (default `utils::zip`) | Case B (`zip::zipr`) |
|---|---|---|---|---|
| Windows, installed (shell's environment) | unset | none | status "ok", warning `'"zip"' not found`, no error, **no archive** | archive, the two expected members, R unzip members equal |
| Windows, case C (folder `Ü-数据`) | | | (not run: case C is case B) | archive, expected members, members equal |
| Linux runner (ubuntu-22.04) | `/usr/bin/zip` | `/usr/bin/zip` | archive written, expected members (754 bytes) | archive, expected members |
| bare `ubuntu:22.04` container (no `/usr/bin/zip`) | `/usr/bin/zip` | none | status "ok", warning `error in running command`, **no archive** | archive, expected members |

Windows' own extractors on G0's case B and C archives: `Expand-Archive` extracted the members equal by SHA-256;
Explorer's engine row was an **instrument fault**: "You cannot call a method on a null-valued expression"
(`Shell.NameSpace` returned null; the archive path R wrote had `/` separators and the engine ran in a background
job; W1 changed both, so which one mattered is not isolated). No result was produced for
that row; G0 was not re-sent for it (W1 extracts the same writer's archives with the repaired mechanics: paths
normalised with `[IO.Path]::GetFullPath`, the engine in an STA Windows PowerShell process). No expected value or
stop rule changed (VAL-54 condition i). `tempdir()` lengths: Windows 47 characters (case folders 61, case C 66);
Linux 15. On Windows `R_RTOOLS45_PATH` names only the missing `no-rtools` folder.

Stop rules: 1 not met (no archive on the installed Windows bundle), 2 not met (case B everywhere, zip present
everywhere), 8 not met (case A failed in the bare container): the NEWS Linux clause stays.

### R1 (run 45583151571, tests-only 89c3f246, no base)

`test-kit241-enrich-download-zip.R`: 9 passed, 3 failed, 3 warnings (utils::zip's). E1, E2 and E3 each fail by
assertion at `file.exists(archive)`: actual FALSE (the value G0 recorded); their status expectations ("ok",
"ok", "error" with "write failed") pass on base; E4 and G1 pass (guards). `r1-tests.tsv`, `r1-failures.txt`.

### R2 (run 45585303533, e4b0f98e against 89c3f246, `check: true`)

`r2-compare.tsv`: the new file `fixed` (head 15 passed, 0 failed); every other file `same`:
`test-omic-enrich-module.R` (1051/19 on both, identical messages), `test-enrichment-significance.R` (191/0),
`test-module-ci-prot-enrich.R` (79/0, 2 skipped), `test-module-ci-manifest.R` (33/0),
`test-e2e-proteomics-enrichment.R` (2 skipped on both: "shinytest2/chromote browser dependencies are not
available"), `test-prot-11-enrichment-module-contracts.R` and `test-prot-11a-…-characterization.R` (0 tests on
both: quarantined Q-163 and Q-164). R CMD check: 0 errors, 3 warnings, 3 notes on both sides; 0 findings only on
head (`r2-check-compare.txt`).

The W1 capture leaf "clusterprofiler download download" in `test-omic-enrich-module.R` fails on base and head with
the same message: the stored summary digest e7a41994… (`test-omic-enrich-module.R:1691`) against the actual
6c76b09e… on **both** sides, and `clusterProfileR_results.tsv` equal to its stored digest 671ae3ed… on both. The
same actual member digests under `utils::zip` (base) and `zip::zipr` (head) are direct evidence that the archive's
members are unchanged by the fix; the stored summary value is pre-existing (it lacks the "Significance:" line the
writer prints; recorded as pre-existing by KIT-226's control) and not this ticket's.

The gprofiler W1 download leaf passes on base (Info-ZIP) and on head (zipr): the summary 13ce7979… and
`gprofiler2_results.tsv` 4182e806… (`test-omic-enrich-module.R:1259-1260`) equal their stored values under both
archivers. Every member of both W1 download captures is therefore observed unchanged by the fix (the summary with
its date line masked).

### W1 (run 45588834064, commit under test e4b0f98e; `w1-head-values.tsv`)

`val_product_tree_ok` in both jobs. The build's package commit and the installed `SourceCommit` are the run-branch head 9619cb16, not e4b0f98e:
`package_commit_rule` (`packaging/package-commit.sh`) names the head when the package source differs from its merge
base with the release line; `val_product_tree_ok` shows 9619cb16's product tree (`PKG_PATHS` and `packaging/`)
equal to e4b0f98e's. Writer default printed as `function(zipfile, files, flags)
zip::zipr(zipfile = zipfile, files = files)`. Closure 385 packages, zip in the plan: equal to G0 (N); the build's plan line reads `closure: 385 packages (roots 137)` at G0 and `(roots 138)` at W1 on
both platforms, every other count equal: one more root, zip, as stop rule 5 predicted.

| Host | Case A (the fix's default) | Windows extractors |
|---|---|---|
| Windows, installed | status "ok", no warning, archive (661 bytes), expected members, R unzip equal | `Expand-Archive` and Explorer's engine: extracted, members equal (SHA-256) for A, B and C |
| Windows, case C (non-ASCII folder) | archive, expected members, equal | as above |
| Linux runner | archive, expected members, equal | — |
| bare `ubuntu:22.04` (no zip) | archive, expected members, equal | — |

Stop rule 5 not met.

## Acceptance criteria

| # | Verdict | Evidence |
|---|---|---|
| A0 | met | G0: Windows case A no archive; case B the two members on every host; case A recorded per Linux host; zip 3.0.1 found everywhere (`g0-base-values.tsv`). Explorer's G0 row an instrument fault, covered by W1 |
| A1 | met | R1: E1 to E3 fail by assertion at FALSE; G1 (and E4) pass |
| A2 | met | R2: new file passes; every other file `same`; check: nothing only on head |
| A3 | met | W1: case A opens on every host; Windows with R's unzip, `Expand-Archive` and Explorer's engine, SHA-256 equal; closure count equal |
| A4 | open | the next checkpoint's full suite (E-3, the per-ticket test rule, `00-KIT-PLAN.md:548`) |
| A5 | met | NEWS (one entry under the unreleased heading, Linux clause kept by G0), ledger row `VAL-<merge>`; `man/`, `NAMESPACE`, `packaging/`, lock file unchanged (`git diff --stat v0.9-staging kit-241`) |

## Inferred, not observed

- What a user sees in the Electron window on the download (F5, handed to KIT-237): G0 and W1 call the writer from
  Rscript under the installed R, not through the window.

## Notes

- N2: `packaging/multischolar.lock.json` keeps `"Declared": "indirect"` for zip; after this change zip is
  declared. Nothing reads the field; the lock is left alone.
- Merge order (coordinator, P03): after KIT-215. KIT-215 was not on staging when the fix was written; it is on
  staging 16112917, which was merged into `kit-241` after the implementation audit (the writer merged without a
  conflict). Neither `test-kit215-not-tested-lists.R` nor KIT-226's `test-kit226-gprofiler-columns.R` was on this
  ticket's base, so neither is in R2; they run in the merge train's run.
- Conventions (I10): the new lines follow the parameter list's alignment and exceed R.md's 100 characters, as 48
  lines of the file do; the tests' injected doubles use `function(...)` as the file's neighbouring tests do.
- Push-triggered Module CI (I12): pushing the two run branches started Module CI runs 45577710878 (166
  job-minutes) and 45585334125 (157), not requested by this ticket; their only failures, `module-ci impacted
  (transcriptomics)` MCI-028/030/031/032, are identical to staging 16112917's run 45586880341. Listed in `runs.tsv`.

## Not checked

- The full suite (A4): at the next checkpoint.
- A source install on Windows with zip 3.0.0 (the laptop's verified library and the lock): such an install may
  write an archive Windows Explorer refuses (r-lib/zip#153, fixed in 3.0.1). Every bundle and runner carries 3.0.1.
  No floor is set in this merge (coordinator's ruling on E-5 and I06: a floor would stop the package loading on a
  library holding 3.0.0); the floor `zip (>= 3.0.1)` is set at the next lock refresh, an item the coordinator owns.
- The browser e2e lane `test-e2e-proteomics-enrichment.R`: skipped by run-tests (no shinytest2/chromote);
  covered by the next stage-`full` run's Linux e2e lanes.
- `test-prot-11` and `test-prot-11a` (G2's injected-fake tests there) are quarantined (Q-163, Q-164) and ran no
  test on either side; G2 is covered by `test-enrichment-significance.R` only.
- Explorer's engine on G0's base archives (instrument fault; W1 covers the fix's archives). The G0 fault's cause is
  not isolated: the archive path had `/` separators and the engine ran in a background job; W1 changed both, so
  which one mattered is not isolated.
- The Explorer step's two weaknesses, for anyone reusing it (I09; KIT-237, KIT-245): with CopyHere flag 1024 (no
  error interface) a refusal by Explorer (the #153 class) would surface as a 120-second timeout, which the step
  classes as an instrument fault rather than a result: class "members missing when the engine finished or timed
  out" as a result ("not extracted: …") and keep the instrument fault for a null namespace or a process that did
  not start. And case C's Windows extractors read an ASCII copy of the archive (`D:\a\_temp\probe\C\download.zip`),
  not the file under the non-ASCII folder.
- Downloads through the Electron window (F5); Push to GitHub without git (F4, KIT-245).

## For the change report

**What was wrong.** The proteomics enrichment tab's "Download All Results" built its archive with R's
`utils::zip`, which runs an external program named by `R_ZIPCMD` (default `zip`). The installed Windows
application gives R a reduced environment with no zip program (PATH = System32, Windows, WindowsPowerShell;
`R_*` removed; Rtools hidden); the Linux bundle's R names `/usr/bin/zip`. Where the program is missing,
`utils::zip` only warns; the writer still returned "ok", no archive was written, and shiny answered the browser's
download with HTTP 404. No message reached the user.

**Who was affected.** Every user of the installed Windows application who pressed "Download All Results" after a
proteomics enrichment (TMT and label-free; DIA is not offered on Windows), and Linux users of the bundle on
computers without `/usr/bin/zip`. The enrichment results themselves, the tab's tables and plots, the Summary and
the publication workbook were unaffected. No numbers were wrong; a feature produced nothing.

**Legacy.** On `main` the same call is `utils::zip(zipfile = file, files = files_to_zip, flags = "-j")`
(`R/mod_prot_enrich.R:2130`, and `:2135`, `:2142`); the default `zipFn = utils::zip` is in tag `v0.5.0-rc1`
(`R/mod_prot_enrich_server_helpers.R:1937`) and tag `v0.9.0` (`R/mod_omic_enrich_panels.R:545`). On `main` and
in v0.9.0 the same call also failed, silently, for package users on Windows without Rtools 4.5 at its default
place (R for Windows puts `${RTOOLS45_HOME:-c:/rtools45}/usr/bin`, where Rtools' `zip.exe` is, at the front of
PATH in its own profile, `src/library/profile/Rprofile.windows:64-86` in R 4.6) and on Linux computers without a
`zip` program; the 0.9 installed application made it fail for every Windows user.

**Code.** Base 3eeaf100: `R/mod_omic_enrich_panels.R:557` (`zipFn = utils::zip`), used at `:614`, `:618`, `:630`.
Fix e4b0f98e: `R/mod_omic_enrich_panels.R:557-560` (`zipFn = \(zipfile, files, flags) zip::zipr(zipfile =
zipfile, files = files)` with its comment), used at `:617`, `:621`, `:633` unchanged; `DESCRIPTION` Imports gains
`zip` (already shipped in both bundles as an openxlsx dependency, 3.0.1; no package added, closure 385 before and
after).

**What changed in results.** Nothing numeric. The archive is now written; its members and their contents are the
same files (R2: every member digest of both W1 download captures the same under both archivers); the archive container's own bytes differ
(compression level, timestamps). On the error path the archive holds only `download_error.txt` (zipr overwrites;
Info-ZIP added to a partial archive). After the fix an archiver failure raises: the writer archives an error
note, and only a second failure leaves the writer (HTTP 500 in a browser).

**Evidence.** G0 run 45577719573 (the failure on the installed Windows bundle and in a bare ubuntu:22.04
container; success with zip::zipr on both and on the runner); R1 run 45583151571 (E1 to E3 failing first); R2 run
45585303533 (regressions base against head, R CMD check); W1 run 45588834064 (the fix on the installed bundles:
Windows with R's unzip, Expand-Archive and Explorer's engine, SHA-256 equal; runner; bare container). Limits:
synthetic one-row tables; Windows Server 2025 runner, not a desktop Windows 11; the Electron window's download
handling not observed (F5).
