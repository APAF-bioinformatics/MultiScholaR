# KIT-134 B-01: Module CI on a Linux runner, measured against the laptop (2026-11-17)

Audit finding B-01 (plan 4.7 audit) has three parts: the workflow had never passed its first job; every stored reference
is bit-exact from the laptop (`-march=native -O3`); CI resolved without `fastduckdb` and `fastlimpadpc`. Its closure
criterion is one full run on a Linux runner, read, and the comparison rule written. The rule is KIT-134 revision 2.0.4,
fixed before any run was read and not changed here: doubles agree when `|a - b| <= 1e-6 * max(|a|, |b|)` or
`|a - b| <= 1e-9`; discrete outputs, outputs whose contract is bit-exactness, and file hashes of copied inputs exactly.

## Verdict: B-01 closes

On run [44054974678](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44054974678), with the library the
bundle ships (the KIT-134 snapshot) on ubuntu-22.04:

1. Every job installs, both internal packages included, and reaches its tests. Four of the six test jobs pass.
2. Every lane runs. DIA, DIA_LIMPA and enrichment pass, as do lipidomics, metabolomics, TMT/LFQ and transcriptomics.
   The only differences from the laptop are the transcriptomics digest oracles. Their values were saved on the runner
   and compared under the rule: every double agrees and every discrete output is identical.
3. Two failures were the workflow's own: a test asserting the workflow's text. They are repaired in `79a03479` and
   verified on the laptop.

No difference is outside the rule and none is unexplained. The laptop's exact comparison stays the gate for every change
to the package. This rule applies to the bundle against the laptop, as 2.0.4 states.

## Workflow (branch `ci-token`)

| Commit | Change |
|---|---|
| `79bc4410` | R 4.6.1 pinned (IR-M3). Cache off and the read token as step env only. The `=?ignore` lines of `5534aa97` removed; the internal packages come from DESCRIPTION's named Remotes (`3ec13e95`). |
| `d2f02aa1` | Audit F2: without the secret, the two packages are ignored. Superseded by `4f0ec49d`'s empty-token branch. |
| `4f0ec49d` | Decision (coordinator, 2026-11-17): CI installs what the bundle ships. This commit makes four changes, listed below. |
| `79a03479` | Each job declares its test tools in `module-ci.yml` again (`extra-packages`). This repairs the two release-gate failures below. |

The four changes in `4f0ec49d`:

- `.github/actions/snapshot-library` (composite) does all the installing, for all 7 jobs that install, now on
  ubuntu-22.04 (the bundle's build platform, D1). It installs, in order:
  - R 4.6.1 from Posit r-builds;
  - the library through `packaging/r/install-library.R`: snapshot 2026-07-30, Bioconductor 3.23, jammy binaries;
  - the internal packages (`MSR_INTERNAL=required` with the token; left out when the secret is empty);
  - MultiScholaR from the checkout;
  - the test tools from the same snapshot.
- The installer and its inputs are copied unchanged from `origin/kit-134` at `ed5c9997`.
- The E2E artifact directories are absolute.
- `.github/scripts/save-trx-products.R` saves the values behind the transcriptomics digest oracles, and
  `packaging/verify/compare-results.R` applies the rule.

Two checks ran on the laptop first (nice 19, the verified library):

- The script reproduces all five KIT-040 digests.
- The laptop's own products agree with the reference: 15 of 15 leaves, maximum relative difference 0.

**No internal package in a cache or an artifact.**

- No `actions/cache` and no setup cache.
- `install-library.R` keeps the internal sources and build logs under `RUNNER_TEMP`; its evidence directory is
  `RUNNER_TEMP/library-evidence`; the library is `/opt/R/4.6.1/lib/R/library`.
- Every upload path is under `tests/testthat/_module_ci_artifacts` or `_e2e_artifacts`.
- Checked after run 44054974678: its 6 artifacts hold 33 files. These are report `.rmd` files, a `.docx`, a DA CSV,
  enrichment zips, the scorecards and the saved transcriptomics values. No file names or contains either package; the
  zip and docx members and the decompressed `.rds` files were checked as well.

## Runs

| Run | Ref | Outcome | Class |
|---|---|---|---|
| [44013986162](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44013986162) | `1e4eb7a7` | every install: "Can't find package called fastduckdb, fastlimpadpc" (pak names an unnamed GitHub Remote after the repository) | (a); fixed by named Remotes |
| [44014345089](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44014345089) | `6767c8ce` | every install: `fastlimpadpc` d4f9cc1d "lazy loading failed" (no `Imports: limpa`) | (a); fixed upstream, 87a90b89 |
| [44040978386](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44040978386) | `8a89123d` | every install: `fastduckdb` 2694b972 `cargo --offline` without vendored crates | (a); fixed upstream, e724f90b |
| [44041893594](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44041893594) | `79bc4410` | first measurement, with current repository versions; below | |
| **[44054974678](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44054974678)** | `4f0ec49d` | **the measurement with the bundle's versions**; below | |

The push-triggered duplicates were cancelled; each was a subset of the dispatched `gate = all` run. They are 44013980819,
44014335748, 44040970699, 44041893072 and 44054974284. A dispatch has no before-SHA, so detect impact routes `none`, and
the push-gate lanes are skipped. The full gates run every test file those lanes select.

## Run 44054974678: the bundle's library (ubuntu-22.04.5, R 4.6.1, snapshot 2026-07-30)

| Job | Result | Install | Wall |
|---|---|---|---|
| detect impact | pass | 8.5 min | 9 min |
| full gate (cross-omic-coexistence) | pass | 10.0 min | 10 min |
| full gate (alternate-launch-browser-smoke) | pass | 8.0 min | 9 min |
| full gate (report-export-fidelity) | pass (failed in 44041893594) | 10.1 min | 20 min |
| full gate (all-module-unit-contracts) | fail: 10 | 9.2 min | 22 min |
| full gate (all-current-e2e-workflows) | fail: 1 | 9.9 min | 40 min |

The library installed as follows:

- `fastduckdb` 0.2.0 (e724f90b, 161 s) and `fastlimpadpc` 0.1.0 (87a90b89, 19 s).
- The DIA manifest was written by configure: supported, native query bound, duckdb 1.5.5.
- The packages on the transcriptomics and DIA paths are at the lockfile's versions: duckdb 1.5.5, edgeR 4.10.1,
  limma 3.68.4, limpa 1.4.2, ruv 0.9.7.1, RUVSeq 1.46.0, statmod 1.5.2.
- The bundle still differs from the lock in 24 packages, none of them on these paths. The 8 in the hard closure are
  arrow, future, gh, httr2, Matrix, plotly, RcppParallel and zip (KIT-134 version table, run 44044598918).

Every E2E file ran without a skip. The unit-contract skips are the declared quarantines (Q-22, Q-23, Q-76, Q-173, Q-174,
Q-182) and the audit-candidate guards. testthat's "maximum number of 10 failures" notice caps only the printed details:
every file ran, and the per-file marks show exactly the 10 failures below.

### Every failure, classified (2.0.4 classes)

| # | Test | What differs | Class | Under the rule |
|---|---|---|---|---|
| 1-2 | `module-ci-release-gate:48`, `:49` | "Missing workflow text: any::shinytest2 / any::chromote": `4f0ec49d` moved the jobs' tool lists into the action | **(a) the workflow's own** | Repaired in `79a03479`; the file passes on the laptop (67 expectations). Not re-run on CI. |
| 3, 5, 7, 8 | `trx-da:7`, `trx-norm:13`, `trx-qc:64`, `trx-summary:40` | sha256 of the regenerated corpus `manifest.json` | **Agree** (b, within the rule) | `counts.tsv`, `design_matrix.tab`, `design.tsv` and `contrast_strings.tab` are byte-identical to the laptop's. Only `truth.json` differs (303,811 against 303,810 bytes): 245 of its 5,400 simulated batch loadings differ in the last bits, maximum relative 9.4e-15 (`rnorm` on two machines). No product reads `truth.json`. The manifest records `truth.json`'s hash, so the manifest's hash differs. |
| 11 | e2e `trx-journey:81` | the same manifest hash | **Agree** | as 3; all other expectations of the journey pass |
| 4 | `trx-da:28` | `topTags` table of the first contrast | **Agree** | Gene order, row names and the FDR < 0.05 calls are identical. logFC, F, PValue and FDR differ in every element (edgeR's iterative QL fit). Every element above the 1e-9 floor is within 5.1e-9 relative; the rest are within the floor. Largest relative figures overall: F 5.2e-6 at a value of 6.3e-8, logFC 9.65e-7 at 8.3e-5, both inside the floor. |
| 6 | `trx-norm:30` | `best_res` of the automatic RUV grid | **Agree** | percentage 10, n_controls 462, k 2, `ctl_idx` and the 462 control genes are identical; separation and composite scores differ by 1.2e-15 relative. |
| 9 | `trx-summary:52` | RUVIII `newY` (12 × 4,622) | **Agree** | 55,157 of 55,464 elements differ, maximum relative 2.1e-13 |
| 10 | `trx-summary:65` | sha256 of `transcriptomics_normalized_data.tsv` | **Agree** (follows #9) | The text that `vroom_write` produces from `newY`; its values are those of #9 |

No failure is class (c). Nothing is outside the rule. Nothing is discrete and different, except the text hashes that
record doubles within the rule.

### Against run 44041893594 (current repository versions, ubuntu-26.04)

| Then | Now |
|---|---|
| DIA, DIA_LIMPA, enrichment E2E (×6) timed out at `precursor_rollup` | **pass** |
| `trx-qc:74` and the control set already matched; manifest, `topTags`, `best_res`, `newY`, TSV differed | the same five digests differ, now measured: all within the rule |
| — | `release-gate:48-49` (the workflow's own; repaired) |

**Correction to the earlier diagnosis.** This record first attributed the DIA timeouts to a DIA manifest missing from
the checkout under `load_all`. That was wrong.

- `pkgload::load_all` runs `configure`: pkgbuild calls `R CMD INSTALL`. On the laptop it stopped in `build_dia_rsum.R`
  when duckplyr was not on the library path; on the runner it writes the manifest.
- The cause is in `newDIARelationalScope` (`R/utils_workflow_dia_relational_rollup.R:76-77`): the rollup stops unless
  `packageVersion("duckdb") == "1.5.5"`, and run 44041893594 had installed duckdb 2.0.0.
- That is class V (a version difference), not a harness defect. With the snapshot's duckdb 1.5.5 the same lanes pass.
- The relative E2E artifact path was a real defect. Its repair is why this run's E2E outputs were uploaded.

### Against the baseline (run 43916003934: R 4.6.2, no internal packages)

| Baseline | Now |
|---|---|
| `lipid-qc:342`, `metab-qc:293` (F-24) | gone since 44041893594 (`fastduckdb` present) |
| `trx-qc:74` and the `trx-norm:30` control set (changed by R 4.6.2) | gone since R was pinned at 4.6.1 |
| trx manifest, `topTags`, `best_res`, `newY`, TSV digests | remain as digests; within the rule |
| DIA and enrichment timeouts | gone (duckdb 1.5.5) |

## What stays open

- **The correction run was not used.** The workflow's own failure is repaired and verified on the laptop. A run of
  `79a03479` would show the release-gate file passing on CI.
- **Pull requests.** No run has exercised the empty-token path (fork and Dependabot pull requests).
- **m-04.** A cache conflicts with the internal-package rule. A concurrency group was outside this item.
- **The digest oracles stay bit-exact.** The laptop's exact comparison is still the gate for changes to the package.
  On any other machine, the five transcriptomics digests are expected to differ and are judged by
  `save-trx-products.R` plus `compare-results.R`.
- **Public logs** quote a few lines of `fastlimpadpc` source in rustc warnings. The user accepted this on 2026-11-17:
  the repositories go public before release. In this workflow the installer now writes the internal build output to
  `RUNNER_TEMP`.
