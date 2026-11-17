# KIT-134 B-01: the measured Module CI run on a Linux runner (2026-11-17)

Audit finding B-01 (plan 4.7 audit) has three parts: the workflow had never passed its first job; every stored reference
is bit-exact from the laptop (`-march=native -O3`); CI resolved without `fastduckdb` and `fastlimpadpc`. Its closure
criterion is one full run on a Linux runner, read, and the comparison rule written. The rule is KIT-134 revision 2.0.4
(fixed before any run was read; not changed here): doubles agree when `|a - b| <= 1e-6 * max(|a|, |b|)` or
`|a - b| <= 1e-9`; discrete outputs, outputs whose contract is bit-exactness, and file hashes of copied inputs exactly.

## Workflow change (branch `ci-token`)

- `79bc4410`: R pinned at `4.6.1` on all 7 `setup-r` steps (IR-M3); `cache: false` and
  `GITHUB_PAT: ${{ secrets.DEPS_READ_TOKEN }}` (step env only) on all 7 `setup-r-dependencies` steps; the `=?ignore`
  lines of `5534aa97` removed. The internal packages come from DESCRIPTION's named Remotes (`3ec13e95`).
- `d2f02aa1` (audit F2): with no secret (fork and Dependabot pull requests) the two packages are ignored again.
  Not yet exercised by a run.
- No internal package in a cache or artifact: the action's cache steps run only for `cache: 'true'`/`'always'`;
  `setup-r` has no cache; all 11 upload paths are under `tests/testthat/_module_ci_artifacts` or `_e2e_artifacts`; the
  library is in `RUNNER_TEMP`. Checked after the run: its 3 artifacts are 6 files (impact routing, release-gate summary,
  the unit-contract manifest and scorecards), none naming either package.

## Runs

| Run | Ref | Outcome | Class |
|---|---|---|---|
| [44013986162](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44013986162) | `1e4eb7a7` | every install: "Can't find package called fastduckdb, fastlimpadpc" (pak names an unnamed GitHub Remote after the repository) | (a); fixed by named Remotes |
| [44014345089](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44014345089) | `6767c8ce` | every install: `fastlimpadpc` d4f9cc1d "lazy loading failed" (no `Imports: limpa`) | (a); fixed upstream, 87a90b89 |
| [44040978386](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44040978386) | `8a89123d` | every install: `fastduckdb` 2694b972 `cargo --offline` without vendored crates | (a); fixed upstream, e724f90b |
| **[44041893594](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44041893594)** | `79bc4410` | **the measured run**, below | |

Push-triggered duplicates (44013980819, 44014335748, 44040970699, 44041893072) were cancelled; each was a subset of the
dispatched `gate = all` run.

## The measured run 44041893594 (gate `all`, ubuntu-26.04 image 20261111.196, R 4.6.1)

| Job | Result | Dependency install | Wall |
|---|---|---|---|
| detect impact | pass (impact `none`: a dispatch has no before-SHA, so the push-gate lanes were skipped) | 38.7 min | 39 min |
| full gate (cross-omic-coexistence) | pass | 29.3 min | 30 min |
| full gate (alternate-launch-browser-smoke) | pass | 32.2 min | 33 min |
| full gate (report-export-fidelity) | **fail**: 2 | 33.1 min | 41 min |
| full gate (all-module-unit-contracts) | **fail**: 8 | 39.2 min | 54 min |
| full gate (all-current-e2e-workflows) | **fail**: 5 | 37.4 min | 67 min |
| release candidate promotion gate | fail (follows the full gates) | | |

Installed from source without a cache: `fastduckdb` 0.2.0 (e724f90, 3 min 28 s), `fastlimpadpc` 0.1.0 (87a90b8, 23 s).
The runner installs current repository versions, not KIT-134's lockfile (the verified library): edgeR 4.10.7 (lock
4.10.1), limma 3.68.5 (3.68.4), ruv 0.9.7.2 (0.9.7.1), duckdb 2.0.0 (1.5.5), duckplyr 1.2.2 (1.2.1); limpa 1.4.2,
RUVSeq 1.46.0, EDASeq 2.46.0 equal. A difference below is therefore machine **or** version (class V of 2.0.4) until a
run installs the lockfile.

Everything else the gates compare matched the laptop exactly: all lipidomics, metabolomics and proteomics module-CI
files (including `module-ci-oracles`, 29 expectations), the transcriptomics import, design, app and QC filter oracles,
the lipidomics, metabolomics, TMT/LFQ and transcriptomics-screen E2E lanes, cross-omic and browser smoke. Skips are the
declared quarantines (Q-22, Q-23, Q-76, Q-173, Q-174, Q-182) and the audit-candidate guards.

### Every failure, classified

| # | Job | Test | What differs | Class | Under the rule |
|---|---|---|---|---|---|
| 1-4 | unit-contracts | `trx-da:7`, `trx-norm:13`, `trx-qc:64`, `trx-summary:40` (`trxDeclReference`) | sha256 of the regenerated corpus `manifest.json`: runner `faa19db7…` (identical under R 4.6.2 in the baseline), laptop `e17d1ed5…` | (b), provisional | Not applicable from this run. The manifest holds per-file hashes of `counts.tsv` (integers) and `truth.json` (doubles from `rnorm`/`runif`/`exp`/`rlnorm`, written with `digits = NA`). The counts are very probably identical: the filtered count matrix (`trx-qc:74`) and the inferred control set (`trx-norm:30` `best_control_genes`) match the laptop bit for bit. The likely difference is last-digit doubles in `truth.json`, which no product reads. Settling it needs the runner's per-file hashes, which the run does not record. |
| 5 | e2e | `trx-journey:81` | the same manifest hash | (b), provisional | as 1-4; every other expectation in the journey passed |
| 6 | unit-contracts | `trx-da:28` | sha256 of the first contrast's `topTags` table (edgeR QL: logFC, logCPM, F, PValue, FDR) | (b) or V | No maximum relative difference can be given: the KIT-040 store holds digests of the products, not values, and the runner's objects are not saved. edgeR 4.10.7 against 4.10.1. |
| 7 | unit-contracts | `trx-norm:30` | sha256 of `best_res` (percentage, n_controls, best_k, separation_score, k_penalty, composite_score, ctl_idx, control_genes) | (b) or V; the discrete fields not separable | The control set alone matches exactly (`best_control_genes`), so k's selection inputs agree; whether k, percentage and n_controls are identical is inside the digest. Needs both sides' values. |
| 8 | unit-contracts | `trx-summary:52` | sha256 of RUVIII `newY` | (b) or V | Needs values. |
| 9 | unit-contracts | `trx-summary:65` | sha256 of `transcriptomics_normalized_data.tsv` | follows #8 | The file is `newY` written by `vroom_write`; a text hash of doubles. |
| 10-13 | e2e | `proteomics-dia:21`, `proteomics-dia-limpa:21`, `proteomics-enrichment:151`, `:155` | "Timed out waiting for proteomics R6 state precursor_rollup; last current was qvalue_filtered" after the full 180 s wait | **(a) harness** | No result produced; see below. |
| 14-15 | report-export | `proteomics-enrichment:151`, `:155` | the same timeout | **(a) harness** | as 10-13 |

No failure is class (c) on the evidence of this run.

**Why the DIA rollup never commits (10-15).** The GUI's in-process rollup (`runDIAReferencePrecursorRollup`, P14)
first checks the native runtime in `diaRollupNativeInfo()`, which reads
`system.file("duckdb", "linux_amd64", "manifest.rds", package = "MultiScholaR")`. That file is generated at installation
by `configure` (`src/build_dia_rsum.R`), is gitignored (`/inst/duckdb/`), and records hashes of the machine's own
`duckdb`, `duckplyr` and `fastduckdb` binaries. The gates run the app under `pkgload::load_all()` from a clean checkout:
the job log shows the package compiled without `configure` (no `build_dia_rsum` output), so the file does not exist, the
rollup aborts inside the observer's `tryCatch` (an error notification), and the state stays `qvalue_filtered` until the
harness times out. The laptop checkout has the file from an installation on 2026-11-15. The KIT-134 token build ran the
same DIA lanes in process from an installed package and matched the laptop at every step (run 44044598918), so the
rollup itself works on a runner. This is an inference: the app's own log is not available, because the failure
artifacts are written under a relative `MULTISCHOLAR_E2E_ARTIFACT_DIR` that `testthat::test_dir` resolves inside
`tests/testthat/` (to `tests/testthat/tests/testthat/_e2e_artifacts/…`), which no upload path covers. That is a second
harness defect; the e2e and report-export jobs uploaded nothing.

### Against the baseline (run 43916003934: performance-v2 at `5534aa97`, R 4.6.2 unpinned, no internal packages, impacted lanes)

| Baseline failure | Now |
|---|---|
| `module-ci-lipid-qc:342`, `module-ci-metab-qc:293` (F-24 `after_click`) | **gone**: both files pass in full (the `fastduckdb` reader is present) |
| `trx-qc:74` filtered count matrix; `trx-norm:30` `best_control_genes` (a discrete control set) | **gone**: equal to the laptop under R 4.6.1 (R 4.6.2 changed them, as recorded when R was pinned) |
| trx manifest ×4, `topTags`, `best_res`, `newY`, the normalized TSV | **remain** (#1-4, 6-9); manifest hash identical in both runs |
| DIA E2E and enrichment timeouts at `precursor_rollup` | **remain** with the packages installed (#10-15): the cause is the missing generated manifest, not the absent package |
| — | **new in scope**: `proteomics-dia-limpa` and `trx-journey` (the full E2E gate runs them; the baseline's impacted lanes did not) |

## Verdict on B-01

**Not closed.** Parts (1) and (3) are answered: every job installs both internal packages from their pinned commits
and reaches its tests, and three of the six test jobs pass. Part (2), the measurement, is incomplete:

- The DIA, DIA_LIMPA and enrichment lanes produced no result (harness, #10-15), so their cross-machine difference is
  not measured.
- The transcriptomics differences (#1-9) are counted (5 assertions in 4 module files and 1 E2E file; 4 products and
  1 corpus file) but cannot be classed under the rule. Every oracle involved is a digest. The runner's objects are not
  saved, and the runner's package versions are not the lockfile's, so machine and version differences cannot be told
  apart.

No result outside the rule has been shown. No result has been shown to be inside it either.

## Next run (needs authority; this budget is spent)

One Module CI run, `gate = all`, on `ci-token`, to finish the B-01 measurement. It needs three harness changes and
one decision:

1. In the jobs that `load_all`, run `"${R_HOME}/bin/Rscript" --vanilla src/build_dia_rsum.R` in the checkout after the
   dependency install (what `configure` does at installation). Alternatively, `R CMD INSTALL` the checkout and test the
   installed package.
2. Make `MULTISCHOLAR_E2E_ARTIFACT_DIR` absolute (`${{ github.workspace }}/tests/testthat/_e2e_artifacts/…`) so
   failure artifacts upload. They hold screenshots and app logs, not a library.
3. Save the runner's values of the four transcriptomics products and the corpus's per-file hashes into the artifact
   directory: `topTags` for contrast 1, `best_res`, `newY`, and the `manifest.json` files. These are package outputs,
   not internal code. Then apply `compare-results.R` on the laptop against the products in
   `MultiScholaR-dia-proof-20260911-VhC1n5/kit/trx-reference/`. This is a `tools/` change and outside this item.
4. Decision for the user: install the KIT-134 lockfile's versions on the runner (P1) so that V can be separated from
   machine, or accept that this run measures both together.

Exit for that run: the DIA lanes either pass or fail with a reason; every remaining difference has a 2.0.4 class with
its maximum relative difference. Budget: 1 run plus 1 correction for harness defects.

## Not addressed

- m-04 (cache and a concurrency group). A cache conflicts with the internal-package rule; a concurrency group was outside this item.
- The 30 to 40 minute uncached install now runs in every job.
- Public logs quote a few lines of `fastlimpadpc` source in rustc warnings. The user accepted this on 2026-11-17: the repositories go public before release.
