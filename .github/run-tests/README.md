# run-tests: MultiScholaR tests, check and proof scripts on GitHub-hosted runners

Branch `run-tests` is run-only: it carries `.github/workflows/run-tests.yml` and the request
`.github/run-tests/params.json`. A push that changes `params.json` starts one run. No test, build,
check or install runs on the laptop.

## Launch a request

1. The commit under test must be on origin: push its branch first (`git push origin <branch>`;
   never main). A local-only commit stops the run at `plan`.
2. In a worktree of `run-tests` (never `git stash`):
   ```
   git -C <repo> worktree add <dir> run-tests      # once
   cd <dir> && git pull --ff-only origin run-tests
   ```
3. Edit `.github/run-tests/params.json` (fields below), then
   ```
   git add .github/run-tests/params.json
   git commit -F <message file>
   git push origin run-tests
   ```
4. Wait and read: `.github/run-tests/fetch.sh` (latest run) or `.github/run-tests/fetch.sh <run id>`.
   It downloads the `results` artifact to `./run-tests-<run id>/` and prints the tables.

A push that does not change `params.json` starts nothing. `workflow_dispatch` reruns the branch's
current params. `gh run rerun <id> --failed` reruns failed shards only (the `results` job reruns too).

## params.json

| field | default | meaning |
|---|---|---|
| `label` | required | free text naming the request (ticket, purpose), one line |
| `package_commit` | required | 40-digit commit under test; its tree is checked out (path `pkg`) and its `tests/testthat` run |
| `tools_ref` | `""` | where the library tooling is read from (`.github/actions/snapshot-library`, `packaging/r/install-library.R`, `packaging/multischolar.lock.json`). Empty: `package_commit` when its tree has all three, else branch `kit-134`. Or a commit or branch name |
| `tests` | `[]` | `test-*.R` names under `tests/testthat`, or `"all"` |
| `shards` | 8 | matrix jobs (1 to 20; never more than the files); files dealt round-robin in the listed order |
| `file_timeout_minutes` | 45 | per file; a file killed at the limit is `timeout` and the shard carries on |
| `extra_packages` | `[]` | packages the library lacks, from the lockfile's snapshot (CRAN and Bioconductor 3.23 jammy binaries) with their hard dependencies, e.g. `["sva"]` (brings genefilter, annotate) |
| `env` | `{}` | environment variables for the test processes (not `GITHUB_*`, `RUNNER_*`, `ACTIONS_*`, `R_LIBS*`, `PATH`); `NOT_CRAN=true` is always set |
| `check` | `false` | a separate job: `R CMD build --no-manual`, then `R CMD check --no-manual --no-tests` with `_R_CHECK_FORCE_SUGGESTS_=false` (the tests are the test jobs' work) |
| `scripts` | `[]` | `{ "path": ..., "args": [...], "env": {...}, "out": <dir>, "timeout_minutes": 120 }`: `Rscript <path> <args>` from the package root after `package_commit` is installed (`R CMD INSTALL`, full git history present); `out` (relative to the package root) is uploaded |
| `os` | `ubuntu-22.04` | Linux only. Windows is refused: the kit-134 workflows install on Windows only inside the distributable bundle build (`packaging/windows/assemble.ps1`), not through a reusable library action |

Example (the KIT-129 RUV files with sva):
```json
{
  "label": "KIT-129 RUV files on kit-129-integration 33cfc99b",
  "package_commit": "33cfc99bdf38605f969a293ad2d6f991e4eb2e0d",
  "tests": ["test-ruv-factor-count.R", "test-trx-ruv-designs.R"],
  "shards": 2,
  "extra_packages": ["sva"],
  "check": false
}
```

## What a run does

Per job (tests shard, check, scripts): R 4.6.1 from Posit r-builds and the KIT-134 snapshot library
(`./.github/actions/snapshot-library` of `tools_ref`), with the two internal packages `fast_duckdb`
and `fast_limpa_dpc` built in that job from the secret `DEPS_READ_TOKEN`, plus `pkgload` and
`testthat`; then `extra_packages`. The tools' own MultiScholaR is removed in the test jobs; each test
file runs in its own process: `Rscript --vanilla .github/run-tests/run-file.R` =
`pkgload::load_all(pkg, export_all = TRUE, helpers = TRUE)` and
`testthat::test_file(<file>, stop_on_failure = FALSE)` with the summary reporter (printed to the
file's log) and a list reporter (the counts), as the laptop's per-file runs
(`tools/testing/test-baseline.R`, load mode `load_all`). Native code is compiled once before the files.

## Results: artifact `results`

| file | content |
|---|---|
| `tests.tsv` | file, passed, failed, errors, skipped, warnings, seconds, status (`ok`, `failed`: an expectation failed or errored, `timeout`, `crashed`: no such file, R died, or its shard job was lost) |
| `failures.txt` | per failing or erroring expectation: file, test, kind, message, the first trace lines; crashes with the log's last lines |
| `logs/` | each file's console |
| `check-summary.txt`, `check/` | the check's counts and text of ERRORs, WARNINGs, NOTEs; MultiScholaR's `00check.log`, `00install.out`, build log |
| `scripts/` | `scripts.tsv` (index, path, exit, seconds, status), per script `console.log` and `out/` |
| `meta.txt`, `library.tsv` | run, request commit, label, package_commit, tools_ref, R version, job results, setup seconds per shard, env, library table |
| `summary.md` | the job summary page's text |
| `withheld.txt` | files not uploaded because they matched a forbidden pattern (only if any) |

The run is red only on infrastructure failure (a job failed or was cancelled), a file that crashed or
timed out, or a failed script. Failing expectations and check findings are results: green run, see
`tests.tsv`.

## Token rules (absolute)

- `DEPS_READ_TOKEN` is used only by the library action's install step. Never print it, never write it
  to a file, never pass it to a test or script.
- No cache of any kind (no `actions/cache`, no setup-action cache).
- Nothing built with the token, and no install or build output of the internal packages, is uploaded.
  The library is rebuilt in every job; the action keeps its logs under `RUNNER_TEMP`.
- Every upload is scanned first (`.github/run-tests/scan.R`, the patterns of the KIT-124 workflow:
  compiler and Rust build output, the internal repositories' names, token shapes; compiled objects and
  archives): a matching file is withheld and named in `withheld.txt`. The assembled `results` is
  scanned again and a match fails the run before upload. The check log may hold MultiScholaR's own
  compiler lines.
- Limits: 300 minutes per test or scripts job (the shard stops starting files 12 minutes before),
  public repository's standard runners (4 cores, 16 GB). The library install is overhead in every job.
