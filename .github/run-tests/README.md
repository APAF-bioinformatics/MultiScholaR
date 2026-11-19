# run-tests: MultiScholaR tests, check and proof scripts on GitHub-hosted runners

Branch `run-tests` is run-only: it carries `.github/workflows/run-tests.yml`, the dispatcher
`.github/workflows/run-tests-requests.yml` and the requests. **Each request is its own file**
`.github/run-tests/requests/<name>.json`. A push that adds request files starts one independent run per
added file, titled `run-tests <name>`. No test, build, check or install runs on the laptop.

## The rule: one file per request

- Start every request from a blank template: `.github/run-tests/new-request.sh <name>`. Never copy,
  edit or reuse another request's file; nothing is inherited from anyone else's request.
- A request file is used once. The dispatcher refuses a changed request file, and a name that the
  branch's history already holds. To run again or correct a request, add a new file with a new name
  (e.g. `-2` suffix), stating every field again.
- `<name>`: a short slug, a-z, 0-9, `.`, `_`, `-`, up to 80 characters, e.g. `kit129-figures-303217ac`.

## Launch a request

1. The commits under test must be on origin: push their branch first (`git push origin <branch>`;
   never main). A local-only commit stops the run at `plan`.
2. In a worktree of `run-tests` (never `git stash`; never force-push: other agents push here):
   ```
   git -C <repo> worktree add <dir> run-tests      # once
   cd <dir> && git pull --ff-only origin run-tests
   ```
3. Write the blank request and fill it in (fields below; `requested_by`, `label`, `package_commit`,
   `base_commit` and `tests` are required):
   ```
   .github/run-tests/new-request.sh <name>
   $EDITOR .github/run-tests/requests/<name>.json
   ```
4. Commit only that file and push (if the push is rejected, `git pull --ff-only origin run-tests` again):
   ```
   git add .github/run-tests/requests/<name>.json
   git commit -F <message file>
   git push origin run-tests
   ```
5. Wait and read: `.github/run-tests/fetch.sh <name>` (the latest run titled `run-tests <name>`; it
   exists about a minute after the push) or `.github/run-tests/fetch.sh <run id>`. It downloads the
   `results` artifact to `./run-tests-<run id>/` and prints the tables. The run's page starts with the
   resolved request (requested_by, head, base, tools, test files, scripts): check it is yours.

The push starts `run-tests requests` (job `dispatch`), which dispatches `run-tests.yml` once per added
file (`workflow_dispatch`, input `request = <name>`); its summary lists the dispatched requests. Runs
are independent and never cancel each other. Several request files in one push are several runs.
`gh workflow run run-tests.yml --ref run-tests -f request=<name>` reruns a request by hand;
`gh run rerun <id> --failed` reruns failed shards only (the `results` job reruns too). Never cancel
a run that is not yours.

The `plan` job checks the request and fails at once, with the reason on the run page, when: a required
field is missing; `requested_by` is empty; `package_commit` is not a full 40-digit hash; `base_commit` is
neither a full 40-digit hash nor `"none"`; `label` is not 1 to 200 characters on one line; any field is
outside the rules below; or a commit is not reachable on origin.

## Deprecated: params.json

For a transition period a push that changes `.github/run-tests/params.json` still starts a run on that
file (and `workflow_dispatch` without a request reads it), with a deprecation warning on the run page.
This is the path where requests inherited each other's base, tests and scripts: do not use it for new
requests. It will be removed once no params.json run has been started for a week.

## Request fields

| field | default | meaning |
|---|---|---|
| `requested_by` | required (request files) | who asks and why: ticket and purpose, one line (at most 300 characters); printed at the top of the run page and in `meta.txt` |
| `label` | required | free text naming the request (ticket, purpose), one line, 1 to 200 characters |
| `package_commit` | required | full 40-digit commit under test (head); its tree is checked out (path `pkg`) and its `tests/testthat` run |
| `base_commit` | required (request files) | full 40-digit commit to compare against, or `"none"` for no base run: the same files (same shards, tools, extra packages, runner image) run on it in separate jobs, and `check` runs on both; see "Base versus head" |
| `tools_ref` | `""` | where the library tooling is read from (`.github/actions/snapshot-library`, `packaging/r/install-library.R`, `packaging/multischolar.lock.json`). Empty: `package_commit` when its tree has all three, else branch `kit-134`. Or a commit or branch name |
| `tests` | required (request files; may be `[]`) | `test-*.R` names under `tests/testthat`, or `"all"` |
| `shards` | 8 | matrix jobs (1 to 20; never more than the files); files dealt round-robin in the listed order |
| `file_timeout_minutes` | 45 | per file; a file killed at the limit is `timeout` and the shard carries on |
| `file_timeouts` | `{}` | per-file timeout overrides, `{"test-gui-entry-proteomics.R": 120}` |
| `extra_packages` | `[]` | packages the library lacks, from the lockfile's snapshot (CRAN and Bioconductor 3.23 jammy binaries) with their hard dependencies, e.g. `["sva"]` (brings genefilter, annotate) |
| `env` | `{}` | environment variables for the test processes (not `GITHUB_*`, `RUNNER_*`, `ACTIONS_*`, `R_LIBS*`, `PATH`); `NOT_CRAN=true` is always set |
| `check` | `false` | a separate job: `R CMD build --no-manual`, then `R CMD check --no-manual --no-tests` with `_R_CHECK_FORCE_SUGGESTS_=false` (the tests are the test jobs' work) |
| `scripts` | `[]` | `{ "path": ..., "args": [...], "env": {...}, "out": <dir>, "out_pattern": ..., "timeout_minutes": 120 }`: `Rscript <path> <args>` from the package root after `package_commit` is installed (`R CMD INSTALL`, full git history present); `out` (relative to the package root) is uploaded; `out_pattern` (extended regex on the path under `out`) keeps only matching files, e.g. `"\\.(tsv|txt)$"` |
| `fetch_refs` | `[]` | branches the scripts read with `git show` (fetched as `origin/<name>` and as local branches), e.g. `["sim-k124-exact"]` |
| `blas` | `reference` | `reference`: the distribution's reference libblas3/liblapack3 in every job (as the laptop; the same results on every runner CPU model). `openblas`: the runner's default. Under OpenBLAS two jobs on different CPU models differ in the last digits (run 44356604347: KIT-073 correlations from the 12th digit), so base-versus-head comparisons need `reference` |
| `os` | `ubuntu-22.04` | Linux only. Windows is refused: the kit-134 workflows install on Windows only inside the distributable bundle build (`packaging/windows/assemble.ps1`), not through a reusable library action |

Example, `.github/run-tests/requests/kit129-ruv-files-33cfc99b.json` (the KIT-129 RUV files with sva;
fields left at their default omitted here, the template states them all):
```json
{
  "requested_by": "KIT-129: RUV design files after the integration merge, base versus head",
  "label": "KIT-129 RUV files on kit-129-integration 33cfc99b",
  "package_commit": "33cfc99bdf38605f969a293ad2d6f991e4eb2e0d",
  "base_commit": "none",
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
| `meta.txt`, `library.tsv` | run, request commit, request file, requested_by, label, package_commit, tools_ref, R version, job results, setup seconds per shard, env, library table |
| `summary.md` | the job summary page's text |
| `withheld.txt` | files not uploaded because they matched a forbidden pattern (only if any) |

The run is red only on infrastructure failure (a job failed or was cancelled), a file that crashed or
timed out, or a failed script. Failing expectations and check findings are results: green run, see
`tests.tsv`.

## Base versus head

Stored hashes recorded on the laptop cannot be judged on a runner by pass or fail. With `base_commit`
the rule is: a change is clean when no test that passes on base fails on head, and every expectation
that fails on both fails with the same message (temporary and checkout paths masked). `compare.tsv`
(file, base and head passed/failed/errors/skipped/status, verdict, tests, note) gives per file one of
`same`, `new-failure`, `message-changed`, `fixed`, `count-changed`, `only-on-head` (the base tree
lacks the file; base status `absent`). `compare-failures.txt` shows base and head messages side by
side for every non-identical failure. A failure identical on both whose message names a digest or
hash but shows no hexadecimal value is noted: its actual value cannot be compared. `check-compare.txt`
lists check findings only on head and only on base. Head results stay in `tests.tsv`/`failures.txt`;
base results are `tests-base.tsv`/`failures-base.txt`/`logs-base/`.

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
