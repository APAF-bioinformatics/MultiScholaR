# KIT-245 records

## G0 declaration (written before the dispatch)

- **Purpose.** Observe, in the installed bundles, what "Push to GitHub" (`pushProjectToGithub`) does today, before any
  package code changes (ticket r1.2, sections 8 to 10, A0): T0 (no token), P1 (no `.Rproj`), P2 (`.Rproj`; the push
  fails), P2n (no network, no redirect: git's own failure text), P2g (git present, redirected to a local bare
  repository: the pushed tree), P3 (two templates), P4 (`git2r::push` with credentials to a local bare repository),
  P5 (the Summary panel of the four omics). Hosts: the installed Windows application (shell's environment); the
  installed Linux bundle on the runner, in `ubuntu:22.04` (no git) and in `buildpack-deps:jammy-scm` (git), both
  `--network none`.
- **Binding commit.** The product tree of `v0.9-staging` at send time, **16112917** (16112917329a63ab062655f04db5bf1dd7df9729),
  checked by the probe step (`val_product_tree_ok`). `R/func_github_management.R`, the GitHub panel of
  `R/mod_omic_summary.R`, `packaging/` and the workflows are unchanged from 3eeaf100 (`git diff 3eeaf100 16112917`).
  Run branch `run-k245-g0` (never merged): the probe file and the dispatch-only stage `probe` of
  `.github/workflows/distributable.yml`, after KIT-241's pattern (KIT-241's own `run-k241-g0` was dispatched as run
  45577719573 and is not touched).
- **Safety.** No real token is used, created, read or stored: the four token variables are empty and asserted so; the
  token is the plainly fake `kit245-dummy-not-a-token`; `gh::gh` is replaced by a recorder that never calls the
  original and stops on any other token; no request reaches github.com: on the runner git is redirected to a local path
  (`insteadOf`), limited to `GIT_ALLOW_PROTOCOL=file` and given a dead proxy, the containers have no network; a local
  bare repository is the only remote. Nothing is uploaded or cached (the stage has no upload or cache step;
  `DEPS_READ_TOKEN` reaches only the install steps). Every probe output is printed with the dummy masked.
- **Ceiling.** 125 job-minutes (8.33 core-hours; ticket section 11), expected about 61 (KIT-241's probe stage took 47).
- **What the result decides.** Stop rule 1: the installed Windows R finds `git`, or P2's push step does not fail: stop,
  no code. Stop rule 7: git on the runner reaching for the network, a token other than the dummy at the recorder, or a
  token-shaped value in the environment: instrument defect, stop. Otherwise A0's observations become the base values of
  Part S's rows and the facts given to the owner for OD-1; stop rule 2's parts that G0 observes (git2r present, its
  HTTPS feature, P4) are recorded for option (a) only. A failure of the probe's mechanics with no result is an
  instrument defect (VAL-54: at most three reruns, recorded here).

## G0 result (runs 45587834788 and 45593577523; product tree 16112917, `val_product_tree_ok` in both)

Windows from run 45587834788; Linux from run 45593577523 (the VAL-54 instrument rerun: in the first run step 0 stopped
on the runner and in jammy-scm because the probe passed a git regex unquoted through `system2`; the fix quotes every
git argument and makes a failing git query stop step 0; no expected value, comparison or stop rule changed). Every
value: `g0-base-values.tsv`. Step 0 passed on every host; every invocation printed its END line.

### Observed

| Case | Windows installed | Linux runner (git 2.56.1) | ubuntu:22.04 (no git) | jammy-scm (git 2.34.1) |
|---|---|---|---|---|
| environment | no git on R's PATH; `interactive()` FALSE; git2r 0.36.2, gh 1.6.1, gitcreds 0.1.2 in the bundle; libgit2 https, ssh, threads TRUE | `/usr/bin/git` | none | `/usr/bin/git` |
| T0, no token | `gh_token()` empty (class `gh_pat`); push: "DEBUG: Failed to retrieve GitHub token from gh::gh_token()"; nothing staged, recorder 0 | same | same | same |
| P1, no `.Rproj` | "R project file does not exist: ..."; recorder 0; **staged copy left in `tempdir()`** | same | same | same |
| P2, `.Rproj` | **`POST /orgs/{org}/repos` called (private TRUE, auto_init FALSE)**, then "Failed to push to repository. Error during push step: 'git' not found"; staged folder left; **its `.git/config` holds the token** | POST, then exit 128 at the redirected local path ("does not appear to be a git repository"); token left in `.git/config` | POST, then "error in running command"; token left | as the runner |
| P2n (jammy-scm, no network, no redirect) | | | | POST, then exit 128; git: "fatal: unable to access 'https://github.com/kit245-probe-no-such-org/k245-p2n.git/': Could not resolve host: github.com"; **no token in git's text** (captured re-run: no; unmasked dummies in the whole invocation output: 0); token left in `.git/config` |
| P2g (git, redirect to a local bare repository) | | **the base push succeeds**; staged folder removed | | same |
| P2g pushed tree | | `.gitignore`, `k245-p2g.Rproj`, `scripts/proteomics/`: `k245-p2g.Rmd`, `design_matrix.tab`, `study_parameters.txt`, `import_sources/x.tsv`, `data_cln_LCMS_Pos.tab`, `filtered_session_x.msdia.tar`, `session_state_x.RDS`; no file of the remote holds the token | | same |
| P3, two templates | "menu() cannot be used non-interactively"; recorder 0; staged copy left | same | same | same |
| P4, `git2r::push` with credentials to a local bare repository | pushed (branch `master`) | same | same | same |
| P5, Summary panel of the four omics | `push_to_github`, `enable_github`, `github_org` and the heading present for proteomics, lipidomics, metabolomics, transcriptomics | same | same | same |

Token in printed output: 0 unmasked dummies in every invocation's output on every host; the one unmasked occurrence in
each job log is the step's own source line that defines the fake string. Stop rule 7: no network phrase in the
runner's output; the only such lines are P2n's (no network by construction, exempt).

### Inferred

- Stop rule 1 does not fire: the installed Windows R finds no git and P2's push step fails. The work may proceed.
- A0 is met as written in r1.2, with one difference: with `GITHUB_PAT` set to the fake string, gh 1.6.1's `gh_token()`
  did not return it (no error; it returned something else, presumably because the string is not shaped like a token),
  so P1 to P3 used the stand-in for `gh_token` that the existing tests use; T0 used the real `gh_token()`.
- Stop rule 2's parts G0 can see hold for option (a) only: git2r is in both bundles, libgit2 reports HTTPS, P4
  pushes with credentials (local transport; the HTTPS transport and its certificate store are not observed).

### Not checked

- An authenticated HTTPS push or any GitHub answer (no token, no network by rule).
- Whether a credential helper would store the URL's credentials after a successful system-git push (no helper existed).
- Whether a forced termination leaves R's temporary folder on Windows.
- Linux with the systemd launcher (the probe runs the installed Rscript directly).
- What gh 1.6.1's `gh_token()` returned for the fake string (only "not the dummy" was recorded).

## Part S (branch `kit-245` from `v0.9-staging` e0747fc1, after the merge train carrying KIT-234)

At the branch point `R/func_github_management.R`, `tests/testthat/test-github-management-direct-shared.R` and the
two `man/` pages are as at 3eeaf100 (every cited line holds); KIT-234's W8 is merged
(`tests/testthat/test-kit234-workbooks-withdrawn.R:194-230`, the md5 assertion at `:228-229`).

### R1s (declared before the request)

- Purpose: show S1 to S5 failing by assertion at their declared base values on the tests-only commit; G1, G2 and
  `test-github-management-direct-shared.R` (with its declared `unlink` edit) passing.
- Binds 009fdb38 (tests only; nothing under `R/`). Ceiling 45 job-minutes. Decides whether the fix is written
  (stop rule 3: an S row passing is dropped; S1 passing stops the work).
- Result, run 45598693322 (`r1s-tests.tsv`, `r1s-failures.txt`): as declared. S1: not classed, no message naming the
  repository, the staged folder left, the clone's `.git/config` holding the token at the push, no environment-borne
  configuration, `k245-s1/.git/config` under `tempdir()` holding the token. S2: POST made, staged folder left. S3 and
  S5: the call returned TRUE after the POST. S4: the sentinel's file beside `tempdir()` deleted. G1, G2 pass; the
  existing file 21 passed, 0 failed. 11.1 job-minutes.

### R2s (declared before the request)

- Purpose: base (the tests-only commit 009fdb38) against head (the fix): the new file passes on head; every other
  regression file `same` except the declared edit of KIT-234's W8 (b) (`fixed` or `same` on head with the edit;
  base runs it unedited and passes); `R CMD check` on both, nothing only on head (the two `man/` pages changed,
  declared). Files: the grep of section 10 at the branch point (12 files: the ten of r1.2, KIT-234's file, the new
  file).
- Binds the fix commit against 009fdb38. Ceiling 300 job-minutes (ticket section 11). Decides AS2; stop rule 4 (any
  other file not `same`, or a new check finding on head) stops the work.
