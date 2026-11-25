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
- Result, run 45600718143 (`r2s-compare.tsv`, `r2s-check-compare.txt`): AS2 met. `test-kit245-github-push-safety.R`
  `fixed` (base 13 passed, 18 failed; head 31 passed, 0 failed). The other eleven files `same`: the edited
  `test-github-management-direct-shared.R` 21/0 on both; KIT-234's file 29/0 on both (W8 unedited on the base,
  edited on the head); `test-module-ci-trx-summary.R` (3), `test-omic-summary-module.R` (28) and
  `test-optional-packages.R` (1) fail identically on both sides (stored references recorded on the laptop, brief
  rule 4); the lipid, metabolomics and proteomics summary files run 0 expectations on both. `R CMD check`: base and
  head 0 errors, 3 warnings, 3 notes, the six findings identical; none only on head. 59.9 job-minutes.

## Product part, option (b) (the owner's decision of 2026-11-25 evening, VAL-244; OD-1b: keep sending, say so)

Step 0: ticket r1.3 records the decision and the (b) budget line. Detection as the ticket names it (E-1, 3.4):
`nzchar(Sys.getenv("MSR_BUNDLE"))`.

### R1p (declared in the request before it was sent)

- Purpose: E1 fails by assertion on the tests-only commit (the push inputs present, no sentence); G1 passes. Binds
  c85267e8. Ceiling 45 job-minutes. Decides whether the (b) fix is written.
- Result, run 45604664357 (`r1p-tests.tsv`): as declared: 12 passed, 28 failed, 0 errors (E1: six inputs present
  and the sentence absent for each of the four omics; the heading present; G1 passes). 7.1 job-minutes.

### R2p (declared before the request)

- Purpose: base (c85267e8) against head (the (b) fix): `test-kit245-github-panel.R` passes on head; every other file
  of section 10's regression list `same` (the frozen UI digests of `test-omic-summary-module.R`, the UI-rendering
  files: tests run without `MSR_BUNDLE`); `R CMD check` nothing only on head. Binds the fix against c85267e8.
  Ceiling 460 job-minutes (section 11). Decides A2; stop rule 4.
- Result, run 45606131371 (`r2p-compare.tsv`, `r2p-check-compare.txt`): A2 met. `test-kit245-github-panel.R` `fixed`
  (base 12 passed, 28 failed; head 40 passed, 0 failed). The twelve other files `same`, among them
  `test-omic-summary-module.R` (1,675 passed and the same 28 stored-reference failures on both sides) with the frozen
  UI digests, `test-omic-gui-config.R`, `test-module-ui-wrapper-direct-shared.R` and `test-module-ci-trx-summary.R`;
  Part S's file 31/0 on both. `R CMD check`: the six findings identical on both sides, none only on head.
  86.9 job-minutes.

### W1 (declared before the dispatch)

- Purpose: A3 for (b) (P5: no push input and the sentence present for the four omics in the installed bundles) and
  AS3 (P2 stops at the git check with no recorder call where git is absent; after P2 and P2g no file under
  `tempdir()` holds the token). Binds the product tree 9efafdad (`val_product_tree_ok`); run branch `run-k245-w1`
  (never merged): G0's stage, the redirect extended to Part S's token-free URL (`insteadOf = https://github.com/`),
  the probe's P5 reporting the sentence. Ceiling 125 job-minutes. Decides A3, AS3 (stop rule 5).
- Result, run 45606206919 (`w1-head-values.tsv`; the probe's lines keep the label "G0"): met on all five hosts.
  Step 0 passed everywhere; every END line printed; 0 unmasked tokens in any invocation's output.

| Case | Windows installed | Linux runner | ubuntu:22.04 | jammy-scm |
|---|---|---|---|---|
| P1 (no `.Rproj`) | "R project file does not exist"; staged folder removed | same | same | same |
| P2 | `multischolar_github_git_missing`, recorder 0, nothing staged | POST, then `multischolar_github_push_failed` naming `kit245-probe-no-such-org/k245-p2`; staged folder removed | as Windows | as the runner |
| P2n (no network) | | | | POST, then the message naming the repository; git: "Could not resolve host: github.com", no token in its text; nothing left |
| P2g | | push succeeds; tree as at G0; nothing holds the token | | same |
| P3 | the `menu` error; staged folder removed | same | same | same |
| P5 (four omics) | push_to_github, enable_github, github_org absent; heading and the sentence present | same | same | same |
| files under `tempdir()` holding the token, after every case | none | none | none | none |

## Acceptance criteria

| # | Criterion | State | Evidence |
|---|---|---|---|
| A0 | G0 on base | met | runs 45587834788, 45593577523; `g0-base-values.tsv` |
| AS1 | Part S tests-only: S1 to S5 fail at their base values; G1 to G3 pass | met | run 45598693322; `r1s-tests.tsv`, `r1s-failures.txt` |
| AS2 | Part S fix against tests-only: new file passes, 10 regression files `same` (declared edits), check nothing only on head | met | run 45600718143; `r2s-compare.tsv`, `r2s-check-compare.txt` |
| AS3 | installed bundles: P2 stops at the git check with no recorder call where git is absent; no file under `tempdir()` holds the token after P2 and P2g | met | run 45606206919; `w1-head-values.tsv` |
| A1 | product tests-only: E1 fails, G1 passes | met | run 45604664357; `r1p-tests.tsv` |
| A2 | product fix: new file passes, every other file `same`, check nothing only on head | met | run 45606131371; `r2p-compare.tsv`, `r2p-check-compare.txt` |
| A3 | (b) P5: no `push_to_github` and the sentence present for the four omics, installed | met | run 45606206919 |
| A4 | full suite at the next checkpoint | open (rides the merge train) | |
| A5 | NEWS, ledger rows, this file with "For the change report"; `NAMESPACE`, `DESCRIPTION`, `packaging/` unchanged | met | `git diff e0747fc1 -- NAMESPACE DESCRIPTION packaging` empty |

## Observed and inferred

Observed: every value in the tables above and in the result files. Inferred, not run: that Windows R users with git
2.31 or later can push (no Windows run of the R function with git exists; see "Not checked"); that "private" means
visible to the organisation's owners and to the members its settings allow (GitHub's model, not re-read).

## Not checked

- An authenticated push to github.com or any real GitHub answer (no token, no network by rule).
- The R function on Windows with git present (run-tests is Linux only; the installed Windows application has no
  git). In particular: R on Windows may unset a variable set to the empty string (`Sys.setenv`, through the C
  runtime), which would leave `GIT_CONFIG_VALUE_0` missing so that git refuses the push after the repository was
  created; INFERRED risk for R users on Windows only, not observed. The alternative carrying the same configuration
  in `GIT_CONFIG_PARAMETERS` would avoid an empty variable.
- The Linux systemd launcher (the probe runs the installed Rscript directly; `MSR_BUNDLE` is passed to its units per
  `R/utils_workflow_execution.R:1441-1445`, read, not run).
- A user's own credential helper being bypassed (no helper existed on the runners).
- The full suite (rides the merge train, A4).

## For the change report

**What was wrong, and who was affected.**
- *Destructive, legacy (on `main`):* the Project ID typed in the 'Push to GitHub' panel was used unchecked as a folder
  name under R's temporary folder, and an existing folder of that name was deleted recursively. An ID that climbs
  out of R's temporary folder and names an existing folder, such as `../../../../Documents` on Windows (R's temporary
  folder is four levels below the user's profile), had that folder deleted recursively; R1s S4 observed it with
  `../<sentinel>`. Whether `..` or `.` alone deleted anything was not observed (R's `?unlink` excludes files named
  '.' and '..'). Anyone with a GitHub token who typed such an ID could trigger it,
  on `main` from proteomics, lipidomics and metabolomics, in 0.9 from every omic. `main` 436cdbe8
  `R/func_github_management.R:122-131` (read with `git show main:R/func_github_management.R`); base e0747fc1
  `:123-131`.
- *The token on disk (legacy):* the token was written into the temporary clone's remote URL
  (`.git/config`, base `:323-325`) after the repository was created, and the folder was removed only on success
  (`:362`): after any failed push it stayed, on Windows indefinitely. Observed at G0 on every host.
- *An empty repository created before a local failure (legacy):* the git program was needed only after the repository
  had been created (base `:308-357`); on the installed Windows application (no git) every push with a token and a
  `.Rproj` created an empty private repository and failed, and the name could not be pushed again until it was
  deleted on GitHub. Observed at G0 on Windows and ubuntu:22.04.
- *Undisclosed upload:* the push sends the omic's whole scripts folder, which in 0.9 holds copies of the input files
  (lipidomics, metabolomics, transcriptomics) and the session archive with its data (proteomics); nothing said so.
  Observed tree at G0 (P2g).
- *A control that could not work in the installed application* (no git on Windows, no `.Rproj` made anywhere).

**What changed.** At the fix commits 37af31e4 (Part S) and 9efafdad (product (b)):
- `R/func_github_management.R`: the Project ID validated before any path is built (class
  `multischolar_github_project_id_invalid`); the staged folder removed on every exit (`on.exit`, `expand = FALSE`);
  the git program and git 2.31 or later checked before anything is staged or created (`..._git_missing`,
  `..._git_too_old`); every staged file within 100 MiB before the repository is created (`..._file_too_large`,
  naming the project path); the remote URL without the token, the credential passed to git through
  environment-borne configuration for the push alone with the user's helpers cleared; a failure after the creation
  names `<org>/<id>` (`..._push_failed`). Help and NEWS state what is sent, where and what private means.
- `R/mod_omic_summary.R` (`mod_omic_summary_ui`): with `MSR_BUNDLE` set (the installed application), the panel shows
  one sentence instead of the push controls; elsewhere unchanged.
- No computed result changes.

**Evidence.** G0 runs 45587834788 and 45593577523 (base behaviour in the installed bundles, dummy token, no network
path to GitHub); R1s 45598693322 and R2s 45600718143 (Part S); R1p 45604664357 and R2p 45606131371 (product (b));
W1 45606206919 (the fix in the installed bundles). Limits: Linux test runners only; no real GitHub call; Windows
only as the installed application (no git).
