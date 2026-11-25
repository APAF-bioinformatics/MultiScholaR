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
