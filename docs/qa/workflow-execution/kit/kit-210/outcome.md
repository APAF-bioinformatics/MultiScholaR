# KIT-210 outcome

## Step A: declaration of the runner job (written before it was sent)

- **Purpose.** A1 and A2 of KIT-210 4.1 in one job: re-fetch the pinned RUVIIIC 1.0.19 source (CRAN archive and
  the GitHub mirror `cran/RUVIIIC` at c534151f), hash every file, compare the archives with each other and with
  ticket 2.2, record the installed package's DESCRIPTION; check that the installed compiled library zero-fills,
  and only in the direction step, with its negative control.
- **Script.** `step-a/a1a2.R` (this branch). Network: only the two archives the ticket names.
- **Commit.** The branch commit that adds this file and the script: v0.9-staging c03cb02a plus records only, so
  its lock (`packaging/multischolar.lock.json`) pins RUVIIIC at c534151f, as the runner's library action installs
  it (`packaging/r/install-library.R:33`).
- **Request.** `tests: []`, `base_commit: "none"`, `blas: "reference"`, one script, `timeout_minutes: 40`,
  `out_pattern` `\.(tsv|txt|md)$`, `OMP_NUM_THREADS=1`.
- **Ceiling.** Ticket section 11 row A: 15 minutes expected, timeout 40: 1.67 expected, 3.33 ceiling core-hours.
  Attempts: once (VAL-54 for an instrument defect only).
- **What it decides.** GA1 (archives equal and equal to 2.2) and GA2 (every deciding part of A2 passes), hence
  which of A-I to A-IV holds together with the reading of A3 to A5: A-I needs GA1 and GA2 and no library form
  that handles a missing value without zeros (already shown from the source, 2.2 S1, S2, S8); A-III if GA2 fails;
  A-II if the archives show such a form; A-IV comes from the reading A4, not from this job.

## Step A: result (run 45351569441, 2026-11-24)

**Outcome A-I.** Details, file:line evidence and the algebra in `step-a.md`; the job's files in `step-a/`.

- **Observed.** GA1 pass: the CRAN archive copy and the GitHub mirror at c534151f are identical file by file (32 of
  32) and equal to ticket 2.2; the runner's installed RUVIIIC is 1.0.19 built from c534151f. GA2 pass: on the
  installed compiled library, a missing value of a feature that is neither target nor control behaves exactly as a
  0 (`identical()`), a 20 in its place changes the directions and the output (negative control), and the corrected
  values are reproduced to 4e-13 from RUV-III on the zero-filled cross-product. GA3 not applied: Poulos et al. 2020
  not read (E4 open; no outside service from this machine; the job may fetch only the pinned source).
- **Read (shown).** The library has no other missing-value handling; the package passes NA on all four routes and
  fills nothing before the call; the legacy route chooses k on NIPALS-imputed data and corrects on zero-filled data.
- **Corrected in the ticket's text by this reading.** `inst/config/config.ini:87` (5 NIPALS components) is the
  `[ruvCovariates]` section, not the legacy route's setting: the proteomics Automatic sweep uses 5 (function default,
  `R/func_prot_norm_optimization_helpers.R:679`), the proteomics `ruvCancor` method 2 (`R/func_prot_s4_norm_methods.R:135`),
  the lipid/metab sweep the object's stored argument else 2 (`R/func_metab_norm.R:298`). Bears on C7's ncomp; not on A-I.
- **Answer to the owner** (shown / inferred / conjecture): `step-a.md`, first section.
- **Cost.** 11.8 job-minutes, 0.79 core-hours, against 1.67 expected and 3.33 ceiling (row A). One attempt.

### Not checked (step A)

- Poulos et al. 2020 (A3; GA3 not applied).
- The library's `version = "R"` forms on the runner (shown from source only).
- Whether k = 0 reaches the legacy lipid/metab method in practice (ticket 2.9).
