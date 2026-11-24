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
