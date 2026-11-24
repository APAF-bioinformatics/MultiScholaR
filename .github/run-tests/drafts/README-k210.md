# KIT-210 blind lifecycle L1: request drafts (not sent)

Three drafts for `requests/` once the bound commit exists (KIT-210 E1). They live here, outside `requests/`, so this
branch starts nothing. To send one: copy it to `.github/run-tests/requests/<name>-<commit8>.json` on `run-tests`,
replace `BOUND_COMMIT` in `package_commit` and `label` with the 40-digit commit of `kit-210-ruviiic-remedy` on the bound
commit, and push (the brief's request rules).

Blindness (KIT-210 6.9; the check of the workflow is in `kit-210/declaration-draft.md` section 8):
- every full row file and intermediate is written under `k210-private/`, which is never a script's `out`;
- each script that uploads names its `out` and an `out_pattern` that matches only its own enumerated files
  (`status-`, `hashes-`, `checks-`, `timing-`, `provenance-<job>`, `check-e1`, `check-e10`, `check-e1pp`, `check-scale`);
- the scripts that write only private files (sim-k185.R and sim-k210.R with the candidates off, for E1 and E10) have
  no `out`, so only their console (unit counts, start and finish lines) is uploaded;
- the scorer runs inside the job only in its comparison and scale modes; it is never run on L1's rows.
