# KIT-111 W6, chain B: install, check, bind, the journeys (M2 and run 7), 2026-11-10

- **Install:** `kit_111` = `browser-smoke.MlQroA`, from a pinned worktree at `c568dfa5` (R code as `341dab60`).
- **`R CMD check --no-tests`:** 3 WARNINGs, 3 NOTEs. Every WARNING/NOTE header equals `kit-108/check.log`'s; the "checking R code for possible problems" section is byte-identical to KIT-108's, so no entry comes from a KIT-111 file.
- **Bindings:** `5096d0fe` (`kit_111` added; 111a..t, 111f, 111q bound).

## The journeys (one attempt each, in sequence, host rule)

| Member | Decision | Exact | Seconds | Dispatches against KIT-108's member |
|---|---|---|---|---|
| 111l | accepted | yes | 139.7 | adds `normalization-ruv`, `pearson-pairs`, `pearson-selection`, each 1/1; nothing else changed |
| 111m | accepted | yes | 250.5 | adds `pearson-pairs`, `pearson-selection`, `protein_ruv_covariates`, each 1/1; nothing else changed |
| 111t | accepted | yes | 136.3 | adds `normalization-ruv`, `pearson-pairs`, `pearson-selection`, each 1/1; nothing else changed |
| 111f | `journey_failed` at `normalisation_ruv` (harness decision `invalid`, failing rule `verifier_accepted`: the expected shape for a failed journey) | — | 69.8 | `normalization-ruv` 1 request, 0 complete; the four QC steps 1/1. Error: "production action failed: … managed action ended: worker_failed". Store head `protein_replicate_filtered` (history initial → … → `protein_replicate_filtered`) |
| 111a | accepted | yes | 196.0 | identical to 108a |
| 111b | accepted | yes | 191.5 | identical to 108b |
| 111c | accepted | yes | 128.9 | identical to 108c |

**M2 (D8's lifecycle test) passes:** the instrumentation sees the moved operations through completion (111l/m/t) and through failure (111f), so M3 (111q) may run.
