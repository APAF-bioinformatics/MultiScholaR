# KIT-183 outcome

Ticket `tickets/kit/KIT-183.md` revision 1.4 (with the coordinator's text edits; VAL-58 (4) superseded by VAL-62;
VAL-64). Branch `kit-183-integration-test-of-record` from `v0.9-staging` 1e627f70; tested head d53a5dd1 (records only
after it). Every run: `runs.tsv`. Simulation: `declaration.md` (committed alone at fae77e5f), `results.md`, `s/`.

## Observed

- **Simulation S (45225854502): outcome O1.** E1 holds at n = 12 (largest upper limit 0.0737) and n = 20 (0.0711);
  E2 holds for gsea (0.0599) and camera_pr (0.0637); E3 holds (33,877 and 34,891 called planted sets, all in the
  expected direction); 2,400 of 2,400 replicates completed; determinism identical. `camera` is the test of record for
  every data set; `camera_pr` is not offered (the code was written to O1).
- Gates: F1 passes (limma 3.68.4 = 3.68.5 in every function read); F3 passes (one seed = five); F4 read (MOFA2 1.22.0).
  L1: three runs, the first two instrument failures of the script's last aggregation (VAL-54); 57 to 62 s per replicate,
  above the 50 s stop rule; re-split accepted (VAL-64); 16 class-size cells below 90 % formed, redefined in the
  declaration's dated section before S. In S, metabolomics C2 size 10 at n = 20 formed 62 % (reported).
- Failing first (R1a2, cab55462, R/ equal to 1e627f70): every new test of behaviour fails or errors on the base code.
- Head against base, integration files (R1b) and the full suite (F, 45237273732, d53a5dd1 against 1e627f70): 447 files
  the same; `test-integration-enrichment-mofa.R` one more skip (T14 without MOFA2); the new file 176 passed and 2
  skipped (T5, T8: their store exists only in the base-store job); no new failure; 39 files fail identically on both
  sides (stored laptop oracles); R CMD check 0 ERROR, 3 WARNING, 3 NOTE on both, nothing only on head. Gate-off pages
  identical to base (R1b).
- Base store in one job (R2b, VAL-64 (2)): the store made by the base code and read by the head in the same job.
  T5 passes (5 of 5: a base KEGG entry replays with gsea, terms identical, its test restored as gsea); T15 passes;
  **T8 fails on head** (0 of 3): the ranking and the gsea product are recomputed, not hits. Cause, read in the code:
  every saved product's scientific identity includes the build's digest (`context$software`,
  `R/utils_workflow_dia_state.R:1655-1665`: the installed binary's SHA-256; under a source load the version and
  `body(newDIAStateAdapter)` with its source references), so a product saved by another build is never a cache hit.
  T8's premise (a base-made product is a hit at head) cannot hold for any change. Mutations (X): m1 (restore removed)
  makes T5 fail; m3 (a changed default of `integrationEnrichmentGsea`) makes T15 fail; m2 (test declared) makes T8 fail,
  but T8 already fails on head, so that is not evidence.
- MOFA2 lane (M, 45237241828): T14 green; camera through the store on a trained model, recorded, replayed offline with
  stand-ins that fail on any call: ranking and terms identical.
- Bundle lane (B, 45237281880, stage suite): the new test file on the installed bundle, ubuntu-22.04 and windows-2025:
  176 passed, 0 failed, 2 skipped each; the cross-platform test within 1e-6 of the Linux reference on both
  (`cross-platform.tsv`). ubuntu-24.04 not run (the suite stage has no 24.04 job).
- T13 (corpus): planted sets found in the expected direction, decoys not; of 40 invented pathways neither test called
  any at BH <= 0.05 on proteomics or metabolomics_GCMS (the measured size of the change, in NEWS).
- Memory case (8.7), not met as worded (VAL-85 (5)): 2,000 sets on 20,000 features x 500 samples in 9.5 to 10.3 s;
  process peak RSS 1,176 MiB (the ticket's 1 GB peak not met as worded); the high-water mark before the case was
  592 MiB, so the case added at least 584 MiB.
- Budget: 2,065.9 job-minutes, 137.7 core-hours by GitHub job times, against 146.36 expected and 262.67 ceiling. The
  L1 item used 4.21 core-hours against its 3.0 ceiling (three runs, two instrument corrections, VAL-54 (iii)); the
  ticket total holds. Module CI runs triggered by the pushes are not counted here (as in earlier tickets' records).
- Manuscripts deliverable (ticket 12: the Methods sentence, the limits sentence and the simulation's summary figures
  with the four statements of 15.2): owed, not done on this branch (drafts go to ~/Projects/APAF/Manuscripts).
- Release rule (DISPATCH-STATE, KIT-184 C02): B's bundles were built from d53a5dd1 and discarded inside their jobs; none
  left development. From this merge until KIT-184 merges, no release or demo installer may be built from staging.

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A1 | met: F1 to F4 recorded; L1's pool counts and stop rule applied (re-split, VAL-64; classes redefined) | `feasibility.md`, `lifecycle.md` |
| A2 | met (T1, T2) | F |
| A3 | met (T3, T6, T7) | F |
| A4 | met (T4) | F |
| A5 | met: T5, T8 (restated, VAL-85 (3)), T15 in the one-job guard run on the merge commit | G, `g/` |
| A6 | met (T10) | F |
| A7 | met (T11; T12 with its base comparison) | F |
| A8 | met (T13; the size of the change in NEWS) | F, R1b log |
| A9 | met (T14) | M |
| A10 | met: declaration alone before S; O1; classes formed or redefined before S (one cell at 62 % in S, reported); determinism identical; comparators reported; texts filled | `results.md`, S |
| A11 | met | R1b scripts job |
| A12 | met: within 1e-6 on the ubuntu-22.04, ubuntu-24.04 and windows-2025 bundles; magnitudes and SHA-256 not observable there | B, B24, `cross-platform.tsv` |
| A13 | met: failing first shown (R1a2); the mutations caught: m1 fails T5, m2 fails T8 (a) and (b), m3 fails T15 | R1a2, G |
| A14 | met: no new failure, check findings unchanged | F |
| A15 | met: ledger rows 1 and 2 (`VAL-<merge>`) | `docs/methods-decision-ledger.md` |

## Inferred

- T8's failure is in the ticket's premise, not in the code: the `integration_enrichment` and ranking declarations are
  unchanged and T5's replay is identical; a same-build cache hit is shown by T9 (camera) and KIT-141's test (STRING).

## Not checked

- ubuntu-24.04 (B); the SHA-256 of the bundle platforms' tables (the lane uploads nothing).
- A cache hit across builds (impossible by design, above).

## After the implementation audit (VAL-85)

- Texts, NEWS, results, ledger and records repaired as AUDIT-KIT-183-impl I01, I02, I04 to I10 give them; NEWS names
  KIT-191's further options.
- T8 restated (a), (b), (c); guard run G (45273640961) on the merge commit 912cd36b: all pass at head, each mutation
  caught. The earlier T8 (R2b) failed because of its premise, kept in `r2b/` with an annotation.
- ubuntu-24.04 bundle (B24, 45273632957): 178 passed, 0 failed, 2 skipped; the cross-platform test within 1e-6.
- Merged v0.9-staging 3ff60190 (merge commit 912cd36b; NEWS and ledger mechanical, the rows stay `VAL-<merge>` after
  VAL-86); the integration test files on the merged tree (G): no failure beyond the stored oracle failing on base too.
