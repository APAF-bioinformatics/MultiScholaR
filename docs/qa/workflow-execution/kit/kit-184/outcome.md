# KIT-184: outcome (running record)

Ticket `tickets/kit/KIT-184.md` revision 1.4. Branch `kit-184-ruv-factors-in-integration`, from `v0.9-staging`
67d38df3 (KIT-183, KIT-187 and KIT-191 merged). All computation on GitHub-hosted runners; no R on the laptop.

## State

| Step | State |
|---|---|
| Step 0: ticket text (D01 to D04; K7 per VAL-85 (3); KIT-191 8.2; citations) | done, revision 1.4 (ticket section 20) |
| Gates G0, G5 (reading) | G0 passed; G5 cannot occur (`feasibility.md`) |
| Tests-only commit a5229552 and its run R1a | committed; run 45321879760 sent |
| The fix feb3d6f5 | committed |
| R0 (base-store job), R1 (base against head, check, scripts), X | not sent yet |
| L1, declaration, S (simulation) | not started |
| M (MOFA2 package lane), B (bundle lane), F (full suite) | not started |

## Commits

| Commit | What |
|---|---|
| a5229552 | tests first (K1 to K16; nothing under `R/`) |
| feb3d6f5 | the fix: writer, reader, import, resolver, kernel, key, replay, record, Methods, tab, workbook, Rd, starter, NEWS, ledger row |
| 4ff76465 | records: base-store capture and job (with run X's mutations), the writer's cost measurement |

## Where the tests are (differs from ticket section 10 in location only)

K1 to K3, K5 to K12 and K14 to K16 are in the one new file `tests/testthat/test-integration-removed-factors.R` (with its
own fixtures), as KIT-191 kept its U1 to U13 in one file; the ticket places K6, K7, K14 and K15 in KIT-183's file and
K10 to K12 in the record, tab and workbook files. K4 is in the per-route writer files through one helper,
`expectIntegrationRemovedFactors()` (`tests/testthat/helper-integration-input.R`). K13 is in the MOFA2 lane's file.

Existing tests restated in the tests-only commit because the ticket changes their behaviour by design:
- `test-integration-enrichment-tests.R`, KIT-183 T10: a result of a test with a design names its form in its row.
- `test-integration-enrichment-tests.R`, KIT-183 T11: the RUV note beside the selector stays with the tests that have
  no design.
- `test-integration-input-matrix.R`: the file list of an input with a corrected matrix gains `removed_proteomics.parquet`.
- `test-integration-record.R` (`intRecRewrite`): its corrected matrix equals the normalised one, so the writer's classed
  warning is expected and silenced.
- `test-integration-gui.R` (the scan of a saved input): the same fixture; the classed warning is asserted.
- `test-workflow-integration-input.R`, the stand-in test of the classic route after the correlation filter: its states
  differ by a constant per protein, so it now expects `not_determined`, the classed warning and no factors file.

K4, not built as the ticket words it (recorded, not silently dropped):
- the lipidomics declared design "with and without a missing value in the fixture": one fixture only (the existing one);
  the principal angle to the DA's stored factors is printed for it;
- the transcriptomics test prints no angle to the state's RUVSeq factors: it asserts instead that the written factors
  span RUV-III's own W recomputed from `fullalpha`, which the RUVSeq factors do not.

## Not checked (so far)

- Everything that needs a run: no result has returned yet.
- The tab in a browser (the module is driven with `shiny::testServer` and `renderTags`).
- The numbers [R3], [I3], [R5], [I5], [R5m], [R6], [I6], [PF], [PC] of the texts: placeholders until the counted
  simulation. No build may be released with them in place.
