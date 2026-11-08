# KIT-106 W5: the two W1 proteomics values that D3 changes, re-frozen once

Lead decision (pre-authorised in the W5 brief, 2026-11-09): re-freeze exactly the proteomics `volcano_glimma` values
that D3 turns into `<req>`, once, with attribution. Any other capture difference is a defect or the stop rule.

## The values

Both are in `tests/testthat/test-omic-da-module.R` (`omicDaW1Frozen[["prot_main"]]`), each with a comment citing
KIT-106 W5, D3 and the lead's decision.

| Case, step, field | Before (W1, frozen at b75561de) | After (W5 tree) | Cause |
|---|---|---|---|
| `prot_main` select `values[["volcano_glimma"]]` | `<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06` (the DA-time page) | `<req>` | D3 |
| `prot_main` volcano_unticked `values[["volcano_glimma"]]` | the same page | `<req>` | D3 |

- **D3, the volcano on demand.** Proteomics' UI box is now unticked (`omicDaPluginVolcanoTab`, mod_omic_da_panels.R),
  and the render's first line is `shiny::req(isTRUE(input$volcano_interactive))` (`registerOmicDaPluginVolcanoOutputs`,
  mod_omic_da.R). W1's driver sets the UI's defaults, so the box is FALSE from the start until the `volcano_ticked` step.
  At the start commit the render read no box, so it gave the DA-time page ticked and unticked (W1's note).
- Every earlier step and every other case was already `<req>` (no contrast selected; F106-01), so nothing else moves.
- From `volcano_ticked` on, the value is the same DA-time page (`<ui> 6c24457f...`) as before, at every step.

## Proof that nothing else moved

Scratch runs, not committed: `run1.R`, `load_all`, `MULTISCHOLAR_STRICT_SKIPS=1`, `TMPDIR` on disk, on a copy of
`test-omic-da-module.R` holding only its helpers, its frozen values and the five proteomics (b) cases, with
`MULTISCHOLAR_KIT106_W1_DUMP` set.

- **Start tree** (fe97aac1 plus the W5 tests): 323/0; dump `20dc37c2...` (69,166 bytes).
- **W5 tree, two separate processes:** each 321/2 (the two values above); dumps byte-identical (`cmp`), `c25c8b33...`
  (69,014 bytes).
- **The dump diff** (start against W5), [w1-refreeze-dump.diff](w1-refreeze-dump.diff): exactly two hunks, each the
  `volcano_glimma` value of `prot_main` at `select` and at `volcano_unticked`. Every other field of every step of the
  five cases (outputs, values, kernel calls, notifications, selector updates, `tab_status`, history, state digests,
  the tab's reactive values, rows, files and errors) is byte-identical.

The two values were set by hand from the W5 dump (`"<req>"`), the only change in the dump diff.
