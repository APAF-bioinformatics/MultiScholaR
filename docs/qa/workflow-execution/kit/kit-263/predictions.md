# KIT-263: predictions and the branch-point re-check (written before any run)

Branch `kit-263` from `origin/v0.9-staging` 127011429f5160e6a1ed57fd0552dc16a90d5aa9. Ticket r1.1.

## r1.2 note: re-read of every cited location at 12701142

Between 0a774e0b (the plan audit) and 12701142 only `docs/methods-decision-ledger.md` changed (`git diff --stat`).
No cited line moved: the writer's signature `:12-17`, the config_list resolution `:83-91`, the merge `:151-167`,
the Summary core's call `:264-266`, the recorder `test-omic-summary-module.R:451-458` (mocked at `:487`), the Rd usage
`:6-17`. The filter goes after `:91` (the `config_list` log line), before the merge at `:151`; nothing between reads
`s4_params` or `config_list` except `:643` (`buildWorkflowRuvOptimizationLines`, RUV sections, none foreign).
The ticket's grep recount: 44 files (41 test files + 3 helpers), the same 41 test files as the ticket's list; with
`test-gui-entry-proteomics.R` and the new file, 43 files.

## Predictions

R1 (tests-only commit, base = 12701142):
- `test-kit263-prot-study-parameters-sections.R`: T1 fails (6 foreign lines against 0), T5 fails (`NULL` against
  "transcriptomics"), T4 errors (`studyParameterForeignSections` not found), T2 and T3 pass.
  Stop rule 1: T1 at 0 on base.

R2 (fix against tests-only commit; 43 files + R CMD check; env MULTISCHOLAR_TRX_MANIFEST_NOTE=1):
- every file `same`, except the new file (`only-on-head`; all five tests pass on head);
- `test-omic-summary-module.R`: `same`, its failures identical on both sides (merge train 11: 1675 passed, 28 failed);
- `test-gui-entry-proteomics.R`: `same` (failures, if any, identical on both sides);
- R CMD check: identical on both sides (the Rd matches the roxygen).

## Diagnostic run D1 (declared 2026-11-26, after R1 run 45676487966; printing only, no expectation, no file under R/)

R1's result stands: T1 3 foreign lines (predicted 6), T2 3 against a computed 6, T3's second assertion absent. The
diagnostic prints (A) the T1 fixture's whole file and the sections each block iterates over, and (B) the file the
Summary save writes for a proteomics project on the module-ci bus (B1 as the fixture is; B2 with the shipped
configuration as `config_list` and `@args`, as the design import makes them). Classification, committed before the run:

1. (B2) shows the foreign block ONCE, in "Parameters from Final S4 Object" only: the ticket's "printed twice, six lines
   today, ten after KIT-217" was a wrong reading of the base. The ticket is amended to the observed counts (T1's base
   value, T2's computed count, T3's second assertion, "twice" in sections 2, 5 and 11, NEWS, the change-report counts),
   recorded as a declared-comparison correction from a diagnostic observation; the fix is unchanged (it filters both
   sources); the coordinator corrects VAL-279's "twice" by a correction row.
2. (B2) shows it TWICE while (A) shows it once: the T1/T2/T3 fixture does not represent the application. Instrument
   defect: the fixture is rebuilt on (B)'s construction with the expectations unchanged, and R1 is re-sent once
   (instrument rule, 1 of 3).
3. Any other pattern: stop and report it as it is.

## Diagnostic run D2 (declared 2026-11-26; the last diagnostic; printing only, no expectation, nothing under R/)

D1 (run 45679194881): fixture A printed the foreign block once (S4 block), the Additional block empty because
`formatConfigList` errors on `plotRle`'s length-2 value (separate defect, recorded in outcome.md, not fixed here); the
module-ci construction gave the writer no S4 object. D2 observes the application's S4 block through a proteomics save
with an owner (the `omicSummaryW1ProtSetup(owner = TRUE)` construction, protein object `@args` and bus `config_list` both
the shipped configuration), printing the whole file with line numbers and the sections each block iterates over; and a
lipidomics owner's save, headers only. Note: with an owner the writer uses the analysis record's formatter for the
Additional block (not `formatConfigList`), so that block may print.
Classification, committed before the run (the coordinator's):
(i) the foreign block appears once, in the S4 block, the Additional block empty: the ticket is amended to the observed
    state (T1's base value the observed count; T2's count for one block from the shipped config; T3's second assertion
    replaced by one on the S4 block; a guard row recording, not asserting as desired, that the Additional block is empty
    on base and head for the separate defect; one test calling the filter on the config_list path directly with a config
    `formatConfigList` can format; "twice / six / ten" corrected), fix unchanged.
(ii) twice in the application: the ticket's counts stand; the T1 fixture is rebuilt on D2's construction (instrument
    rule, R1 re-sent once).
(iii) no foreign line at all through this path: STOP and report; VAL-279's premise needs the owner.
Lipidomics: no proteomics section expected.
