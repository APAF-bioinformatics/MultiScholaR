# KIT-108 W3: the one W1 re-freeze (F108-03)

**Authority:** the lead's decision on F108-03 (`kit-108/findings.md`): after a refused integration write there is no
success notice. Assigned to W3; this is the only W1 capture value that moves.

**The value:** `tests/testthat/test-omic-summary-module.R`, `omicSummaryW1Frozen[["trx_main"]]$save_before_export$notes`.

| | Value |
|---|---|
| Frozen at the start commit | `"[warning 10 -] The integration export needs the RUVIII export: run it in the panel above, then save again."`, `"[message 5 -] Saved Integration S4 Object"`, `"[message 5 -] Study parameters saved successfully"` |
| Re-frozen at W3 | `"[warning 10 -] The integration export needs the RUVIII export: run it in the panel above, then save again."`, `"[message 5 -] Study parameters saved successfully"` |

Every other field of that step (the session summary's `Integration Object: not written (no integration export yet)`,
`workflow_args_saved = TRUE`, the calls, the files, no integration file) is unchanged. The code: the core's
`saveWorkflowSummaryParameters` shows `integration_saved` only when `startWorkflowSummaryIntegrationExport` reports the
file written (`done(TRUE)`); the gate's refusal reports `done(FALSE)` (R/utils_workflow_summary.R).

**Test first:** `(W3, F108-03) trx: the save before the RUVIII export ...` (test-omic-summary-module.R) failed on
c2702c21, where the wrappers were still the start commit's modules (`kit-108/w3/test-first.txt`), and passes on the W3
tree. `test-omics-summary-core.R`'s `(W3, F108-03)` block tests the entry directly.

**Nothing else moved:** two dump runs of the whole W1 file on 1b3b9205 (the wrappers on the frame), in separate
processes (`MULTISCHOLAR_KIT108_W1_DUMP`; 2026-11-09, no `multischolar*` unit active):

```
dump A: COUNTS test-omic-summary-module.R pass 1702 fail 1 error 0 skip 0 (the one failure: (b) trx_main), 192 s
dump B: COUNTS test-omic-summary-module.R pass 1702 fail 1 error 0 skip 0 (the one failure: (b) trx_main), 187 s
cmp dump-A.R dump-B.R: byte-identical, 270,388 bytes, sha256 6fa8a11986133c28d2e7fbb940b9f802a438dbbce444c3a91c7596f620a6b55f
```

`diff <frozen text between the markers> dump-A.R` before the re-freeze (the frozen text is the UI block and the cases,
without the `omicSummaryW1Frozen <- list()` line; the two `a` hunks are the dump's blank line after its last entry of
each block, which the marker extraction drops):

```
141a142
>
1802c1803
<             "[message 5 -] Saved Integration S4 Object", "[message 5 -] Study parameters saved successfully"),
---
>             "[message 5 -] Study parameters saved successfully"),
3715a3717
>
```

After the edit (the two comment lines above the value left out), the diff is the two blank-line hunks only. The
re-frozen file's run: `kit-108/w3/targeted.txt`.
