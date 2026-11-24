# KIT-205 R2b dump control (run 45404531700; base 8cb536cc, head db09ddfa). Nothing re-frozen.

| Capture | Base dump = stored block | Two head dumps byte-identical | Leaves that differ, base to head |
|---|---|---|---|
| `omicEnrichW1FrozenUi` (`test-omic-enrich-module.R:138-207`) | yes (text) | yes (block SHA-256 c0a5316240a7...) | `html` (digest); `help` (the default rule's note; the table legend); `choices$enrichment_analysis-significance_rule` (the default's label): "measured" to "tested" in each (`ui-base-head.diff`) |
| `omicEnrichW113Frozen[["gprofiler"]]` | yes (text and value) | yes | `provenance/rest` (digest); `record/parameters$enrichment_background` (new: "fixed ; fixed ; tested") |
| `omicEnrichW113Frozen[["kernel_gprofiler"]]` | yes | yes | the same two |
| `omicEnrichW113Frozen[["restore"]]` | yes | yes | the same two |
| `omicEnrichW113Frozen[["organism_filter/no_mapping"]]` | yes | yes | the same two |
| `omicEnrichW113Frozen[["organism_filter/mapped"]]` | yes | yes | the same two |

Instrument note: `dumps-k205.R`'s UI row in `dumps.tsv` reads FALSE because MULTISCHOLAR_KIT107_W1_DUMP also receives the
W1 (b) captures the W113 test drives (9 objects in the file), so the script compared whole files. The UI row above is
computed from the same artifact files, comparing the `omicEnrichW1FrozenUi` assignment only (no rerun). The two head UI
files differ only outside that block, at `omicEnrichW1Frozen[["kernel_gprofiler"]]`'s `gprofiler_plot` digest (line 2009),
a capture this control does not cover. The W113 dump files are byte-identical as wholes.
