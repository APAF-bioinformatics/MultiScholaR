# KIT-264 wording (drafts for the owner's pass, VAL-41)

No sentence of the report changes: the stale cases now print the sentences the never-exported case prints today.

1. NEWS (`### Fixes`, as committed): "Transcriptomics report and analysis record: a RUVIII integration export made from an
   earlier normalisation (before the normalisation was run again, with another k or without RUV) is no longer described
   as part of the analysis. Before, the report still said RUVIII had been applied, gave the earlier k and showed its
   figure, next to a caption saying no correction had been applied. The earlier export figure is now removed when you run
   the normalisation again, so the results package no longer carries it. A report of an analysis saved with an earlier
   version keeps the earlier statement until Save Workflow Arguments is pressed again. No result changes: with RUV applied
   the integration input already refused such an export. This was in the first 0.9 release (v0.9.0). (KIT-264)"
2. The run's log line (existing text, now also naming the earlier export's figure): "Removed the earlier run's files this
   run replaces: ruviii_before_and_after.png, ruviii_before_and_after.pdf". The figure is the earlier export's, not the
   run's: accept as is or reword.
3. Open for the owner: whether the report should tell the reader that an earlier export exists (ticket option (d), not
   taken), and X1 (an export after a run without RUV): see outcome.md.
