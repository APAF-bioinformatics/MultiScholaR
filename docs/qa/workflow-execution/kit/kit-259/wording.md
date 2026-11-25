# KIT-259: user-facing text (drafts for the owner's pass, VAL-41)

Every text below is a draft. The owner sets the final words. "Where" names the function that shows it.

| Id | Where | Draft text | Note |
|---|---|---|---|
| W1 | `buildProtNormFigureUnavailableNotice` (warning notice, one per run; in process and managed) | "The normalisation finished, but the figure '<figure>' could not be drawn: <reason>. The normalised data are not affected. In this tab, the composite QC figure and the report the figure is replaced by a note giving the reason." | Ticket section 7. `<figure>` is the figure's name, e.g. "Post-normalisation Pearson correlation" (`protNormFigureTitle`). `<reason>` is the drawing error's message as R gives it. |
| W2 | the same, when more than one figure failed | "The normalisation finished, but these figures could not be drawn: '<figure 1>' (<reason 1>); '<figure 2>' (<reason 2>). The normalised data are not affected. In this tab, the composite QC figure and the report each figure is replaced by a note giving the reason." | Not in the ticket: the plural of W1 (one notice per run, ticket section 3). |
| W3 | `writeProtNormFigurePlaceholder` (the placeholder image) | "<Figure>: could not be drawn." followed by the reason, wrapped | Ticket section 7. |
| W4 | `buildProtPearsonUnavailablePlot` (proteomics histogram, no within-group pair) | "No two samples share a value of '<Color by>' (Pool and QC groups are left out), so there are no within-group pairs to correlate. Choose a 'Color by' with at least two samples per value." | Ticket section 7 (A2). |
| W5 | the same, when within-group pairs exist but none has a correlation | "None of the <n> within-group sample pairs has a correlation (too few shared values, or a sample with no variation), so there is nothing to draw." | Not in the ticket. A-E5 (every correlation NA) is covered by A2, but W4 would be false there (pairs exist), so this case has its own sentence. |
| W6 | lipidomics and metabolomics `plotPearson` subtitle (only when a pair has no value) | "<n> sample pairs shown; <k> have no correlation (too few shared values, or a sample with no variation) and are not shown." | Ticket section 3 (A3) as amended by Re-check 1 R1-03. "1 sample pairs" is not made singular: the owner's pass. |
| W7 | possible plain translations of ggplot2's reasons | e.g. for 16 or more groups as point shapes: "more than 15 groups to show as point shapes; choose another 'Shape by'" | Ticket section 7, last bullet: not implemented; the reason is shown as R gives it. For the owner. |
