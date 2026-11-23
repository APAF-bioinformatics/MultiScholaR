# KIT-177 wording (ticket section 14, D3; implemented as proposed, VAL-41; the owner edits before release)

## In the app (`omicGuiProtLog2Texts()`, `R/utils_omic_gui_config.R`)

- Recorded transform: `log2(x); 0 missing`.
- Reopen notice (type "message", once): "This design was saved by an earlier version, which added 1 to every value
  before taking log2. Your saved results are unchanged. If you save the design again, values are log2-transformed
  without adding 1, and later results will differ slightly, most for small values."
- Counts notification (type "warning", each clause only above 0): "n values of 0 [and k values that were not a
  number] were read as not measured (missing). They are not used in the analysis."
- Step log: "Values equal to 0 set to missing: n"; "Values that were not a number set to missing: k".
- Refusal (`prot_design_quantity_invalid`): "Run '{run}' contains negative or infinite values (for example {values}).
  Quantities must be zero or positive; correct the file and import it again."
- `scale` warning (type "warning"): "Scale normalisation divides each sample by the median of its absolute log2
  values. In this data, [n] values are below 1, so their log2 values are negative, and taking absolute values treats a
  log2 value of −2 like +2. The result then depends on the unit of your table. Cyclic Loess and Quantile do not have
  this problem." Step log: "Scale normalisation: log2 values below 0 in the input: n".

## Methods (for the manuscript; drafts go to ~/Projects/APAF/Manuscripts)

- Label-free: "Protein intensities were log2-transformed; a value of zero was read as not detected and was not imputed,
  so a protein was tested only where it was measured, and a fold change can be underestimated when values fall below
  the detection limit in one group."
- TMT: "Reporter abundances were log2-transformed; values missing in a plex were not imputed."
- Both (limitation), with R0's numbers (VAL-80; 20 simulated data sets per setting): "No variance-stabilising transform
  was applied. With the default test (limma's moderated t-test with TREAT at a log fold-change threshold of 0), each
  protein's variance is moderated toward one prior variance shared by all proteins; no abundance-dependent trend is
  fitted, so the larger variance of log values at low abundance is not modelled. In simulations with typical
  intensities near 10, fold changes in the lowest third of abundance were shrunk to 0.75 of their size (0.55 with
  log2(x + 1)), 8.4 % of true changes there were detected (14.5 % with log2(x + 1)) and 71 % of the calls there were
  false (65 %), with RUV at its automatic setting; from typical intensities of about 150 on, the two transforms gave
  the same results."
- The test is described as run (TREAT statistics without trend or robust; `R/func_prot_da_model_stats.R:306-313`), not
  as configured (`eBayes_trend`, `eBayes_robust`, `inst/config/config.ini`).

## NEWS

`NEWS.md`, `## MultiScholaR (unreleased; numbered at merge)`, `### Scientific changes`: the D3 text with 9.1's measured
numbers and F2's finding ("some Proteome Discoverer exports report signal-to-noise values").
