# F107-05 triage (KIT-107 W3)

**Question.** On W1 (b)'s `kernel_*` cases (the live `processEnrichments` over the fake `gost`/`enricher`), the
`enrichment_completed` commit is refused with "payload contains active or lazy bindings" while the run reports complete.
What holds the binding, and does it come from the test fakes or from what the real kernel returns?

**Probe** (scratch `kit-107-w3/probes/f107_05.R`, start tree `eacd4a5b`): the commit's payload (`saveDIAOwnerState`'s
`s4_data_object`) walked as `checkSavedWorkflowPayload` walks it, reporting every environment with an active or lazy binding
and every live pointer.

- `kernel_gprofiler`: 17 environments, all under `results@enrichment_plots[[contrast]][[up|down]]`, the kernel's static
  ggplot. Path: the layer's ggproto `super` closure, its environment (ggplot2's `layer()` frame), `user_env` (the frame that
  built the plot: `generate_enrichment_plots`, `func_general_enrichment_engine_helpers.R:199-280`, which also holds
  `interactive_plot <- plotly::ggplotly(static_plot)`), then the built plot's scales and guides: ggproto method closures
  whose `self` is an unforced promise, and a scale transform's `reverse`. No active binding.
- `kernel_clusterprofiler`: 19 environments, the same kinds (`self`, `reverse`), under the static ggplot the kernel builds
  inline (`processEnrichments`, `func_general_enrichment_orchestration.R:690-740`, then `ggplotly` and `ggsave`). After
  forcing them, the check refuses one more object: through `user_env`'s parents, the `processEnrichments` frame's
  `da_results@da_data[[contrast]]` (the DA table `createDAResultsForEnrichment` read with readr) carries a live parser
  pointer in its `problems` attribute ("payload contains a nonportable resource").

**Cause: real in live use.** Every object involved is made by the kernel with ggplot2, plotly and readr: the static plots
and the frames they capture (ggplot2's layer records its caller's environment), the ggproto methods the build leaves, and
the DA table readr returns. The fakes only supply the services' result tables; they hold none of these bindings. With real
gprofiler2/clusterProfiler results the kernel builds the same plots in the same frames, so a live enrichment run never
saved `enrichment_completed` (the save's warning is console-only; the tab reported success). The deterministic double
stores no ggplot (its `enrichment_plots` are NULL, its plotly figures are `plot_ly`), which is why the test-mode lanes commit.

**Fix (core, no method change).** `portableWorkflowEnrichmentResults` (`R/utils_workflow_enrichment.R`), called by
`saveWorkflowEnrichmentState` immediately before the provenance and the commit: every lazy binding reachable from the
results is forced (it then holds the value it would have read), and every readr table's diagnostics are materialised
with the existing `prepareDIAParserTable`, in place. No value, method, parameter or result changes; the portability check
itself is unchanged and still refuses anything the helper cannot make portable (an active binding, a binding whose
promise fails). Probe after the fix (`probes/save_warning.R`): `saved` TRUE, no warning, for `gprofiler` (the double: a
no-op) and `kernel_clusterprofiler`.

**Tests.** `test-omic-enrich-module.R` "(W3, F107-05) a Run on the live kernel (…) commits enrichment_completed with the
tab's results", for both kernel cases: the state commits, is an `EnrichmentResults`, and its `enrichment_data`,
`enrichment_summaries`, `contrasts` and plot classes are identical to the tab's results. Failing on the start tree, passing
on the W3 tree (`test-first.txt`, `targeted.txt`). The two kernel cases' W1 history, state and store-count values are
re-frozen once (`w1-refreeze.md`).

**Not changed:** the kernels (`func_general_enrichment_*`), `createDAResultsForEnrichment`, the saved-entity check.
Recorded for KIT-113 (D10's memory note): the committed state now carries the static plots' captured frames, as the check
would always have allowed had they been portable.
