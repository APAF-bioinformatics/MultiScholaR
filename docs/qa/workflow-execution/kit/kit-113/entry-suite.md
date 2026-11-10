# KIT-113 entry-point suite

- **Simulation Setup (before run 3):** see [proof/chain-a.md](proof/chain-a.md). `prot_dia`, run alone with traces, is identical to its oracle.
  - Called once (one Run): the runner, its helpers, `runDIAReferenceEnrichment` and its two halves.
  - Called 0 times (in-process session): the managed request half and the worker handler.
  - Trace: [proof/sim-prot_dia.txt](proof/sim-prot_dia.txt).
- **Run 3:** the four files ran in parallel under strict skips.
  - trx 29, lipid 35, metab 105 and proteomics 193 passed;
  - 0 failed, 0 errors, 0 skipped;
  - all ten lanes show 0 step differences against their oracles.
- **No re-freeze.**
