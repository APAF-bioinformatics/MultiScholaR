# KIT-112 entry-point suite

- **Simulation Setup (before run 3):** see [proof/chain-a.md](proof/chain-a.md). `prot_dia` and `prot_dia_limpa` were each run alone with traces, and both are identical to their oracles. The runner, request, compute and accept halves were each called once (one contrast per lane). The worker handler was called 0 times, as expected in an in-process session. The traces are [proof/sim-prot_dia.txt](proof/sim-prot_dia.txt) and [proof/sim-prot_dia_limpa.txt](proof/sim-prot_dia_limpa.txt).
- **Run 3:** the four files ran in parallel under strict skips. Results:
  - trx 29, lipid 35, metab 105 and proteomics 193 passed;
  - 0 failed, 0 errors, 0 skipped;
  - all ten lanes have 0 step differences against their oracles.
- **No re-freeze:** KIT-112 changes no lane value.
