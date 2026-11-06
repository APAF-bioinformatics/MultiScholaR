# GUI workover scoping (2026-11-06)

Nine read-only scoping passes, one per ticket (KIT-101 to 109), at c4839d25/c783e97a. Each ticket has a "Scoping" section, and the full note is in this folder. Line counts are measured with `wc`. Every "after" figure is an estimate.

## The finding

The intended architecture was a Shiny layer that only observes and displays, over a shared workflow core. The core was built (bindings, declarations, worker, store). The GUI layer was not, for any omic, proteomics included:
- **The GUI:** 55k lines across 164 module files: proteomics 24.0k (75 files), metab 15.0k (40), lipid 14.8k (44), the transcriptomics interim 1.1k (5).
- **Computation in module files:** file parsing, format detection, design validation and construction, ITSD analysis, figure and composite building, session-export assembly, result reshaping, UniProt fetch, FASTA processing.
- **Core-to-module back-references:** the core calls 34 module functions (25 proteomics, 5 lipid, 4 metab).
- **Backend in the GUI process:** proteomics normalisation, correlation, enrichment and the DA limma kernel run in the app process; lipid, metab and trx send them to the worker.
- **Proteomics is the pattern only in parts:** its QC sub-tab spec tables, its one-call step functions, and the DA run handler shape. Its import, design, normalisation and enrichment are not light.

## Size by family (lines, before → estimated after)

| Family | Before | After (estimate) | Note |
|---|---|---|---|
| Shell and import (102) | 9,485 | about 3,900 | ~1,830 computation lines move to the core |
| Design builder (103) | 9,633 | about 4,500–5,000 | every design-table algorithm is the same across omics; run order stays per omic and recorded |
| QC (104) | 5,888 | about 1,500 | lipid and metab are one design written twice; proteomics joins through its spec tables |
| Normalisation (105) | 6,774, plus proteomics 4,505 | about 2,100–2,300, plus proteomics | trx has its own panel set; proteomics largely separate |
| DA (106) | 7,038 | about 3,600 | proteomics shares the UI and renders; its run is in process |
| Enrichment (107) | 4,982 | about 2,800–3,300 | omic-agnostic; trx later as configuration; lipid and metab need new kernels |
| Summary and report (108) | 4,596 | about 880 | behaviour differences kept as per-omic settings |
| Removal and closure (109) | — | — | 41 seam-test files quarantine candidates; skips compared per file |

## Decisions collected for the user

Each scoping note lists its own. The cross-cutting ones are below, with the coordinator's recommendation.
1. **Proteomics joins the shared layer:** UI frames and panels shared, its route-specific pieces as plugins, through a DIA binding adaptor. *Recommended.*
2. **Proteomics' in-process backend** (normalisation, correlation, enrichment, the DA limma kernel) moves to the worker, as its own tickets, measured. *Recommended* under "no backend in the frontend".
3. **Base:** new shared omic-neutral files, rather than generalising the lipid modules in place. *Recommended.*
4. **Exported module functions** stay as thin wrappers (public API). *Recommended.*
5. **Per-omic tab ids, module ids and status keys** kept in the configuration, so the journeys and e2e lanes stay valid. *Recommended.*
6. **Old copies are deleted in each family ticket** once nothing calls them, with their seam tests quarantined there. This avoids an interim size increase. *Recommended over leaving it all to KIT-109.*
7. **Behaviour differences outside the parity list:** keep them per omic, except three to unify:
   - metab import opens the design step as the others do;
   - metab's design-save failure becomes a notification, as lipid's is;
   - the interactive volcano is off by default and built on demand for every omic.

   Kept per omic: the summary's state choice, record writing and copy-failure handling; lipid's zero-row refusal; assay-name rules.
8. **Split KIT-105** into 105a (core moves plus the lipid/metab frame) and 105b (transcriptomics).
9. **Test infrastructure in KIT-101:**
   - committed fixtures for the entry-point suite;
   - an offline UniProt double;
   - a guard sentinel that turns the silent file-level skips (13 files, 490 tests) into failures;
   - a per-file baseline of tests run and skips.
