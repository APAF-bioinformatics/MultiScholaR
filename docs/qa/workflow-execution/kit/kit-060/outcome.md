# KIT-060 outcome: decisions for the user, the re-estimate, and what remains

Audit of 2026-11-17 at `015d6577`. Nothing was run, installed or committed. The register, the inputs contract, the GUI notes, the corpus specification and the MOFA2 source check are in `AUDIT-integration-touchpoints.md` (same folder); section and row numbers below point into it.

## 1. The main findings

1. **Each omic already writes an "integration file", and the four are not alike.** Proteomics writes its whole session archive (`.msdia.tar`), lipidomics and metabolomics write a whole S4 object, transcriptomics writes one TSV with samples in rows. Three of them carry a free-text label in the name. None says whether its matrix is corrected or only normalised, and none is hashed in the record (section 2.1; rows I2 to I4, I9, I10).
2. **The workbook reads from a folder no omic writes to**, under a file name the package does not produce (`ruv_normalized_...` against the package's `ruv_normalised_...`) (rows I1, I2).
3. **"Corrected" names five matrices, not four.** Transcriptomics has two: the `removeBatchEffect` display matrix and a separate RUVIII export made for integration, which is the one the workbook reads. Revision 1.3 names the first (section 2.2).
4. **Transcriptomics cannot run with RUV off**, and its normalised matrix (TMM logCPM) is never written to the project (section 2.1).
5. **On DIA_LIMPA in its default mode the corrected matrix is a side product**, not the state the Summary exports (section 2.2, row C2).
6. **No key ties a sample in one omic to the same biological sample in another.** The workbook aligns by sorting sample ids and matching positions, and removes the dropped sample by an id typed in the code, with two different ids in two places (rows A1, A2).
7. **The shipped workbook does not train a model.** The training chunk is commented out and a ready `model.RDS` is read instead, so there is no run to be exact against until the verified workbook arrives (row M5).
8. **The per-feature Z-scaling is the workbook's own step.** MOFA2's documentation speaks of scaling views, not features (section 7.2).
9. **The GUI frame can host integration by name with little change** (the dispatch, the folders and the home tile exist), but the shared shell's validator requires a single-omic pipeline's groups, and the import scan, the alignment and the model tabs are new (section 5).
10. **Several downstream functions read global variables** (`plotMofaWeights`, the metabolomics enrichment helpers) and cannot run inside a module as they are (rows DS4, E4).

## 2. Decisions for the user

Each is stated with what the audit found and, where there is one, a proposal. A proposal is not a decision.

| # | Decision | What was found | Proposal |
|---|---|---|---|
| D1 | Where each omic's integration input lives. | The ticket's point 4 says "under its own results folder". The code already writes to the shared `<project>/integration/`. The workbook reads `data/integration/`. `results/integration/multiomic_inputs/` exists unused. | `<project>/integration/<omic>/`: the place the Summary tabs already write to, one subfolder per omic. This departs from the ticket's wording and needs the user's word. |
| D2 | The form of the input, and what happens to today's four files. | A tar archive, two S4 objects, one TSV. | Plain tables with a JSON manifest (section 4.1), fixed names, no label. Today's files either stay beside it for compatibility or stop being written; the user chooses. |
| D3 | Which transcriptomics matrix is "the corrected matrix". | C3a: `removeBatchEffect` on TMM logCPM with the group protected; display only today. C3b: `ruv::RUVIII` on TMM logCPM with the group (default) or technical replicates as replicate structure; it is what the workbook reads and what the Summary writes today, after a separate export step the user must run. | None. This is a methods choice. It also settles whether the separate RUVIII export step stays. |
| D4 | A sample that one omic's QC removed. | The workbook removes it from every view. MOFA2's FAQ says missing values are ignored in the likelihood. | None; record the rule in the integration design either way. |
| D5 | The biological-sample key across omics, and technical replicates. | Every omic has `Run`, `group`, `replicates` and `tech_rep_group` (`<group>_<replicates>`). Nothing guarantees the replicate numbers mean the same sample in two omics. Trx, lipid and metab may carry technical replicates; the demo has none. | The Alignment tab proposes a mapping from `tech_rep_group` and the user confirms or edits it; the confirmed table is the integration design. How technical replicates enter the model waits for the workbook. |
| D6 | MOFA2's Python runtime. | `MOFA2` and `reticulate` are suggested packages. The workbook calls `run_mofa` twice, the first expected to fail. With basilisk the first run downloads an environment. KIT-130 and KIT-134 plan to bundle a Python with `mofapy2`. | An optional capability, probed; absent, the Model step refuses with a plain message and nothing else in the application changes; no in-R fallback. In the distributable it is bundled and `use_basilisk = FALSE` points at it. From a source install it is the user's to provide. |
| D7 | The minimum for a run. | "Two or more omics present" in the ticket. Two assays of one omic are two views. MOFA2's FAQ asks for more than 15 samples; the demo has 11. | Two or more views from two or more omics. Below 16 aligned samples: a visible warning, not a refusal. |
| D8 | Feature filtering before the model. | MOFA2 "strongly" recommends highly variable features per assay and views of similar size. The workbook filters nothing. | None; waits for the workbook. |
| D9 | Scaling. | The workbook Z-scales each feature and leaves `scale_views` off, sets `center_groups` off and fixes two factors. | Record as is; waits for the workbook. |
| D10 | The claim limit's wording and where it is shown. | Revision 1.3: factors are descriptive; no test of group separation from factors computed on a matrix corrected with the groups as replicate structure. MOFA2's pages neither state nor contradict it. | Shown on the Results tab beside any factor-against-group plot and printed in the report, with each view's RUV design from the manifests. |
| D11 | Display names and enrichment for any organism. | The workbook's lookups (a DIAMOND table, two UniProt downloads, a UniParc file, a GFF3) and the STRING species id are for one organism. An API key is written in the workbook. | Display names from each omic's own annotation (section 4.1); the stable id stays the model's feature id. The STRING key is typed by the user and never recorded; the key in the shipped workbook is removed when KIT-063 rewrites it. |
| D12 | The reference for exactness. | No trained reference exists in the repository; the shipped workbook reads a model file. | The verified workbook arrives with its inputs, its options and seed, and one trained model; KIT-062 compares numerically (the cross-machine rule). Without it, KIT-062 can be checked on the synthetic corpus only. |
| D13 | Splitting KIT-062. | With the GUI it holds the owner and model, the corpus, and three new tabs. | Two tickets: the owner, declarations and corpus; then the GUI tabs. |
| D14 | The DIA_LIMPA input in covariates mode. | The state is uncorrected; the corrected matrix is a product. | The input takes the product (section 4.3). Confirm. |
| D15 | Corrected matrices after KIT-126 and KIT-127. | Both move RUV into the model for the DA. | Both keep producing the corrected matrix as a product for integration and display. This is a requirement on those tickets, to be written into them by the planner. |

Carried over from the ticket and answered above: the Python capability (D6) and the minimum set of omics (D7).

## 3. Re-estimate

The stubs' estimates (6h, 10h, 8h; 24h in all) were written before the GUI form, before four omics were on the owner, and before the corpus was required. The science is still unknown, so each figure carries roughly plus or minus a third.

| Ticket | Was | Now | Why |
|---|---|---|---|
| KIT-061, the per-omic input | 6h | 14h | Seven routes instead of one. The contract's schema and its validator (2h). Four writers on the existing Summary hook: proteomics classic and DIA_LIMPA, the latter from a product (2.5h); lipid and metab shared (1.5h); transcriptomics, including the outcome of D3 (1.5h). The record entry and file hashes, and the session export's members (2h). Targeted tests per omic and the dependents (3h). Check, install, evidence (1.5h). More if D2 keeps today's files working beside the new ones. |
| KIT-062, the owner, the model and the GUI import | 10h | 32h | The corpus: one additive argument in four generators with byte-identical defaults, the composer, truth and both tiers (6h). The integration owner and the declarations for alignment, scaling and the model, with the model file and the three tables as saved products (7h). The Python probe and its absence path (2h). The module shell, the Import scan, the Alignment tab and the Model tab (10h). Module-CI and one e2e lane on the corpus (4h). The numerical comparison against the workbook and the evidence (3h). Proposed as two tickets of 15h and 17h (D13). |
| KIT-063, downstream, record, report | 8h | 24h | Results tab with the variance, factor, weight and data plots as saved assets, for any view and factor rather than the demo's (6h). Enrichment as declarations with the external resources recorded, and its tab; three network services and a mapping file (6h). The integration record, its replay and a report template (5h). The workbook rewritten over the module (3h). Tests and closure (4h). |
| Total | 24h | 70h | |

What would lower it: D2 answered as "replace" rather than "keep both"; D3 answered as "keep the existing RUVIII export"; enrichment left out of the first GUI version.
What would raise it: a workbook whose methods differ much from the starter's (feature selection, groups, several factors with per-factor enrichment); RUV becoming optional for transcriptomics.

If the workbook has not arrived by the 7 December checkpoint, KIT-061 and the corpus and owner half of KIT-062 can still be built and tested on synthetic data, because neither depends on the model's options. The Model tab's options, KIT-063's plots and the enrichment cannot be fixed without it.

## 4. Amendments proposed for the stubs

This run was told not to edit ticket files, so the three stubs are unchanged. The planner would carry these in when authoring them in full.

- **KIT-061.** Implement the contract of section 4 on the existing `summary$integration` hook. Name the matrix by role with the per-omic table of section 4.3. Depends on D1, D2, D3, D14, and on how KIT-126 and KIT-127 land (D15). The stub's `workflowIntegrationInput` and "a sample-alignment table emitted per omics" become the manifest and `samples.tsv`. Estimate 14h.
- **KIT-062.** Build the corpus of section 6 first. The stub's "samples dropped by any omics' QC dropped from all" and "Z-scales each feature within its view" are the starter workbook's behaviour, pending D4, D9 and the verified workbook. Add the Import and Alignment tabs (section 5.2). Exactness is numerical. Estimate 32h, or split (D13).
- **KIT-063.** The functions it places as declarations read global state today (rows DS4, E4, E6); "kernels unchanged" needs a decision on wrapping against repairing. Add the claim limit's display (D10) and the removal of the key from the shipped workbook (D11). Estimate 24h.

## 5. Acceptance against the ticket

| Criterion | Status |
|---|---|
| The register published before any repair, each row with a gap and a closing ticket | Done (section 3; 46 rows in six groups). Function and chunk names were read from the tree at `015d6577`. |
| The inputs contract written | Done (section 4). |
| KIT-061 to KIT-063 authored as planning-gated stubs | The stubs exist from 2026-10-26. Their amendments are proposed in section 4 above and not applied, by instruction. Remaining for the planner. |
| The user's decisions listed | Done (section 2 above; fifteen). |
| The four "corrected" matrices listed, one named per omic (revision 1.3) | Done: five found (section 2.2), one per omic named where the code allows, D3 open for transcriptomics. |
| The synthetic corpus specified (revision 1.2) | Done (section 6). |
| The MOFA2 source cited or the sentence dropped (revision 1.3) | Cited from MOFA2's FAQ; its limits stated (section 7.1). |
| No code changed | Two new documents under `docs/qa/workflow-execution/kit/kit-060/`. Nothing else touched. Not committed. |

## 6. Limits of this audit

- Read from code and the workbook only. No file in a real project was opened; the file names and layouts in section 2 are what the code writes, not what was observed on disk.
- The demo data was not opened, so the two differing sample ids (row A1) and the metabolomics feature column `metabolite` (row I4) are reported as written in the workbook, unresolved.
- The MOFA2 pages were read on 2026-11-17; the vignette is for MOFA2 1.22.1. The installed MOFA2's defaults (the seed among them) were not queried, because no R was run.
- Whether the owner's store under `scripts/<omic>/` outlives a session was not traced; the contract does not read it.
