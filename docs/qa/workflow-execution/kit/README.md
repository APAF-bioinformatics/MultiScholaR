# The workflow kit (KIT series)

Status 2026-10-09: series created; nothing built yet. Series authority: `tickets/kit/00-KIT-PLAN.md` (gitignored, local). Measurement protocol: [protocol.md](protocol.md). Diagram of the target architecture along the DIA workflow: `https://claude.ai/artifact/VZaUi4mPYvmqf2qSR2jsAd` (private artifact). Predecessor record: [perf-retirement-20261009.json](perf-retirement-20261009.json).

## What the kit is

The CORE series ([closeout](../core/closeout/outcome.md)) proved, on DIA proteomics only, that workflow state can live on disk as Parquet components and be hydrated per kernel; that an unchanged R kernel can be wrapped (identity, cache, projected hydrate, checkpoint) with bit-identical output; that per-sample partitioning with adaptive batches lowers peak memory where a step is per-sample independent; and that the three-process namespace floor was the largest fixed cost (CORE-018 cut it from 2.4 to 1.2 GiB). The kit turns those into functions that any workflow uses through declarations, so that a workflow is kernels plus declarations and the memory behaviour comes from the kit.

Two axes bind the functions:

- **The class declaration** binds storage. `workflowClassDeclaration(class)` says which slot is the assay (long or wide; a single table or a named list of tables), which is the design, which slots are tables, which columns are identity roles, and how to describe a state. The component checkpoint saves and reads any declared class. DIA's `PeptideQuantitativeData` and `ProteinQuantitativeData`, the lipidomics and metabolomics assay classes, and the result classes are declarations.
- **The independence unit** binds execution. There are few units, so there are few partitioners: `partitionByRows` (row-independent steps), `partitionBySamples` (complete-sample batches, the proven P14 engine), `partitionByFeatureGroup` (batches of whole feature groups), `partitionByPlex` (TMT), and `runProjected` (whole-matrix kernels on one projected object under the cap). A partitioner presents each batch in the shape the kernel already expects and appends the results; the kernel never learns it ran on a slice.

A **step declaration** is six facts: input class, fields and columns read, independence unit, kernel, output class, parameter names, plus a declared **global-statistic pre-pass** where the kernel takes a dataset-wide value (such as a quantile threshold): computed once over the whole input before any partition and recorded in the identity. A kernel that silently recomputes such a value per batch would change results, which is why the pre-pass is part of the contract and not left to the partitioner. The **operation runner** does the rest: identity (source references, operation, parameters, declaration version, software identity), cache lookup, projected hydrate, kernel call, checkpoint, release. The **operation registry** makes any declaration dispatchable through the managed worker (controller, monitor, worker in one slice under the admission cap) or in-process. The **on-disk state manager** is the default `state_manager` of every workflow module, so history is an immutable chain of references and modules hold references, not objects.

## Independence units

| Unit | Partitioner | Examples |
|---|---|---|
| row | `partitionByRows` | q-value filter, protein replicate filter (wide), log transforms |
| complete sample | `partitionBySamples` | min peptides per sample, P14 precursor-to-peptide rollup |
| feature group | `partitionByFeatureGroup` | min peptides per protein, peptide replicate filter, imputation, duplicate resolution, IQ (checked per batch); the peptide intensity filter after a declared global-statistic pre-pass |
| plex | `partitionByPlex` | TMT |
| whole matrix | `runProjected` | protein missing-value filters (dataset-wide quantile), between-sample normalisation, RUV, correlation, limma, Limpa, edgeR, enrichment |
| store | none | import adoption, design, reports, Save/Open |

A unit is declared only after the science owner confirms it, and the partitioner's exactness test checks it (per-batch results must hash-identical to per-item results). Whole-matrix steps do not reach constant memory; they reach one projected object under a cap with nothing else resident.

## Rules the kit keeps

- Kernels are never rewritten. A replacement kernel (such as the P14 SQL aggregate) needs an exactness proof against the R kernel.
- Memory decides; time is reported. Measurement works with the host as normally used.
- One mechanism per ticket; a valid unfavourable result is recorded, never retried.
- A new workflow writes one class declaration, one declaration per step, its kernels and its modules on the request/poll pattern. Nothing else. Lipidomics (KIT-031) is the acceptance test: if it needs new infrastructure, that is a kit gap.

## Stages

1. DIA proof (KIT-001–017): driver and step map, class declarations, runner, harness, partitioners, registry, production dispatch, on-disk state manager by default, reference-holding modules, app defects, step declarations, closure.
2. LFQ/TMT (KIT-020–021).
3. Lipidomics and metabolomics (KIT-030–032).
4. Transcriptomics (KIT-040), after TRX-005 extracts the notebook kernels.
5. Gates (KIT-050–052): the PERF-128/129/130 cases, carried verbatim.

Evidence per ticket lives in `kit-0NN/` beside this file; bindings in `bindings.json` (KIT-001).
