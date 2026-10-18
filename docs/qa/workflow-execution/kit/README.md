# The workflow kit (KIT series)

Status 2026-10-09: KIT-001 harness built (journey driver, coordinator, bindings, step map); the control build is `kit_001` (core_018 source plus the KIT journey launch record, the RUV-III technical-replicate structure and the proteomics quantile-ceiling removal, user decisions 2026-10-09). Series authority: `tickets/kit/00-KIT-PLAN.md` (gitignored, local). Measurement protocol: [protocol.md](protocol.md). Diagram of the target architecture along the DIA workflow: `https://claude.ai/artifact/VZaUi4mPYvmqf2qSR2jsAd` (private artifact). Predecessor record: [perf-retirement-20261009.json](perf-retirement-20261009.json).

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
| whole matrix | `runProjected` | protein missing-value filters (dataset-wide quantile), between-sample normalisation, RUV, correlation, limma, Limpa, edgeR, enrichment |
| store | none | import adoption, design, reports, Save/Open |

A unit is a property of the kernel's code and is determined by reading it, with the evidence recorded in the tracked step map (`kit-001/step-map.json`); only a step where the code is batch-safe but the intended science might not be, or where the code is ambiguous, is flagged with a question for the science owner. The partitioner's exactness test then checks every determination against the control build's output on every component the step writes (assay, design, support tables, ledgers, summaries, products). Partitioned declarations carry merge rules for the non-assay components (design recomputed once on the merged output, ledgers concatenated, counts summed). Whole-matrix steps do not reach constant memory; they reach one projected object under a cap with nothing else resident.

## Rules the kit keeps

- Kernels are never rewritten. A replacement kernel (such as the P14 SQL aggregate) needs an exactness proof against the R kernel.
- Memory decides; time is reported. Measurement works with the host as normally used.
- One mechanism per ticket; a valid unfavourable result is recorded, never retried.
- A new workflow writes one class declaration, one declaration per step, its kernels and its modules on the request/poll pattern. Nothing else. Lipidomics (KIT-031) is the acceptance test: if it needs new infrastructure, that is a kit gap.

## Kernel contracts the kit adapts to (KIT-022, for the KIT-017 checklist)

Kernels are never rewritten, and three of their idioms are therefore adapted to, not removed. They are named here so the closure checklist reads them as adapters, not as kit design:

- **Legacy globals at publication: none (KIT-023).** The figure writer (`updateProteinFiltering`) takes `project_dirs` as an argument and falls back to the global only when none is passed; nothing in its chain reads `config_list`. The kit always passes the directories, so the publication phases bind nothing into the global environment (the controller is the application process, shared by every session), and the test at the figure boundary fails if a binding returns. The kernels' own global fallbacks stay in the kernels and are never reached from the kit.
- **Call-stack parameter resolution.** The kernels resolve their parameters with `checkParamsObjectFunctionSimplify` (a `dynGet` over the calling frames, falling back to the object's `args` under the calling function's name). The kit always calls a kernel with explicit arguments and through the kernel's own name (declaration bodies call the legacy `…Fn` aliases so the recorded args keep their name; `uniprot_tbl` must be a table, an empty one when there is none).
- **The active owner.** Kit-level readers (`readDIAContrastOutput`, the DA display and product functions) find the store through `getDIAActiveManager()`, which resolves the innermost open `withDIAStateScope` owner first and the Shiny session's owner second (KIT-022). Kit code, workers and scripts open the scope; modules rely on the session.

## Before a counted member exercises a new operation role or class

KIT-010 needed three launches because its pre-execution checks covered the compute, the controller's completion
and publication and the verifier, but not the adapters vehicle's input check nor the worker entry point's dispatch:
plain code paths that only the installed member reached. Before a counted member runs an operation with a new
compute role, input class or variant, exercise the whole chain in process on the fixture, and keep the exercise as
a test where one fits:

1. the registry entry validates and names handlers that exist (`test-workflow-operation-registry.R`);
2. the controller's compute helper reaches the handler (`computeWorkflowPeptideAction` by compute role);
3. the worker entry point dispatches the role from a frozen request (`workflowWorkerMain`; the KIT-010 test in
   `test-workflow-partitions.R` is the template);
4. the adapters vehicle's checks hold for the population kind (the input component is the class declaration's
   assay; the variant is in the registry entry; the input state is the population's);
5. the completion and publication phases run on the fixture;
6. the verifier's branch compares the member's state with its oracle (the in-process check on the real population).
7. (KIT-017 post-hoc audit PH-03) before a pair rule is bound to a memory mechanism, the mechanism's slice cost at the
   step it changes is measured on the screen corpus and the owner of that step's peak is named from the trace (a
   child process such as pandoc shows as a ramp-and-drop; the worker's baseline as a plateau step). A rule bound to
   a resident-memory claim is still judged on the slice's step peaks.

## Stages

Stage 1 (DIA) closed 2026-10-15 at `ced86a66` (KIT-017, [outcome](kit-017/outcome.md), [checklist](kit-017/kit-checklist.md)): the pair's required resident-end reduction met, the per-step guard failed on the DA step. The [post-hoc audit](kit-017/KIT-017-POSTHOC-AUDIT-20261016.md) named the cause (pandoc's self-contained pass over the Glimma page); KIT-036 ([outcome](kit-036/outcome.md)) removed it and was accepted on 2026-10-16 (DA step 0.675, resident end 0.921), so stage 1 ends on an accepted build and `kit_036` is the control.

1. DIA proof (KIT-001–017): driver and step map, class declarations, runner, harness, partitioners, registry, production dispatch, on-disk state manager by default, reference-holding modules, app defects, step declarations, closure.
2. Stage 2, re-cut on 2026-10-19 ([review](REVIEW-20261019-stage2-recut.md)), in order: the record states what ran and reads for a reviewer (KIT-038, KIT-039; the 2026-10-19 audit's IM findings); LFQ/TMT imports on the owner, their sessions, reports and one measured LFQ pair plus a TMT exactness member (KIT-020, KIT-042, KIT-043, KIT-021); Limpa as three whole-matrix declarations and `fastlimpadpc` as a declared, recorded backend (KIT-046, KIT-047; the acceptance rule for a tolerance-aligned kernel is a user decision at KIT-047's closure); the owner on every platform by capability, not platform identity (KIT-044); the in-memory proteomics path retired (KIT-048); omic-agnostic `workflow*` names for the owner, scope, commit, component and display primitives with the DIA names kept as aliases (KIT-045). `partitionByPlex` is dropped: no kernel is per-plex.
3. Lipidomics and metabolomics (KIT-030–032), on the single-path proteomics pattern and the KIT-045 names.
4. Transcriptomics (KIT-040), after TRX-005 extracts the notebook kernels.
5. Gates (KIT-050–052): the PERF-128/129/130 cases, carried verbatim.

Evidence per ticket lives in `kit-0NN/` beside this file; bindings in `bindings.json` (KIT-001).

## Checks that need the proof environment

- `R CMD check --no-manual --no-tests --no-build-vignettes` only configures with an installed dependency set on `R_LIBS` (the configure step locates the exact-sum artifacts through installed packages): set `R_LIBS=<proof>/core/browser-smoke.<control>/library`. Two suggested packages (GSEABase, MOFA2) are not installed on the measurement host, so set `_R_CHECK_FORCE_SUGGESTS_=false`. Build the tarball from a clean worktree under `~/.cache` (a working-tree build carries ignored local artifacts and adds NOTEs) and compare against a baseline worktree built the same way. Baseline after KIT-002 remediation: 3 warnings, 3 notes.
- Installed checks (`install-corrected`, `writer-r-preflight`) run under `<proof>/core`, where the exact-sum artifacts live; the kit proof root holds only runtime stages and corpora.
