# KIT-112 W1: does the DA click hydrate the protein matrix before the kernel, and would D3's request half? (2026-11-10)

The lead's question from W0 (`memory/baseline-large.md`). On `111q`, the controller rises +0.31 GiB within about 2 s of the
DA click and +0.17 GiB about 28 s in. Which statements in the click path load the protein matrix or the design in the
controller? Would D3's request half (the gate, the validation and the DA-arguments commit) load anything?

Method:
- Read the code at the start commit `64e7c2ae` (R/ = `c568dfa5`).
- Measure at fixture scale with [request-half-hydration.R](request-half-hydration.R). It builds a synthetic DIA state,
  accepts it at `correlation_filtered`, and runs each statement of `runDIAReferenceDaAnalysis` before the contrast loop
  on its own, then the kernel's projected read.
- Each statement runs twice, cold and warm. For each run the script records the R heap's peak over the statement
  (`gc(reset = TRUE)` before, then `gc()`'s "max used" after; this is an upper bound that includes garbage not yet
  collected) and the size of the value the statement returns.
- No large run.

## Answer

**By code, before the kernel the click path reads metadata and the design only. The protein matrix is first loaded
inside the compute half.** D3's request half needs nothing it does not read today:
- the state records (bounded RDS records of at most 8 MiB, `utils_workflow_dia_state.R:249`);
- the design table;
- the zero-row schema of the other tables.

None of these grows with the number of proteins. The +0.31 GiB at the click and the +0.17 GiB at about 28 s are the
kernel's: the projected read of the protein table, the matrix copies, the fits and the per-sample result tables. All of
that moves to the worker under D3. So W0's "better" target is not ruled out by the request half.

## The click path, statement by statement (start commit)

The observer `mod_prot_da_handlers_run_analysis.R:25-76` calls `runDIAReferenceDaAnalysis` (`:47-48`), which calls the
statements below (`utils_workflow_dia_da.R:341-433`). The last column gives the measured heap peak, warm, at 4,000 and
16,000 proteins × 300 samples.

| Line | Statement | What it reads | Data-scaled? | Peak, 4k → 16k proteins |
|---|---|---|---|---|
| :345 | `assertDIAExportReady(manager, allow_action = "da")` | the head snapshot's `action_outputs` / `completed_outputs`; a journal JSON only if outputs are pending; the QC-pending check runs only on QC states (`utils_workflow_saved_exports.R:280-308`) | no | +9.4 → +9.4 MiB |
| :353 | `readDIAStateMetadata(manager, source)$data` | the source's state record (component references, slots, dimensions) | no | +10.2 → +10.2 MiB |
| :354 | `readDIAExecutionComponents(manager, source, fields = "design_matrix")` | the design table; every other table becomes a zero-row schema (`rebuildDIACheckpointObject`, `utils_workflow_dia_state.R:436-490`: `prototype()` opens a metadata-only read scope) | no (the design is one row per sample) | **+152 → +153 MiB**, value 0.10 MiB |
| :355 | `validateProtDaModelAndContrasts(prototype, …)` | the design and the protein table's column names, from the zero-row schema (`utils_workflow_dia_da_validation.R:112-186`) | no | +17.3 → +17.3 MiB |
| :366 | `workflow$da_ui_params <- …` | nothing | no | — |
| :367 | `readDIAStateMetadata(manager)` | the head's record | no | +0.3 → +0.3 MiB |
| :383 | `commitDIAReferenceState(…)` | validates the component references without reading payloads (`payload = FALSE`), saves a record of at most 8 MiB and the head (`utils_workflow_dia_state.R:848-895`) | no | +14.1 → +14.1 MiB (cold +45) |
| `utils_workflow_operations.R:151` | `workflowComponentScientificIdentity` (the contrast product's identity) | the source's record; component entity ids only | no | +6.9 → +6.9 MiB |
| `utils_workflow_operations.R:174` | `runProjected` → `readWorkflowExecutionComponents(protein_quant_table, design_matrix, protein_id_table)` | **the protein table**, then the kernel: `as.matrix(column_to_rownames(...))` (`func_prot_da_model_methods.R:348`), the base-fit identity's digest of the matrix, `duplicateCorrelation` and `lmFit`, the contrasted fit, and the result tables with per-sample columns (`norm_counts`, `da_proteins_wide`; another `as.matrix` at `:565`) | **yes** | **+74 → +145 MiB** for the read alone (value 10 → 40 MiB) |

The same measurement at 4,000 proteins with 100, 300 and 600 samples:
- the design read: +149, +152 and +157 MiB;
- the commit: +14.1 MiB at every size;
- the projected read: +37, +74 and +131 MiB.

The raw output is in [request-half-hydration.txt](request-half-hydration.txt).

## Reading the +0.31 and +0.17 GiB on `111q`

These are inferences from the code and the timings, not measured on the corpus.
- **The rise at the click (+0.31 within about 2 s, then flat).** The request half is metadata-sized and takes well under
  a second. The first statement that grows with the data is `runProjected`'s read. It is followed by the kernel's matrix
  copy and the base-fit identity's digest of that matrix. All three are held while `duplicateCorrelation` runs on 300
  samples. That explains a level that holds from about 2 s to about 28 s.
- **The rise at about 28 s (+0.17).** It starts when the fit finishes. Then come `lmFit`, `contrasts.fit`, `eBayes`,
  `topTable`, the full-population q-values, and the per-sample result tables (`norm_counts`, `da_proteins_wide`, a
  second `as.matrix`), followed by the product's serialisation.

Both are in the compute half (D3: `computeDIAReferenceDaContrast`, run on the worker). None is in the request half or
the accept half.

## What D3's request half would load

D3's request half is `utils_workflow_dia_da.R:345-389`, unchanged. It runs the statements in rows :345 to :383 above,
on both paths. It needs:
- the design (for the validation);
- the protein table's column names (for the validation's sample check), which come from the zero-row schema;
- the source's and the head's records (for the manifest skeleton, the arguments comparison and the commit).

It needs no assay values. The accept half (`:407-423`) reads references and the base-fit reference's factors only.

**One cost to keep in view, though it is not the matrix.** The design-only read (`:354`) peaks at about 0.15 GiB of R heap
at every size probed (100 to 600 samples, 1,000 to 16,000 proteins), against a returned value of 0.1 MiB. By reading,
the cost is in the read scopes that build the zero-row schemas of the protein and identifier tables. It is transient and
does not grow with the data in the range probed. On the 300-sample corpus, the controller's heap has already been far
higher earlier in the session (1.7 to 2.9 GiB). So a transient of this size is likely absorbed by heap the controller
already holds, and is unlikely to set a new anonymous-memory peak. That is not proven at corpus scale. If `112q`'s DA
peak lands near the "better" threshold, this is the first place to look, and W2 should not add another prototype read
to the request half.

## Hand-off to W2

Keep the request half to the statements above, with no added read. The managed path then holds, in the controller:
- the request half's metadata-sized reads;
- the transient design read;
- the accept half's record writes.

The matrix stays on the worker.
