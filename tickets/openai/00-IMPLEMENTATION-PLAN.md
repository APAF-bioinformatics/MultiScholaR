# Peptide QC Scientific Integrity Implementation Plan

## Objective

Make the DIA peptide-QC pipeline scientifically explicit, reproducible, and
auditable from imported precursor rows through protein quantification. Preserve
the corrected identification rule from commit `a2ecddc`:

- active protein identity: `Protein.Group`;
- identification gate: at least one distinct `Stripped.Sequence` and at least
  two distinct `Modified.Sequence` values per protein group;
- evidence population: q-valid, proteotypic identifications across the whole
  experiment, independent of sample recurrence;
- evidence remains frozen after the identification stage while quantitative
  survivor counts remain separately dynamic.

## Locked Scientific and Product Decisions

1. **Q thresholds are fixed at 0.01.** The GUI must not imply that a different
   value was applied. `Q.Value`, `Global.Q.Value`, and the global protein-group
   metric `Global.PG.Q.Value` are distinct fields and must never be silently
   substituted for one another. Legacy/pre-filtered inputs require an explicit
   provenance mode rather than fabricated confidence.
2. **`Protein.Group` is the DIA-NN quantitative identity.** `Protein.Ids` is
   accession provenance/annotation and must not replace, split, or merge the
   active group key downstream.
3. **Modified and unmodified forms are merged for quantitative protein input.**
   Charge states and peptidoforms may be summed into one stripped-peptide
   quantity. Identification peptidoform identities remain preserved in frozen
   evidence fields and the audit ledger.
4. **Technical-replicate mean imputation is retained.** It is only valid within
   an explicitly declared technical-replicate group and must count distinct run
   IDs. It must remain flagged with `is_imputed` and must not masquerade as an
   observed value in audit summaries.
5. **The IQ compatibility behavior is retained except for the protein key.**
   Synthetic compatibility q-values, conversion of remaining missing IQ input
   values to zero, and `normalization = "none"` are intentional. The IQ
   `primary_id` must change from hardcoded `Protein.Ids` to the active
   `Protein.Group`, with the compatibility transformation recorded.
6. **Raw imported rows are immutable evidence.** Filters produce analysis views
   and removal ledgers; they do not destroy the import snapshot.
7. **No hidden defaults.** GUI values, config values, S4 arguments, executed
   values, and saved audit values must agree. Both the historical `*_cutoff`
   aliases and canonical `num_*_thresh` names remain accepted, with one resolved
   value recorded.

## Audit Finding Coverage

| Audit finding | Planned correction | Ticket |
|---:|---|---|
| 1 | Apply the fixed `0.01` gate to correctly named DIA-NN confidence fields; make legacy confidence provenance explicit. | [PQC-001](PQC-001.md) |
| 2 | Replace percentage-derived peptide missingness decisions with direct observed-run and passing-group counts. | [PQC-002](PQC-002.md) |
| 3 | Count distinct `(Protein.Group, Stripped.Sequence)` identities per sample and align the default threshold across config/UI/methods. | [PQC-003](PQC-003.md) |
| 4 | Count distinct run support per replicate group, while preserving the approved “supported in any group” policy explicitly. | [PQC-003](PQC-003.md) |
| 5 | Correct technical-replicate eligibility, the 50% boundary, mean calculation, exclusion metadata, and imputation provenance. | [PQC-004](PQC-004.md) |
| 6 | Validate precursor uniqueness/scale/missingness and intentionally sum modified, unmodified, and charge forms into stripped-peptide quantities. | [PQC-005](PQC-005.md) |
| 7 | Treat DIA-NN decoys and contaminant targets separately, default-exclude confident classifications, and preserve an exact ledger. | [PQC-006](PQC-006.md) |
| 8 | Preserve the approved IQ compatibility transformations but carry active `Protein.Group` through every downstream seam. | [PQC-007](PQC-007.md) |
| 9 | Add immutable parent-linked state records, deterministic summaries/digests, identity-level removal ledgers, and a headless E2E audit. | [PQC-008](PQC-008.md) |

The additional silent threshold-plumbing regression is covered by PQC-001:
config-only `peptides_per_protein_cutoff` and
`peptidoforms_per_protein_cutoff` values must exercise the historical alias
path and resolve to the same executed canonical `num_*_thresh` values as the
GUI and direct-call paths.

## DIA-NN Decoy and Contaminant Policy

DIA-NN generates decoys for FDR estimation. Its main documentation states that
decoy identifications are exported for downstream use only when
`--report-decoys` is requested, in which case the `Decoy` column identifies
them. Therefore decoy filtering is normally a no-op for a standard main report,
but MultiScholaR must exclude `Decoy == 1` whenever that signal exists.

Common-contaminant/cRAP proteins are different: when their sequences are in the
FASTA/library, they are target sequences. DIA-NN exposes
`--cont-quant-exclude [tag]` for explicitly tagged proteins, but it cannot infer
an arbitrary downstream cRAP policy. MultiScholaR will therefore:

- classify explicit `Decoy`, contaminant columns, and recognised decoy or
  contaminant ID tags;
- support a versioned packaged cRAP accession manifest plus an optional user
  manifest;
- match accessions token-wise, never by an unbounded description substring;
- default-exclude confidently classified decoys/contaminants from biological
  identification and protein counts;
- retain every excluded row and the exact classification reason in the ledger;
- warn, rather than guess, when no reliable contaminant signal is available.

Primary reference: <https://github.com/vdemichev/DiaNN>. Relevant documented
behaviour includes `--report-decoys`, the `Decoy` output field,
`--cont-quant-exclude`, and the recommended main-report q-value filters.

## Pipeline Contract

```text
immutable DIA-NN import
  -> fixed 0.01 precursor/global/global-PG q gate + Proteotypic == 1
  -> decoy/contaminant classification and biological-analysis exclusion
  -> frozen distinct identification evidence (Protein.Group)
  -> >=1 stripped peptide AND >=2 modified peptidoforms
  -> precursor rollup (linear sum; forms intentionally merged)
  -> count-based intensity/missingness filtering
  -> distinct peptide identities per sample
  -> distinct runs per technical-replicate group
  -> technical-replicate mean imputation with flags
  -> IQ handoff keyed by Protein.Group
  -> ProteinQuantitativeData retaining Protein.Group plus Protein.Ids provenance
```

Each arrow must produce a parent-state reference, resolved parameters,
before/after summaries, removed-identity ledger, and deterministic data digest.

## Shared Invariants

- One design row per run; duplicate or unmatched run mappings fail with the
  offending IDs.
- Peptide quantitative identity is
  `(Run, Protein.Group, Stripped.Sequence)` after rollup.
- Precursor quantitative identity is
  `(Run, Protein.Group, Precursor.Id)` before rollup when `Precursor.Id` exists.
- Counts use `n_distinct()` over declared identity fields, never incidental row
  counts.
- Quantities must be finite, non-negative, and linear before summation.
- A partially observed sum ignores missing components; an all-missing group
  remains `NA`, not zero.
- Boundary language is literal: “at least” uses `>=`; “maximum missing” uses
  `<=`.
- Frozen identification counts and current quantitative-survivor counts have
  different names and are never presented as the same measure.

## File-by-File Implementation Map

### PQC-001 — q/FDR and parameter contracts

- `R/func_prot_qc_peptide_replicate_filters.R`: explicit fixed-metric q gate;
  remove cross-metric substitution paths; preserve exact evidence annotation.
- `R/func_prot_qc_peptide_methods.R`: resolve schema/provenance mode and reject
  silent confidence fallbacks.
- `R/mod_prot_qc_peptide_qvalue.R`: read-only 0.01 presentation and exact metric
  labels.
- `inst/config/config.ini`: one canonical 0.01 contract.
- Existing protein-evidence method/module tests: add the previously missed
  alias-only config regression and conflicting GUI/config cases.

### PQC-006 — decoys and contaminants

- New `R/func_prot_qc_exclusion_helpers.R`: tokenisation, classification,
  manifest provenance, exclusion reasons, and analysis/raw split.
- `inst/extdata/contaminants/`: versioned cRAP accession manifest and source
  metadata, subject to source/licence verification during implementation.
- Q-value module/method: invoke classification after confidence filtering and
  before frozen evidence calculation.
- Import and q-filter tests: explicit `Decoy`, `CON__`, cRAP accession, mixed
  protein groups, user manifest, and unclassified cases.

### PQC-005 — rollup integrity

- `R/func_prot_rollup_methods.R`: assert precursor uniqueness and linear scale;
  use all-missing-preserving sums; intentionally merge charge/modification
  forms; preserve group/accession/evidence provenance.
- `R/mod_prot_qc_peptide_rollup.R`: report aggregation and missingness effects.
- Rollup tests and documentation: duplicate, partial-NA, all-NA, modified plus
  unmodified, multiple charge, and `Protein.Group != Protein.Ids` cases.

### PQC-002 — intensity and missingness

- `R/func_prot_qc_peptide_group_filters.R`: count distinct observed runs and
  evaluate per-group minimum replicate counts directly.
- `R/func_prot_qc_peptide_methods.R`: prefer the object's normalised peptide
  quantity for the cross-sample percentile filter while preserving an explicit
  quantity-column override.
- `R/mod_prot_qc_peptide_intensity.R`: pass `min_reps_per_group` and
  `min_groups` directly; retain old percentage parameters as deprecated API
  aliases only.
- Group-aware tests: unequal groups, duplicated input, absent group/peptide
  combinations, exact boundaries, strict mode, and non-finite values.

### PQC-003 — sample and replicate identity counts

- `R/func_prot_qc_peptide_replicate_filters.R`: count distinct
  `(Protein.Group, Stripped.Sequence)` per run and distinct runs per
  technical-replicate group/feature.
- `R/func_prot_qc_peptide_methods.R`: resolve one canonical sample threshold.
- Sample and replicate modules: initialise from the resolved config, display the
  exact executed count definition, and record support by group.
- `inst/config/config.ini`: align the packaged sample minimum to 500, matching
  the existing GUI intention; user/config overrides remain authoritative.
- Preserve the current “supported in any group” feature-retention policy, but
  make it explicit and auditable.

### PQC-004 — technical-replicate imputation

- `R/func_peptide_qc_imputation.R`: count distinct run IDs; calculate means over
  all observed technical-replicate measurements; use `<=` for maximum missing;
  remove the unconditional print and default HEK regex exclusion.
- `R/mod_prot_qc_peptide_impute.R`: expose/configure an explicit metadata-based
  exclusion policy and report observed/imputed counts.
- Preserve raw-quantity mean imputation for the intended IQ input and preserve
  `is_imputed` on every row.
- Tests: equal-valued replicates, duplicate runs, 50% boundary, all-missing
  groups, explicit QC exclusion, and no technical replicates.

### PQC-007 — `Protein.Group` lineage

- `R/func_pept_s4_core.R`: do not relabel the active group identity as
  `Protein.Ids` while building peptide matrices; retain a reversible feature-key
  map.
- Peptide normalisation/missingness/accession methods: preserve the active
  protein-key metadata through matrix round trips.
- `R/mod_prot_qc_protein_rollup.R`: use `peptideS4@protein_id_column` as IQ
  `primary_id`, parse the resulting column dynamically, and construct the
  protein object with `Protein.Group` as its active key.
- Preserve the approved IQ compatibility q-values, zero conversion, and
  `normalization = "none"` unchanged, while recording them in audit metadata.
- Update DIA workbooks at the peptide-to-protein seam and add tests where
  `Protein.Group` and `Protein.Ids` deliberately differ.

### PQC-008 — complete audit trail and E2E closeout

- New `R/func_prot_qc_peptide_audit_helpers.R`: deterministic summaries,
  identity ledgers, parameter capture, parent-state references, and hashes.
- `R/utils_workflow_state.R`: accept optional audit metadata without breaking
  existing callers.
- All peptide-QC modules: emit a common audit record and prerequisite/data-shape
  contract at each applied step.
- New `tools/audit/run_peptide_qc_audit.R`: headless reproduction against an
  arbitrary DIA-NN TSV without embedding private data paths.
- New focused tests plus module-CI/E2E assertions: exact state order, no hidden
  parameter differences, deterministic ledgers, and rollback/revert behaviour.
- Real-data closeout, when the two local files are available, must reproduce
  their input hashes and the accepted post-identification counts before the
  ticket can be completed.

## Verification Strategy

Each ticket carries an exact focused command. Closeout additionally runs:

```bash
Rscript --vanilla -e 'devtools::test(filter = "prot-(02|03|04).*peptide|prot-qc-peptide|module-ci-prot-peptide")'
Rscript --vanilla tools/ci/run-module-ci.R --modules proteomics-peptide-qc
```

If supported by the local E2E harness:

```bash
Rscript --vanilla tools/ci/run-e2e-ci.R --scenario proteomics-dia
```

Real-data verification is opt-in and never commits source data:

```bash
Rscript --vanilla tools/audit/run_peptide_qc_audit.R \
  --input /home/doktersmol/Downloads/MS_tester/proteomics/cotton_report.tsv \
  --output /tmp/multischolar-pqc-cotton
Rscript --vanilla tools/audit/run_peptide_qc_audit.R \
  --input /home/doktersmol/Downloads/MS_tester/proteomics/KV_DIANN_report.tsv \
  --output /tmp/multischolar-pqc-kv
```

Expected three-gate baseline evidence before later quantitative QC changes:

- cotton SHA-256 `b3bd682e1e34e7bb494be8d4162d25abb9cd6162f8f149818f6fd0e5ed0fba0a`:
  18,662 rows and 1,482 protein groups after all three mandatory 0.01
  gates, with 519 passing `>=1/2`;
- KV SHA-256 `d18f74c8a98e8ec9ec4ac0e9a4b32b6a1a0719e2b9ed46125e9e3f2aa267698e`:
  159,658 rows and 2,160 protein groups after all three mandatory 0.01
  gates, with 1,763 passing `>=1/2`.

The superseded 19,551/1,673 and 160,121/2,311 figures are retained as an
explicit two-gate (`Q.Value` plus `Global.Q.Value`) diagnostic. They omit the
required `Global.PG.Q.Value <= 0.01` gate and therefore cannot be the accepted
post-confidence baseline. The final `>=1/2` result is unchanged by correcting
that stale label (519 cotton; 1,763 KV).

Contaminant exclusion will intentionally create a separately labelled
biological-analysis count; it must not overwrite the all-target evidence count.

## Compatibility and Rollback

- Existing public helper arguments remain accepted; deprecated aliases produce
  one warning and resolve to the same recorded canonical value.
- Non-DIA workflows must remain unchanged. Every shared-helper change requires
  LFQ/TMT-focused regression coverage or must be isolated to a DIA helper.
- Existing state objects without audit metadata remain readable.
- Every ticket is independently revertible. Schema/data migrations are additive
  until the integrated E2E ticket proves round-trip compatibility.
- Ticket status remains `pending` until its simulation, acceptance, regression,
  backward-compatibility, and rollback checklists have concrete evidence.
