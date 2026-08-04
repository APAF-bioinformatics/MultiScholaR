# Peptide QC Scientific Integrity Board

Implementation is intentionally not started. All tickets are `pending` and all
checklists remain unchecked until supported by implementation or test evidence.

| Order | Ticket | Scope | Depends on |
|---:|---|---|---|
| 1 | [PQC-001](PQC-001.md) | DIA-NN q/FDR semantics and threshold plumbing | — |
| 2 | [PQC-006](PQC-006.md) | Decoy/contaminant classification and exclusion | PQC-001 |
| 3 | [PQC-005](PQC-005.md) | Precursor-to-peptide rollup integrity | PQC-001, PQC-006 |
| 4 | [PQC-002](PQC-002.md) | Intensity/missingness direct-count rules | PQC-005 |
| 5 | [PQC-003](PQC-003.md) | Sample and replicate distinct-identity counts | PQC-002 |
| 6 | [PQC-004](PQC-004.md) | Technical-replicate imputation | PQC-003 |
| 7 | [PQC-007](PQC-007.md) | `Protein.Group` lineage into protein quantification | PQC-004, PQC-005 |
| 8 | [PQC-008](PQC-008.md) | Audit ledger, stage contracts, and E2E evidence | all prior tickets |

The controlling design is [00-IMPLEMENTATION-PLAN.md](00-IMPLEMENTATION-PLAN.md).
