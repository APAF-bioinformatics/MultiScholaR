# KIT-015 pre-execution audit remediation (2026-10-14, user: "fix findings all")

Register: `KIT-015-PREEXEC-AUDIT-20261014.md` (unchanged). Package code changed (PE-01, PE-04), so the tree is reinstalled and rebound before launch; the harness and generator changed (PE-03), so the bindings are regenerated.

| ID | Closure | Evidence |
|---|---|---|
| PE-01 | Repaired. The IQ declaration's call passes `pdf_out = NULL`: IQ writes no QC plot into the process working directory (the module's in-memory call keeps its default, untouched). | `test-workflow-prot-declarations.R` publication test asserts no `qc-plots.pdf` in the working directory after the declaration runs; IQ step rerun in process against `journey_iq` after the change: exact, 198 s (`in-process-protein-qc-check-iq.json`). Reinstall `kit_015` = `core/browser-smoke.hXymGJ` from `3ad41274`; stage `kit-015` validates (4 cases). |
| PE-02 | Repaired. The four-state publication loop (IQ files, accession results RDS, figures for every state, the replicate's table) is a kept test. | `test-workflow-prot-declarations.R`. |
| PE-03 | Repaired. The generator records the bound FASTA table's sha256 beside its path (`inputs$fasta_table`); `kitCases` refuses a single member whose bound input is missing or changed, so `validate` and `run` stop before any member. | Harness test (a changed file refuses with "missing or changed"); bindings regenerated. |
| PE-04 | Repaired. The replicate compute carries its product field and display text like the KIT-015 values; `completeDIAProteinAction` has one shape and refuses a value without a field. | `test-workflow-prot-declarations.R`; KIT-010 completion and publication tests unchanged and green. |
| PE-05 | Claim limit, recorded in the outcome; no change. | — |
| PE-06 | Pre-existing baseline failures, recorded; no change. | — |
