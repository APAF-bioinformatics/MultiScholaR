# KIT-216 R2d: the re-frozen files after the re-freeze (predictions committed before the run)

The run compares:
- **Head:** the kit-216 tip, with the re-freeze 28a42a75 on staging 80f26158.
- **Base:** `kit-216-r2c-base` 923cc370, which is tests-only plus staging 80f26158.

Files: the three re-frozen files, their two readers, KIT-213's new file and the KIT-216 file. Expected cost: about 80
job-minutes.

| File | Verdict | Where |
|---|---|---|
| `test-omic-enrich-module.R` | `message-changed` | only W1 `kernel_clusterprofiler` run `states` (not re-frozen: base fails on 1167a2b0 vs d8f57c99, head on "<unchanged>"); everything else as on base |
| `test-omics-enrich-core.R` | `message-changed` | only the `kernel_clusterprofiler` reader at run `states`; everything else as on base |
| `test-omic-summary-module.R`, `test-omics-summary-core.R` | `same` | — |
| `test-gui-entry-proteomics.R` | `same` | `prot_dia` fails on both sides at the pre-existing artifact lines; the state digests other than `enrichment_completed` were left as stored |
| `test-kit213-study-parameters-labels.R` | `same` | — |
| `test-kit216-prot-session-after-enrichment.R` | `fixed` | — |

Any other difference: stop and report.
