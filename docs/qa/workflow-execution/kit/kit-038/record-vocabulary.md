# The analysis record's provenance vocabulary (KIT-038, schema `multischolar.analysis_record.v1`)

Each parameter row of a step, of a `preceded_by` entry and of an analysis carries one of these values in `provenance`. The set is closed: `validateWorkflowAnalysisRecord` rejects any other string. The legend rendered in the supplementary document (KIT-039) is this table.

| Value | Meaning |
|---|---|
| `default` | The value equals the shipped `config.ini` default for the parameter's mapped key (exact numeric comparison; a comma-separated shipped list is compared element-wise only for a multi-valued recorded value). |
| `user_set` | The parameter has a shipped default and the recorded value differs from it. |
| `not_set` | The parameter has a shipped default but the recorded (requested) value was `NULL` or `NA`; the module's own fallback applied (for example an automatic RUV `k` chosen by the optimisation, whose realised value is in the step's `outcomes`). |
| `no_shipped_default` | `config.ini` ships no key for this parameter; the module's own default applied whether or not the user changed it. KIT-027 wrote this as `no_default`; records carrying that spelling validate (read as `no_shipped_default`) and are never written by builds from KIT-038 on. |
| `fixed` | A kernel constant the operation always passes (declared `fixed` in the registry; for example the P14 rollup's `backend` and the IQ rollup's `filter_double_less` and `normalization`). |
| `input` | A value taken from the data or an earlier step (a list, a table, a hash), never a setting. |
| `unknown` | No registry mapping exists for the name: a `legacy` step's arguments listed whole, or a declared name with no recorded value (`source = missing`). |

`source` names where the value was read: `args` (the kernel's `slots$args` section), `provenance` (the values the operation received, recorded on the head), `config` (the state's commit config), `manifest` (the DA result manifest), `fixed` (the head's fixed constants), `registry` (a fixed constant the registry declares but the head entry predates, for example the P14 backend on stores written before KIT-038), or `missing`.

`outcomes` (every step): `rows`, `columns`; plus whatever the commit recorded under `config$outcomes` or `config$skipped`; for the normalisation step `ruv_applied` (from the state name), `k`, `negative_control_percentage` (the realised values from the commit config) and `negative_control_count` (`sum(controls)` of the saved optimisation product; `NA` when RUV was skipped).

`software` (record level): `package`, `version`, `sha`, `build` (the installed `MultiScholaR.rdb` sha256, the key of the head's `software_environments`), `r_version`, `platform`, `blas`, `lapack`, `omp_threads`, `blas_threads`, `duckdb_threads` (`rollup`, `codec_read`, `codec_bulk`), `os_release`, `locale`, `packages`; each environment field reads `"not recorded"` when the build that wrote the last step predates KIT-038. Each step's `software` carries `version`, `sha` and `build`.
