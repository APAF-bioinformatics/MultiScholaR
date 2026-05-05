# E2E Test Fixtures

Synthetic, deterministic fixture data for the MultiScholaR end-to-end test suite.
All data is generated from `set.seed(20260428)` — no real experimental data is used.

---

## Fixture Architecture

Fixtures are **manifest-driven**: `manifest.json` is the single source of truth for
all 8 workflow lanes. Each lane defines its format, sample layout, and expected outputs.

```
tests/testdata/e2e/
├── manifest.json             # Lane catalogue (schema_version 1.0.0)
├── cross_omic_packs.json     # Multi-omic integration pack definitions
├── generate_fixtures.R       # Regeneration script (developer tool)
├── report_templates/         # Stub Rmd templates for report rendering tests
│   ├── DIANN_report.rmd
│   ├── TMT_report.rmd
│   ├── LFQ_report.rmd
│   ├── metabolomics_report.rmd
│   └── lipidomics_report.rmd
├── prot_dia/                 # DIA-NN proteomics (long format, no LIMPA)
├── prot_dia_limpa/           # DIA-NN proteomics (long format, with LIMPA)
├── prot_tmt/                 # TMT proteomics (Proteome Discoverer format)
├── prot_lfq/                 # LFQ proteomics (MaxQuant proteinGroups.txt)
├── metab_lc/                 # LC-MS metabolomics (MS-DIAL format, Pos+Neg)
├── metab_gc/                 # GC-MS metabolomics (MS-DIAL format, GCMS)
├── metab_combined/           # Combined LC/GC metabolomics (multi-platform)
└── lipid_canonical/          # Lipidomics (LipidSearch results)
```

Each lane directory contains:
- **seed file** — import-tool-specific TSV/TXT (e.g. `report.tsv`, `proteinGroups.txt`)
- **`design.tsv`** — two-column sample-to-group mapping (`sample`, `group`)
- **`proteins.fasta`** — required proteomics FASTA input for browser import workflows

---

## How to Regenerate

Fixtures are idempotent: running the script twice produces byte-identical files.

```r
# From the R console (project root):
source("tests/testdata/e2e/generate_fixtures.R")

# Or via Rscript:
Rscript tests/testdata/e2e/generate_fixtures.R
```

**Requirements:** `jsonlite` must be installed. No other package dependencies.

After regeneration, verify integrity:

```r
devtools::load_all()
source("tests/testthat/helper-e2e-fixture-integrity.R")
validate_e2e_fixtures()
```

---

## Manifest Schema

`manifest.json` top-level structure:

```json
{
  "schema_version": "1.0.0",
  "lanes": [ ... ]
}
```

Each lane object:

| Field                | Type            | Description |
|----------------------|-----------------|-------------|
| `lane_id`            | string          | Unique identifier, used as directory key |
| `omic_type`          | string          | `"proteomics"`, `"metabolomics"`, `"lipidomics"` |
| `workflow_type`      | string / null   | `"DIA"`, `"TMT"`, `"LFQ"`, or `null` |
| `import_tool`        | string          | `"diann"`, `"pd_tmt"`, `"maxquant"`, `"msdial"`, `"lipidsearch"` |
| `use_limpa`          | boolean         | Whether the lane uses the LIMPA DA method |
| `fixture_dir`        | string          | Subdirectory name under `e2e/` |
| `seed_file`          | string          | Filename of the primary import seed file |
| `fasta_file`         | string / null   | Optional FASTA filename for proteomics lanes |
| `assays`             | array           | Platform assay names (`null` for single-assay lanes) |
| `expected_contrasts` | array\<string\> | DA contrasts the E2E suite exercises |
| `report_template`    | string          | Rmd stub filename in `report_templates/` |
| `enrichment_backend` | string          | `"gprofiler"` or `"clusterprofiler"` |
| `expected_exports`   | array\<string\> | Output filenames the E2E suite asserts on |
| `sample_count`       | integer         | Total number of samples across all groups |
| `group_count`        | integer         | Number of experimental groups |
| `groups`             | array\<string\> | Group labels (e.g. `["WT", "KO"]`) |
| `cross_omic_packs`   | array\<string\> | Pack IDs from `cross_omic_packs.json` this lane participates in |

---

## Synthetic Signal

Each lane has **2 DA-significant features** injected:

| Feature | Direction | Effect Size |
|---------|-----------|-------------|
| Feature 1 | KO up-regulated | ~2x fold change |
| Feature 2 | KO down-regulated | ~0.5x fold change |

Remaining features have near-null group differences. This ensures E2E tests can
assert on the presence of known DA hits without hardcoding arbitrary values.

---

## Cross-Omic Pack System

`cross_omic_packs.json` defines multi-omic integration test scenarios (MOFA):

```json
{
  "schema_version": "1.0.0",
  "packs": [
    {
      "pack_id":         "prot_metab",
      "description":     "Proteomics + Metabolomics (LC-MS) MOFA integration",
      "lanes":           ["prot_dia", "metab_lc"],
      "integration_type": "mofa",
      "expected_factors": 2,
      "expected_views":  ["proteomics", "metabolomics"]
    }
  ]
}
```

Lanes declare which packs they participate in via `cross_omic_packs` in `manifest.json`.
The `cross_omic_packs.json` pack `lanes` array must reference only valid `lane_id` values.

Active packs: `prot_metab`, `prot_lipid`, `triple_omic`.

---

## Relationship to E2E Tests

E2E helper files in `tests/testthat/`:

| File | Purpose |
|------|---------|
| `helper-e2e-manifest.R` | `read_e2e_manifest()` — parses and validates manifest paths |
| `helper-e2e-fixture-integrity.R` | `validate_e2e_fixtures()` / `expect_e2e_fixtures_valid()` |
| `helper-e2e-cross-omic.R` | Cross-omic pack loading utilities |
| `helper-e2e-enrichment-doubles.R` | gprofiler2 / clusterProfiler test doubles |
| `helper-e2e-report-doubles.R` | Report template resolution doubles |
| `test-e2e-fixture-integrity.R` | Structural integrity test (runs in CI) |
| `test-e2e-manifest.R` | Manifest schema validation tests |

The integrity check (`test-e2e-fixture-integrity.R`) **does** run in CI and will fail
if any seed file is missing or malformed. The generation script does **not** run in CI.
