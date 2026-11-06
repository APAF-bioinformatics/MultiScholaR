# E2E Test Fixtures

Synthetic, deterministic fixture data for the MultiScholaR end-to-end test suite.
All data is synthetic; no real experimental data is used. `generate_fixtures.R` writes the generated lanes and
checks every committed fixture against the sha256 recorded in `manifest.json` (KIT-101 F-06); the lanes it does not
write are recorded there as hand-made, with the reason.

---

## Fixture Architecture

Fixtures are **manifest-driven**: `manifest.json` is the single source of truth for
all 10 workflow lanes (5 proteomics, 3 metabolomics, 1 lipidomics, 1 transcriptomics). Each lane defines its format,
sample layout, and expected outputs.

```
tests/testdata/e2e/
├── manifest.json             # Lane catalogue (schema_version 1.0.0)
├── cross_omic_packs.json     # Multi-omic integration pack definitions
├── generate_fixtures.R       # Regeneration script (developer tool)
├── report_templates/         # Stub Rmd templates for report rendering tests
│   ├── DIANN_report.rmd
│   ├── DIANN_limpa_report.rmd
│   ├── TMT_report.rmd
│   ├── LFQ_report.rmd
│   ├── metabolomics_report.rmd
│   ├── lipidomics_report.rmd
│   └── transcriptomics_report.rmd
├── prot_dia/                 # DIA-NN proteomics (long format, no LIMPA): lane prot_dia
├── prot_dia_limpa/           # DIA-NN proteomics (long format, with LIMPA): lane prot_dia_limpa
├── prot_tmt_pd/              # TMT proteomics (Proteome Discoverer format): lane prot_tmt
├── prot_lfq_maxquant/        # LFQ proteomics (MaxQuant proteinGroups.txt): lane prot_lfq
├── prot_lfq_fragpipe/        # LFQ proteomics (FragPipe combined_protein.tsv): lane prot_lfq_fragpipe
├── prot_tmt/, prot_lfq/      # retained six-sample proteomics fixtures (no lane; see below)
├── metab_lc/                 # LC-MS metabolomics (custom format, Pos+Neg)
├── metab_gc/                 # GC-MS metabolomics (custom format, GCMS)
├── metab_combined/           # Combined LC/GC metabolomics (multi-platform)
├── lipid_canonical/          # Lipidomics (multi-assay LipidSearch results)
└── trx_screen/               # Transcriptomics counts (KIT-040 trx_screen recipe)
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

**Requirements:** `jsonlite`, `withr` and `digest`. The proteomics lanes read
`tools/profiling/kit_protein_corpus.R` and trx_screen reads `tools/profiling/kit_omics_corpus.R` from the same
source tree. The script stops with an error if any recorded sha256 differs after writing.

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
| `import_tool`        | string          | `"diann"`, `"pd_tmt"`, `"maxquant"`, `"fragpipe"`, `"msdial"`, `"lipidsearch"` |
| `use_limpa`          | boolean         | Whether the lane uses the LIMPA DA method |
| `fixture_dir`        | string          | Subdirectory name under `e2e/` |
| `seed_file`          | string          | Filename of the primary import seed file |
| `fasta_file`         | string / null   | Optional FASTA filename for proteomics lanes |
| `assay2_file`        | string / null   | Optional secondary assay upload for multi-assay metabolomics/lipidomics lanes |
| `assays`             | array           | Platform assay names (`null` for single-assay lanes) |
| `expected_contrasts` | array\<string\> | DA contrasts the E2E suite exercises |
| `report_template`    | string          | Rmd stub filename in `report_templates/` |
| `enrichment_backend` | string          | `"gprofiler"` or `"clusterprofiler"` |
| `expected_exports`   | array\<string\> | Output filenames the E2E suite asserts on |
| `sample_count`       | integer         | Total number of samples across all groups |
| `group_count`        | integer         | Number of experimental groups |
| `groups`             | array\<string\> | Group labels (e.g. `["WT", "KO"]`) |
| `cross_omic_packs`   | array\<string\> | Pack IDs from `cross_omic_packs.json` this lane participates in |
| `fixture_recipe`     | object          | How the lane's files are made: the generator, seeds and shape; `"generator": "none"` with `hand_made` (the reason) for a lane the script does not write |
| `files_sha256`       | object          | sha256 of each fixture file of the lane |

Top level, `retained_fixtures` records files that belong to no lane: the earlier six-sample proteomics fixtures in
`prot_lfq/` and `prot_tmt/` (`regenerated_sha256`, rewritten by their original writers in `generate_fixtures.R`; read
by `tests/testthat/test-workflow-native-session.R` and by the KIT-021 lifecycle members) and the hand-made `seed_*`
files (`hand_made_sha256`).

---

## Synthetic Signal

**Proteomics (KIT-101 revision 2.2).** The five proteomics lanes share one synthetic universe of 400 proteins
(P00001–P00400, genes GENE1–GENE400, tryptic peptides, one FASTA). Each lane's protein abundances come from the KIT-021
protein corpus recipe (`kitProteinMatrix`) at the lane's seed: 4 (DIA) or 5 (TMT, LFQ) replicates per group, and
P00001–P00040 (10%) with a KO effect drawn from N(0, 1.2). With this the DA's FDR can be estimated and some proteins
are truly changed. The DIA-NN reports carry every column the import, the Q-value filter and the precursor rollup read.
See each lane directory's README.

**Lipidomics and metabolomics (hand-made).** Each lane has **2 DA-significant features**:

| Feature | Direction | Effect Size |
|---------|-----------|-------------|
| Feature 1 | KO up-regulated | ~2x fold change |
| Feature 2 | KO down-regulated | ~0.5x fold change |

Remaining features have near-null group differences. This ensures E2E tests can
assert on the presence of known DA hits without hardcoding arbitrary values.

The lipidomics canonical lane is intentionally multi-assay: `LCMS_Pos` and
`LCMS_Neg` are uploaded in the full browser workflow, while `lipidsearch_gcms.txt`
is available for the E2E-009 GCMS-named import smoke.

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
