# Wave 1 Handover

## Goal

Systematically demonolith `R/` into smaller files using manifest-driven source
extraction rather than manual read/rewrite.

Wave 1 scope is the proteomics DA slice only:

- `R/mod_prot_da.R`
- `R/func_prot_da.R`
- `R/da_analysis_function_wrapper.R`

## Invariants

- Do not hand-rewrite moved code bodies.
- Use `tools/refactor/manifest-wave1.yml` as the source of truth.
- Extract exact source ranges.
- Keep the split staged outside `R/` until manual merge review is done.
- Do not rewrite package sources until staged output has been inspected.

## Canonical Artifacts

- Workflow/tooling overview:
  [tools/refactor/README.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/README.md:1)
- Extractor:
  [tools/refactor/extract_blocks.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/extract_blocks.R:1)
- Verifier:
  [tools/refactor/verify_refactor.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/verify_refactor.R:1)
- Wave 1 manifest:
  [tools/refactor/manifest-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-wave1.yml:1)
- Staged Wave 1 output:
  [tools/refactor/staging/wave1/R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/wave1/R)
- Staged collate list:
  [tools/refactor/staging/wave1/collate-wave1.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/wave1/collate-wave1.txt:1)

## What Was Changed

- Added manifest-driven extraction tooling.
- Added `--output-root` staging support so split files can be generated outside
  package `R/`.
- Rewrote filename-coupled test/script logic:
  - [tests/testthat/test-protein-da-golden-master.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-protein-da-golden-master.R:1)
  - [tests/run_mock_de.R](/home/doktersmol/Documents/MultiScholaR/tests/run_mock_de.R:1)

## Current State

- Verifier passed before application:

```bash
Rscript tools/refactor/verify_refactor.R --manifest tools/refactor/manifest-wave1.yml
```

- Staged extraction succeeded:

```bash
Rscript tools/refactor/extract_blocks.R \
  --manifest tools/refactor/manifest-wave1.yml \
  --write-targets \
  --output-root tools/refactor/staging/wave1 \
  --emit-collate collate-wave1.txt
```

- Wave 1 has now been applied into package `R/`:

```bash
Rscript tools/refactor/extract_blocks.R \
  --manifest tools/refactor/manifest-wave1.yml \
  --write-targets \
  --rewrite-sources
cp tools/refactor/staging/wave1/R/*.R R/
```

- `DESCRIPTION` `Collate:` has been updated for the new Wave 1 files.
- The duplicate legacy `writeInteractiveVolcanoPlotProteomicsMain()` was
  removed from `R/da_analysis_function_wrapper.R`.
- Live package `R/` files parse cleanly after the split.

## Known Manual Review Item

One manifest entry was intentionally marked `manual_merge`:

- `compat_write_interactive_main`

Reason:

- `writeInteractiveVolcanoPlotProteomicsMain` exists in both:
  - [R/func_prot_da.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da.R:3621)
  - [R/da_analysis_function_wrapper.R](/home/doktersmol/Documents/MultiScholaR/R/da_analysis_function_wrapper.R:598)

Current decision:

- Treat the `func_prot_da.R` version as canonical.
- Review the compatibility version only for `de_` vs `da_` backward-compat
  semantics before applying the split into package `R/`.

Status:

- Resolved in staged output and applied to live package `R/`.
- The canonical helper in
  [R/func_prot_da_volcano.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_volcano.R:749)
  accepts legacy `de_*` aliases and normalizes to the current `da_*` path.
- The duplicate legacy implementation was not staged into
  `func_prot_da_compat.R`, and it has been removed from
  `R/da_analysis_function_wrapper.R`.

## Environment Blocker

Full `devtools::load_all()`-based test execution is blocked in this environment
 because required package imports are missing:

- `arrow`
- `circlize`
- `ComplexHeatmap`
- `Glimma`
- `limma`
- `qvalue`

This affects:

- [tests/testthat/test-protein-da-golden-master.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-protein-da-golden-master.R:1)
- [tests/run_mock_de.R](/home/doktersmol/Documents/MultiScholaR/tests/run_mock_de.R:1)

The mock script now fails with a clear message instead of relying on dead
 filename-based `source()` calls.

## Exact Next Step

1. Run `devtools::document()` in a dependency-complete environment and inspect
   `NAMESPACE` changes.
2. Run the Wave 1 proteomics DA tests in a dependency-complete environment.
3. Decide whether the now-empty `R/func_prot_da.R` and
   `R/da_analysis_function_wrapper.R` placeholders should remain as migration
   breadcrumbs or be removed in a later cleanup pass.

## Zero-Loss Resume Assessment

For Wave 1: yes, with this file plus the manifest and staged output.

Without this file: almost, but not quite. The tooling and manifest capture the
mechanics, but not the unresolved manual merge, the environment blocker, or the
recommended next action.
