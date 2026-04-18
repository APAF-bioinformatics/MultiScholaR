# Import Wave 1 Staging Handover

## Goal

Stage the first proteomics import helper split without touching live `R/`.

This wave intentionally leaves `mod_prot_import_server()` frozen in live code.
The only live import code change kept outside this staging wave is the verified
fallback-config bugfix in `R/mod_prot_import.R`.

## Canonical Inputs

- Refactor protocol:
  [PLAYBOOK.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/PLAYBOOK.md:1)
- Import backlog entry:
  [GOD_MODULE_STABILIZATION_BACKLOG.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/GOD_MODULE_STABILIZATION_BACKLOG.md:128)
- Current import handover context:
  [HANDOVER-wave1.1.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/HANDOVER-wave1.1.md:1)
- New staged manifest:
  [manifest-import-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-import-wave1.yml:1)

## Intent

Split [func_prot_import.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_import.R:1)
first because it is already a top-level helper collection and fits the
manifest/extractor workflow cleanly.

Do not live-split `mod_prot_import_server()` until this staged helper wave is
reviewed and the import wrapper remains covered by:

- [test-prot-01-import.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-01-import.R:1)
- [test-prot-01b-import-detection-characterization.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-01b-import-detection-characterization.R:1)
- [test-prot-01c-import-module-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-01c-import-module-contracts.R:1)
- [test-format-diann.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-format-diann.R:1)

## Staged Output

Manifest verification passed:

```bash
Rscript tools/refactor/verify_refactor.R --manifest tools/refactor/manifest-import-wave1.yml
```

Exact source extraction was written to staging only:

```bash
Rscript tools/refactor/extract_blocks.R \
  --manifest tools/refactor/manifest-import-wave1.yml \
  --write-targets \
  --output-root tools/refactor/staging/import-wave1 \
  --emit-collate tools/refactor/staging/import-wave1/collate-import-wave1.txt
```

Staged files:

- [func_prot_import_detection.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/import-wave1/R/func_prot_import_detection.R:1) `153`
- [func_prot_import_readers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/import-wave1/R/func_prot_import_readers.R:1) `272`
- [func_prot_import_tmt.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/import-wave1/R/func_prot_import_tmt.R:1) `258`
- [func_prot_import_diann_format.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/import-wave1/R/func_prot_import_diann_format.R:1) `173`
- [func_prot_import_organisms.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/import-wave1/R/func_prot_import_organisms.R:1) `201`

All staged `R/*.R` files parse.

Collate output was emitted with the expected contents, but the extractor wrote
the file under a nested path when given a repo-relative destination:

- actual emitted file:
  [collate-import-wave1.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/import-wave1/tools/refactor/staging/import-wave1/collate-import-wave1.txt:1)

Contents:

- `R/func_prot_import_detection.R`
- `R/func_prot_import_readers.R`
- `R/func_prot_import_tmt.R`
- `R/func_prot_import_diann_format.R`
- `R/func_prot_import_organisms.R`

## Next Step

Review the staged helper files against live
[func_prot_import.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_import.R:1).

If the staged split looks correct, the next apply step should be:

1. promote the reviewed staged helper files into live `R/`
2. update `DESCRIPTION` `Collate:` with the reviewed order
3. decide whether `R/func_prot_import.R` becomes a breadcrumb stub or remains as
   the source until source rewriting is applied in a later step
4. rerun the four import tests above
5. only then plan a staged breakup of `mod_prot_import_server()`
