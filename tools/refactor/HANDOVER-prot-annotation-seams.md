# Proteomics Annotation Seam Map

## Goal

Stabilize `func_prot_annotation.R` one bounded checkpoint at a time while preserving the public annotation entry points and existing annotation test behavior.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/func_prot_annotation.R`
- classification: `direct-extraction-ready`
- next step: `Treat the proteomics annotation target as complete and archived; keep the focused annotation gate as the regression surface and do not reopen func_prot_annotation.R unless a real regression appears.`
- latest completed checkpoint: `annotation_wave12_applied`
- latest completed checkpoint summary:
  reviewed and applied the twelfth exact-source UniProt batch-helper wave from
  `tools/refactor/manifest-prot-annotation-wave12.yml` into live `R/` by
  materializing `R/func_prot_annotation_uniprot_batch_helpers.R` (`134`
  lines), removing `subsetQuery()`, `.prepareUniprotBatchInput()`,
  `batchQueryEvidenceHelper()`, `batchQueryEvidence()`,
  `batchQueryEvidenceHelperGeneId()`, and `batchQueryEvidenceGeneId()` from
  `R/func_prot_annotation.R`, emitting live collate artifact
  `tools/refactor/collate-prot-annotation-wave12.txt`, updating
  `DESCRIPTION` collate order so the new helper loads before the public
  annotation file, shrinking `R/func_prot_annotation.R` to `458` lines, and
  rerunning both `tools/refactor/check_wave_apply.R` and the focused
  annotation gate green with `19` passes
- prior completed checkpoint summary:
  reviewed and applied the eleventh exact-source UniProt ID helper wave from
  `tools/refactor/manifest-prot-annotation-wave11.yml` into live `R/` by
  materializing `R/func_prot_annotation_uniprot_id_helpers.R` (`84` lines),
  removing `getUniprotRegexPatterns()`, `normalizeUniprotAccession()`,
  `cleanIsoformNumber()`, and `.cleanProteinIds()` from
  `R/func_prot_annotation.R`, emitting live collate artifact
  `tools/refactor/collate-prot-annotation-wave11.txt`, updating
  `DESCRIPTION` collate order so the new helper loads before the public
  annotation file, shrinking `R/func_prot_annotation.R` to `586` lines, and
  rerunning both `tools/refactor/check_wave_apply.R` and the focused
  annotation gate green with `19` passes
- earlier completed checkpoint summary:
  drafted, verified, and staged the eleventh exact-source UniProt ID helper
  wave from `tools/refactor/manifest-prot-annotation-wave11.yml` without
  rewriting live `R/` by materializing staged
  `tools/refactor/staging/prot-annotation-wave11/R/func_prot_annotation_uniprot_id_helpers.R`
  (`84` lines), emitting staged collate artifact
  `tools/refactor/staging/prot-annotation-wave11/collate-prot-annotation-wave11.txt`,
  keeping `R/func_prot_annotation.R` at `666` lines, and rerunning the
  focused annotation gate green with `19` passes
- oldest retained checkpoint summary:
  reviewed and applied the tenth exact-source full-annotation wave from
  `tools/refactor/manifest-prot-annotation-wave10.yml` into live `R/` by
  materializing `R/func_prot_annotation_uniprot_full_helpers.R` (`352`
  lines), removing `.extractProteinIdFromHeader()` and
  `getUniprotAnnotationsFull()` from `R/func_prot_annotation.R`, emitting
  live collate artifact `tools/refactor/collate-prot-annotation-wave10.txt`,
  updating `DESCRIPTION` collate order so the new helper loads before the
  public annotation file, shrinking `R/func_prot_annotation.R` to `666`
  lines, and rerunning both `tools/refactor/check_wave_apply.R` and the
  focused annotation gate green with `19` passes
- next bounded stop point:
  none; `func_prot_annotation.R` is now the archived breadcrumb shell for the
  completed proteomics annotation split, so future work should stay in the
  dedicated helper files unless the focused annotation gate reports a real
  regression
- active unattended run at documentation update time:
  `6-proteomics-annotation-iter-221`
- active runtime phase at documentation update time:
  `executor`

## Existing Safety Net

- `Rscript tools/test_with_renv.R tests/testthat/test-prot-10-annotation.R`

## Notes

- `R/func_prot_annotation.R` is now `458` lines and acts as the archived
  breadcrumb shell for the completed annotation split rather than a live
  helper surface.
- The first exact-source wave is now live in:
  `R/func_prot_annotation_go_helpers.R` (`85` lines) and
  `R/func_prot_annotation_matching_helpers.R` (`318` lines).
- The second exact-source wave is now live in:
  `R/func_prot_annotation_fasta_helpers.R` (`104` lines) with manifest
  `tools/refactor/manifest-prot-annotation-wave2.yml` and live collate
  artifact `tools/refactor/collate-prot-annotation-wave2.txt`.
- The third exact-source wave is now live in:
  `R/func_prot_annotation_fasta_processing_helpers.R` (`265` lines) with
  manifest `tools/refactor/manifest-prot-annotation-wave3.yml` and live
  collate artifact `tools/refactor/collate-prot-annotation-wave3.txt`.
- The fourth exact-source wave is now live in:
  `R/func_prot_annotation_accession_ranking_helpers.R` (`305` lines) with
  manifest `tools/refactor/manifest-prot-annotation-wave4.yml` and live
  collate artifact `tools/refactor/collate-prot-annotation-wave4.txt`.
- The fifth exact-source wave is now live in:
  `R/func_prot_annotation_protein_cleaning_helpers.R` (`226` lines) with
  manifest `tools/refactor/manifest-prot-annotation-wave5.yml`, live collate
  artifact `tools/refactor/collate-prot-annotation-wave5.txt`, and staged
  reference artifacts retained under
  `tools/refactor/staging/prot-annotation-wave5/`.
- The sixth exact-source wave is now live in:
  `R/func_prot_annotation_phosphosite_helpers.R` (`98` lines) with manifest
  `tools/refactor/manifest-prot-annotation-wave6.yml`, live collate artifact
  `tools/refactor/collate-prot-annotation-wave6.txt`, and staged reference
  artifacts retained under `tools/refactor/staging/prot-annotation-wave6/`.
- The seventh exact-source wave is now live in:
  `R/func_prot_annotation_processing_helpers.R` (`163` lines) with manifest
  `tools/refactor/manifest-prot-annotation-wave7.yml`, live collate artifact
  `tools/refactor/collate-prot-annotation-wave7.txt`, and staged reference
  artifacts retained under `tools/refactor/staging/prot-annotation-wave7/`.
- The eighth exact-source wave is now live in:
  `R/func_prot_annotation_uniprot_helpers.R` (`481` lines) with manifest
  `tools/refactor/manifest-prot-annotation-wave8.yml`, live collate artifact
  `tools/refactor/collate-prot-annotation-wave8.txt`, and staged reference
  artifacts retained under `tools/refactor/staging/prot-annotation-wave8/`.
- The ninth exact-source wave is now live in:
  `R/func_prot_annotation_ensembl_helpers.R` (`285` lines) with manifest
  `tools/refactor/manifest-prot-annotation-wave9.yml`, live collate artifact
  `tools/refactor/collate-prot-annotation-wave9.txt`, and staged reference
  artifacts retained under `tools/refactor/staging/prot-annotation-wave9/`.
- The tenth exact-source wave is now live in:
  `R/func_prot_annotation_uniprot_full_helpers.R` (`352` lines) with manifest
  `tools/refactor/manifest-prot-annotation-wave10.yml`, live collate artifact
  `tools/refactor/collate-prot-annotation-wave10.txt`, and staged reference
  artifacts retained under `tools/refactor/staging/prot-annotation-wave10/`.
- The eleventh exact-source wave is now live in:
  `R/func_prot_annotation_uniprot_id_helpers.R` (`84` lines) with manifest
  `tools/refactor/manifest-prot-annotation-wave11.yml`, live collate artifact
  `tools/refactor/collate-prot-annotation-wave11.txt`, and staged reference
  artifacts retained under `tools/refactor/staging/prot-annotation-wave11/`.
- The twelfth exact-source wave is now live in:
  `R/func_prot_annotation_uniprot_batch_helpers.R` (`134` lines) with
  manifest `tools/refactor/manifest-prot-annotation-wave12.yml`, live collate
  artifact `tools/refactor/collate-prot-annotation-wave12.txt`, and staged
  reference artifacts retained under
  `tools/refactor/staging/prot-annotation-wave12/`.
- `DESCRIPTION` now collates the live annotation helper files, including
  `R/func_prot_annotation_uniprot_batch_helpers.R`, before
  `R/func_prot_annotation.R`, and `tools/refactor/check_wave_apply.R` passed
  after the wave-12 rewrite.
- Treat `tools/refactor/HANDOVER-prot-annotation-seams.md` as complete and
  archived; future annotation maintenance should run against the focused
  annotation gate and touch the dedicated helper files instead of rebuilding
  `R/func_prot_annotation.R`.
- No loop-control artifacts were changed in this iteration.
- April 14, 2026 control-plane hardening note:
  the overnight annotation failure was not in annotation code. The failing
  executor run died during `codex exec` session initialization with
  `Failed to create session: Read-only file system`. The hardened runner now
  launches child Codex sessions from an isolated temp cwd with a sanitized env
  and one retry for known session-init markers, and the annotation lane was
  reopened cleanly on the same target afterward.
