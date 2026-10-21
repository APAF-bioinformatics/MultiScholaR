# KIT-055: Proteomics report templates: no canned findings, conditional sections, format-aware methods, correct file and folder names

Status 2026-10-21: **completed** (ticket revision 1.0; user approval the same day after reading the KIT-054 real-data report). Implemented at `b6fc60c9`; install `kit_055` = `core/browser-smoke.Po9hwP` (source `acc2ea85`); the user's real-data report re-rendered on that build and placed in their KC folder.

## What changed (templates and tests only; no kernel, record, store or journey output)

All four proteomics templates (`inst/reports/proteomics/{LFQ,TMT,DIANN,DIANN_limpa}_report.rmd`):

- **R-01** The fixed QC narrative, the Comments and the Opinions sections (sentences such as "all samples passed initial QC criteria", "batches merged effectively", "technical reproducibility meets industry standards") are gone. `FormatResults` now describes what Figure 1 contains (the stages and the four panel types) and states that the assessment is the analyst's; `FormatComments`/`FormatOpinions` are replaced by `FormatAnalystInterpretation`, one section "Analyst interpretation" whose body is an explicitly marked field ("[To be completed by the analyst: …]").
- **R-02** "SOP xxx", the sample-receipt bullets, the trailing method bullets and the author "Your_fancy_self" are marked fields ("[to be completed by the analyst]"); the author line reads "Prepared with MultiScholaR; analyst: [to be completed]". Report input fields for these are a later ticket if wanted.
- **R-03** The significance figure reads `Publication_figures/NumSigDaMolecules` (the name the pipeline and the copy helper write) and the DA summary sentence reads `da_<omic>_num_sig_da_molecules.tab` (falling back to the old name). The copy helper's remaining `NumSigDeMolecules` path (`func_general_filemgmt_results_summary_helpers.R:574`) is the lipidomics/metabolomics branch and is untouched.
- **R-04** The enrichment section renders only when enrichment outputs exist (`Publication_figures/Enrichment_Plots` images or files under the pathway directory); otherwise one sentence: "No enrichment analysis was run for this project." Warnings in every figure chunk go to the render log (`message`), never into the report body.
- **R-05** "for the comparisons defined in the contrast table"; "A linear model was fitted"; the supplementary table named `DA_results_<omic>.xlsx` as the copy writes it; the volcano caption "one per comparison".
- **R-06** The LFQ import paragraph names the recorded format (`params$workflow_name`: MaxQuant proteinGroups with LFQ intensities and the removal of contaminant and decoy rows; FragPipe MaxLFQ; Spectronaut; else the format named); the DIA-NN templates name the DIA-NN main report.
- **R-07** The version block is the last section, "Software and build".
- **R-08** The filter rule is stated once: intensities below the recorded percentile of the whole table are treated as missing; within each group a protein fails at the recorded missing percentage; retained in every group or in at most the recorded share of groups; the percentile carries an ordinal suffix.

## Evidence

- Static tests `tests/testthat/test-prot-summary-report-templates.R`: 111 pass (no canned sentence or placeholder in the printed text of any template, the output names, no `cat` of a warning, the conditional enrichment, the heading order, the LFQ format switch).
- Render check [render-check.R](render-check.R) on copies of finished kit projects (the KIT-021 LFQ and TMT fixture lifecycle members, the KIT-036 screen lifecycle DIA-NN candidate), the record and `study_parameters.txt` written from each copy's store, the publication copy run, the working-tree templates rendered twice (without and with an enrichment placeholder image): [render-check-kit-055-tree.json](render-check-kit-055-tree.json), rendered texts beside it. Every assertion true for all three: no canned sentence, no placeholder, no printed warning, Figure 2 caption present with the figure embedded (media 3 without enrichment, 4 with), the enrichment sentence without and the methods with, the format named, `DA_results_proteomics.xlsx`, the analyst field, "Software and build" last, the source commit, the q-value sentence and the percentile sentence still printed from the record.
- The user's real-data report (KIT-054 new lane) re-rendered on the installed build (`kit-054/fill-summary.R`): 383 KB, placed in `results_summary/proteomics/Study_report/` of the copied outputs; not read by the agent (data handling rule).
- `R CMD check --no-manual --no-tests` at `b6fc60c9`: 3 WARNINGs, 3 NOTEs, the baseline ([log](check-b6fc60c9/r-cmd-check.log)). Dependents: `test-general-filemgmt-config-*` (24, 33), `test-lipid-da-export-shared.R` (18), `test-metab-02c-*` (17) pass.
- The KIT-043 render check (`kit-043/report-check.R`) rerun on the installed build with the pre-ticket templates (`88df9314`) as the baseline and the repaired templates as current ([kit043-rerun/](kit043-rerun/), the texts and its JSON; KIT-043's own committed evidence restored unchanged): all eight renders succeed; the with-record texts print the recorded non-default percentile with its ordinal ("below the 2nd percentile of all intensities"), the recorded q-value (0.1) and the source commit; no canned sentence, no printed warning. Two of that script's flags read FALSE because their patterns are the old template's ("Version Information" as the section name, "%dst" as the ordinal), not because a value is missing; the script is KIT-043's evidence and is left as it was.

## Deviations and claim limits

- Interpretation fields are placeholders by design; no report input UI was added (out of scope).
- The DIA-NN limpa template kept its STRING section unchanged apart from the shared repairs.
- The rendered texts kept as evidence are from the non-sensitive fixtures; the real-data render is the user's to read.

## Hand-off

KIT-046 inherits a `DIANN_limpa_report.rmd` free of the canned narrative. The copy helper's `NumSigDe` lipidomics branch is a note for the stage-3 tickets.
