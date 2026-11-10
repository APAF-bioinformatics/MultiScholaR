# KIT-102: the reader survey and the conversion rule

Committed on 2026-11-11 by the audit of the GUI workover series (finding A102-04): the ticket's acceptance box and P-05 cite this file, and it was never committed. The text under "As written in the ticket" is the parity section of `tickets/kit/KIT-102.md` (a git-ignored file), unchanged.

## Corrections (audit A102-01, A102-02; remediation R4)

- In a managed GUI session a character quantity column was converted at the design step with `suppressWarnings(as.numeric())` for every reader, before KIT-102. The statement below that such a column "fails" is wrong for the GUI.
- KIT-102 changed PD-TMT without declaring it: a text "0" became NA; a token in some columns became a pivot error; a token in every column still imported silently.
- Measured through vroom by R4: decimal-comma columns were read 100 times too large, silently, in every reader; a token in rows vroom does not sample became a silent NA in DIA-NN and Spectronaut.
- Since `7bc24709` every proteomics reader checks its quantity columns at import: `NaN`, `Filtered` and `n.q.` are read as NA with a counted notice, and other text refuses the import, naming the column.

## As written in the ticket

**The rule** is metab's today, and becomes shared:
- **Which columns:** a resolved sample column that is character or factor, and whose every value, after `trimws`, is empty or NA or a finite number (`metabImportValuesAreNumericLike`, `mod_metab_import_processing_helpers.R:95`).
- **What happens:** the column is converted with `as.numeric` (`coerceMetabImportWorkflowPayloadSamples`, `:197`).
- **Other columns are not touched,** and each omic's own validation then refuses them as today:
  - lipid: "non-numeric sample columns" (`mod_lipid_import_processing_helpers.R:531-537`);
  - metab: its numeric-like refusal;
  - trx: "counts must be numeric", class `trx_counts_invalid` (`func_trx_counts.R:28`).
- **Sample-column detection does not change.** It is the same in lipid and metab (`resolveLipidImportSampleColumns`, `:101`; `resolveMetabImportSampleColumns`, `mod_metab_import_column_helpers.R:16`).
- **Trx:** the conversion runs before `checkTrxImportCounts`. So `"12"` passes as 12, `"-1"` is refused as negative, and a column with `"x"` is still refused as `trx_counts_invalid` (MCI-026.2 keeps passing).

**Does any proteomics reader leave text-numeric quantity columns? Yes, two of the five.** Evidence, by reading each reader:

| Reader | Today | Leaves text? |
|---|---|---|
| DIA-NN (`importDIANNData`, `func_prot_import_readers.R:10`) | `Precursor.Normalised` is converted with `as.numeric` (`:41-44`). Any other quantity column that is not a plain double is recorded as `needs_numeric` (`mod_prot_import_processing_helpers.R:126-127`) and converted at the design step (`mod_prot_design_builder_server_helpers.R:225-227` → `transformDIADesignTable`, `utils_workflow_dia_state.R:1068`) | no, converted already |
| FragPipe (`:107`) | `as.numeric` after the pivot (`:203`, F-18) | no |
| MaxQuant (`:258`) | `as.numeric` after the pivot (`:330`, F-18) | no |
| Spectronaut (`importSpectronautData`, `:80-97`) | `vroom` guess only. `PG.Quantity`/`PEP.Quantity` are never converted | **yes**: a character quantity column reaches the protein-route design step |
| Proteome Discoverer TMT (`importProteomeDiscovererTMTData`, `func_prot_import_tmt.R:27`) | `readxl::read_excel` (`:51`) or `vroom` (`:60`), then a pivot (`:159`). No conversion; the F-26 zero-to-NA step runs only `if (is.numeric(...))` (`:172`) | **yes**: Excel cells stored as text arrive as character. A mix of text and numeric abundance columns also breaks the pivot |

A character quantity column fails later on the protein routes, in the design step's `log2(.x + 1)` (`mod_prot_design_state_helpers.R:69-76`), with "non-numeric argument". **Decision:**
- **Spectronaut:** the reader applies the shared rule to its quantity column after the read.
- **PD-TMT:** the reader applies it to each `Abundance` column before the pivot, so the F-26 step then always applies.
- **DIA-NN, FragPipe and MaxQuant are unchanged.** They already convert. They turn any non-numeric text into NA silently; that difference is parked (P-05), not changed.
- **The proteomics lanes are not affected:** every committed proteomics fixture has numeric quantity columns, so no proteomics lane or journey changes.
