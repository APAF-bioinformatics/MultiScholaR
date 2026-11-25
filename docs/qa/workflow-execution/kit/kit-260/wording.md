# KIT-260: user-facing wording as implemented (VAL-41: implemented as proposed; the owner edits before release)

All texts live in `omicGuiItsdTexts()` (`R/utils_omic_gui_config.R`), except the NEWS entry (`NEWS.md`, Scientific
changes). `%s`/`%d` are filled as stated.

| Key | Where it shows | Text |
|---|---|---|
| `import_help` | Import tab, under "Internal Standard Pattern (Regex)" (lipidomics and metabolomics) | Regular expression, case-sensitive. The normalisation matches it in the annotation column chosen above (the integration also in the identifier column). With 'Apply ITSD normalization' ticked and no standard selected in the ITSD Selection tab, every matching row is used as an internal standard to normalise each sample and is then removed from the data. Leave it empty for none. |
| `used_log` | Normalisation tab log, one line per assay (both routes) | ITSD, assay %s: %d rows used as internal standards, then removed: %s |
| `row` | each row in that line and in the record outcome `itsd_features_per_assay` | %s (%s %s), i.e. annotation (identifier column identifier), e.g. "IS_Caffeine (Alignment ID 41)" |
| `row_no_annotation` | a row without an annotation | %s %s, e.g. "Alignment ID 41" |
| `used_notification` | Normalisation tab notification (type "message"), once per run | Internal-standard normalisation used %s rows as internal standards and removed them; the log names them. (%s e.g. "5 (A1) and 5 (A2)") |
| `legacy_note` | Analysis record (step outcome `itsd_note`, printed by the record viewer), study-parameters file (`* Note:` in the ITSD block), lipidomics and metabolomics reports (Methods, after the ITSD paragraph), integration Sample Matching tab facts | This project was imported with an internal-standard pattern that also matches features MS-DIAL annotated without MS/MS, or names containing 'IS-' in capitals; such features were used as internal standards and removed (%s). (%s e.g. "assay A1, 9 rows; assay A2, 9 rows"; on the integration tab "data set metabolomics_A1, 9 rows") |
| `integration_names` | Integration Sample Matching tab facts, per data set with marked rows | %s: internal standards marked: %s. (rows as "display name (identifier)") |

Notes for the owner's wording pass:
- `legacy_note` starts from the ticket's proposed sentence (OD-3 (a)) and adds the `IS-` case and where; it is shown
  only where the rows a stored pattern matched included a name containing `w/o MS2:` or a name matched only through
  the unanchored `IS-`.
- The NEWS entry is the ticket's section 11 draft with the decided options (OD-1 (c)) written in.
