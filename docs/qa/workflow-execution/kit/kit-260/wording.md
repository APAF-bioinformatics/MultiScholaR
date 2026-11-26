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

## Points from the implementation audit (AUDIT-KIT-260-impl-20261126), for the owner's wording pass

- **I08 (what the note counts).** The `IS-` case counts a name that contains 'IS-' after its first character and that
  no other alternative of the stored pattern matches; a name that starts with 'IS-' (a label the corrected default also
  matches, as `^IS-`) is not counted. Read literally, 'matched only through the unanchored IS-' would count it too (E8's
  object: 10 rows per assay instead of 9).
- **I09 (the import help's last sentence).** "Leave it empty for none." With 'Apply ITSD normalization' ticked (the
  default), an empty box stops the normalisation with an R error ("missing value where TRUE/FALSE needed"; the slot
  becomes NA). Proposed: "Leave it empty for none, and untick 'Apply ITSD normalization' on the normalisation tab." The
  error itself is KIT-242's (its ticket 2.1 (4)).
- **I11 (the note on a new project).** The note fires for any used name containing `w/o MS2:`, whatever the pattern.
  Under the corrected default a labelled standard MS-DIAL prefixed, such as 'w/o MS2:PC 15:0_18:1(d7)' or
  'w/o MS2:Caffeine_d9', is still used (through `\(d[0-9]+\)` and `_d[0-9]+$`) and then gives the sentence "This project
  was imported with an internal-standard pattern that also matches features MS-DIAL annotated without MS/MS ..." on a
  new project: true of the row, within the decided rule, but it reads as a problem of an older project.
- **I12 (wording points).** `used_log` writes "1 rows" for one row; the normalisation tab writes the identifier as
  "IS_Caffeine (Alignment ID 41)" and the integration tab as "Caffeine_d9 (42)"; the record viewer prints outcomes by
  their internal key ("itsd_note = ..."); NEWS's "every such feature" holds for metabolomics, but for lipidomics only
  where the user picked the name column (a real MS-DIAL lipid export gets no annotation column proposed, G0 D2n);
  standards written with a hyphen ('Cholesterol-d7') are still not recognised (the hyphen form was not decided) and no
  text says so.
- **I13 (standards the corrected default no longer finds; decided, informational).** Labels with "IS-" after the start
  ('Leu-IS-1') and MS-DIAL-prefixed standards ('w/o MS2:IS_Caffeine'); an assay whose only standards are written so is
  skipped (KIT-242 shows and records the skip).

## The merged block, as the normalisation tab now logs it (KIT-242 + KIT-260; AUDIT-KIT-260-impl M04)

After the merge with KIT-242 (666ace2e) one reporter, KIT-242's `omicsItsdReportResult`, writes these lines after the ITSD
step, for example with standards matched in assay A1 and none in A2:

```
ITSD normalization complete
Internal standards used per assay: A1 5; A2 0
ITSD, assay A1: 5 rows used as internal standards, then removed: IS_Caffeine (Alignment ID 41); Caffeine_d9 (Alignment ID 42); Caffeine-ISTD (Alignment ID 43); IS-Leucine (Alignment ID 44); PC 15:0_18:1(d7) (Alignment ID 53)
Warning: Internal-standard normalisation was not applied to assay A2: no internal standard matched the pattern '...' in column 'Metabolite name'.
```

and two notifications: KIT-242's warning (only where an assay or a sample was not normalised to standards) and
KIT-260's message "Internal-standard normalisation used 5 (A1) rows as internal standards and removed them; the log names
them." The count of each assay appears twice (KIT-242's count line and KIT-260's rows line); the counts agree by
construction. A single line per assay (the ticket's planned shape, 8.3) is the owner's choice in the wording pass.
