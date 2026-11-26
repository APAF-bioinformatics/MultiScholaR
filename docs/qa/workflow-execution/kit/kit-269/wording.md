# KIT-269 wording for the owner's pass (VAL-41)

DRAFT, ticket option (b). The owner edits before release.

| Where | Text (exact, as implemented) |
|---|---|
| `study_parameters.txt`, the line after the "Additional Configuration Parameters" dashes when the listing fails (`R/func_general_filemgmt_study_parameter_writer.R`, the handler) | The configuration settings could not be listed here: <reason> |
| NEWS (unreleased, Fixes) | as in `NEWS.md`, from the ticket's draft with the failure sentence unchanged |

`<reason>` is the condition's message with ANSI codes stripped and white space collapsed (one line). For the failure
this ticket removes it would read (by reading, not run): "i In index: 17. Caused by error in `purrr::map_chr()`: i In
index: 1. Caused by error: ! Result must be length 1, not 2." (about 115 characters; in a UTF-8 console the "i" is
"ℹ"). After the fix the line is expected only for a hand-built configuration with an unnamed entry ("subscript out of
bounds"). The owner may prefer the reason without rlang's bullet marks and index lines; that is a one-line change.
