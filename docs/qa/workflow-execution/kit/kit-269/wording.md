# KIT-269 wording for the owner's pass (VAL-41)

DRAFT, ticket option (b). The owner edits before release.

| Where | Text (exact, as implemented) |
|---|---|
| `study_parameters.txt`, the line after the "Additional Configuration Parameters" dashes when the listing fails (`R/func_general_filemgmt_study_parameter_writer.R`, the handler) | The configuration settings could not be listed here: <reason> |
| NEWS (unreleased, Fixes) | as in `NEWS.md`, from the ticket's draft with the failure sentence unchanged |

`<reason>` is the condition's message with ANSI codes stripped and white space collapsed to single spaces (one line). For a
condition raised by `purrr::map_chr` (the case this ticket removes) that message begins with the whole call's source text
before the cause, so a listing that fails inside `purrr` would print a very long line. The plan re-check (R03) specified
this form; shortening it (for example to the text after the last "Caused by error:") is the owner's call and a one-line change.
