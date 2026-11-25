# KIT-218: user-facing texts (drafts for the owner's wording pass, VAL-41; OD-2)

| Where | Text | Status |
|---|---|---|
| NEWS, `### Fixes`, the `daAnalysisWrapperFunction()` entry | as committed in `NEWS.md` (KIT-218) | draft; the owner sets the final wording |
| Refusal of a `contrasts_tbl` that is not a data frame (`R/func_prot_da_model_wrapper.R:23-26`, class `prot_da_contrasts_invalid`) | "`contrasts_tbl` must be a data frame with one comparison per row in its first column" / "✖ it is <type>" | **open (audit I09):** this wrapper reads `full_format` first when that column is present, so "in its first column" is not exact here. The text is shared with KIT-221 by decision (VAL-193 (1)); change it in both or neither |
