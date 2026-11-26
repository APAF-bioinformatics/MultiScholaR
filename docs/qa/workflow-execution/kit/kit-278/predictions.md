# KIT-278 declaration (committed before the run)

Ticket: tickets/kit/KIT-278.md (git-ignored; main checkout). Test-only declared edit, ruled by the coordinator 2026-11-26.

- Base: 75acd4219cd2d2db9cab9a470147c0e8a479ef73 (origin/v0.9-staging)
- Head: 07e07e67124c81289e37d97feb634ba45a18b886 (the edit: "mod_omic_norm.R" added to allowed_files in tests/testthat/test-omics-module-globals.R, with a KIT-261 comment)
- Files: tests/testthat/test-omics-module-globals.R (the one file)
- Check: none (test-only; R CMD check rides the merge train)
- Ceiling: 30 job-minutes; one run

| File | Base | Head | Verdict |
|---|---|---|---|
| test-omics-module-globals.R | 13 passed, 1 failed: "the super-assignments ... are the listed local closures", setdiff(files_hit, allowed_files) is "mod_omic_norm.R" | 14 passed, 0 failed | fixed |

A missed prediction stops the work and is reported.

## Correction 1 (coordinator ruling on audit K01, 2026-11-26): the allowance is the one line, not the file

The file entry in `allowed_files` (07e07e67) is replaced by the helper's own `allow` argument for exactly the line
`shown_idle <<- idle`: `omics_global_hits(sources, "<<-", allow = "shown_idle <<- idle")`. A future `<<-` elsewhere in
`R/mod_omic_norm.R` (or any scanned file outside `allowed_files`) still fails. Tests only; nothing under `R/`.

- Base: 75acd4219cd2d2db9cab9a470147c0e8a479ef73 (staging, the branch's base)
- Head: the commit that carries this paragraph (the request names its full hash)
- Files: tests/testthat/test-omics-module-globals.R (the one file); no check
- Ceiling: 40 job-minutes; one run (expected about 17)

| File | Base | Head | Verdict |
|---|---|---|---|
| test-omics-module-globals.R | 13 passed, 1 failed (setdiff = "mod_omic_norm.R") | 14 passed, 0 failed | fixed |

A missed prediction is a result: recorded, no repair, stop.
