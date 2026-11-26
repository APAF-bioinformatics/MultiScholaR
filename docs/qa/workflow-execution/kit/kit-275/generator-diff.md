# generator-diff.md: sim-k275.R against sim-k188.R at 0866982f

Pinned original: `git show 0866982f:docs/qa/workflow-execution/kit/kit-188/sim-k188.R` (2,283 lines). Bound commit for the package: the first staging commit holding KIT-255 with its I01 repair, `ec27f69e0` (merge aa119308f), checked in `g0.txt` and `step0-lines.txt`. This file is produced from `python3 -I e0-diff.py --patch generator-diff.patch` (the check of E0: the committed patch is the complete diff of the retained text; `e0-diff.py` stops when a hunk of it carries no `K275-` mark outside the header comment) and lists every hunk of the retained text, old lines at 0866982f against new lines of `sim-k275.R`.

## What "copy" means here, stated plainly

The ticket (5.1, 8.1) has the generator and the journeys copied unchanged except the marked lines. KIT-188's file also holds the proteomics journeys, KIT-188's own estimands, checks, units, plants and job mechanics, which this ticket does not use. They are **not carried**: five line ranges of 0866982f are removed (table B) and the added K275 section replaces them (table C). What remains of the original is byte-equal to 0866982f outside the hunks below. This is one reading of "copy" the ticket does not spell out; it is listed so the audit can accept or reject it.

## A. The ticket's marks and the marks added beyond the ticket's list

| Mark | In the ticket's 5.1/5.3 list? | Hunk (old lines at 0866982f) |
|---|---|---|
| K275-seed | yes (:166) | 166 |
| K275-sig | yes (:164) | 164 |
| K275-cells | yes (cellDef :62, :64, CELLS after :80) | 62, 64, 80 |
| K275-C | yes (:231) | 231 |
| K275-noIS | yes (:183) | 183 |
| K275-plant | yes (after :285, 5.3) | insert after 285 |
| K275-ret | yes (:303-308; here the list's last line, :308) | 308 |
| K275-mask | yes (`maskOf` :315, :319) | 315, 319 |
| K275-W | yes (:690-691) | 691 |
| K275-markers | yes (`ARM_MARKERS` :137-139) | 139 |
| **K275-guard** | **no: the lead's instruction (lesson 2): the partial-match warning is an error, also inside `quiet()`** | header (the shared `guardWarning`, the global net and the child's error writer, new lines 18-43), `quiet()`, and (smoke 3 repair) `journeyL`'s `tryCatch(withCallingHandlers(..., warning = guardWarning))` and its handler |
| **K275-prefix** | **no: the same lesson. `out$da` read when absent matched `da_reached`; `out$e1` could match `e12`, `e13` (KIT-188's smoke 2 found it). The journey's fields `e12`, `e13`, `da_reached` are renamed `itsd_check`, `norm_check`, `reached_da`, as KIT-188 did** | 558-559, 582, 659-661, 677-679, 708 |
| **K275-pairs** | **no, but 8.2 of the ticket names these reads ("the S4 states before and after RUV passed to `pearsonCorForSamplePairs`") and E5 needs the filter's own removed set: three added blocks in `journeyL`** | inserts after 672, 695 and 704 |
| **K275-REFUSALS** | **no: two refusal classes of the control rule that exist at the bound commit (`multischolar_refusal`, `multischolar_ruv_grouping_unusable`) added to `REFUSALS` (the comment sits at the end of the statement, audit S15), as KIT-188's later version did, so a refused journey is not recorded `crashed`** | 83 |
| **K275-assaycols** | **no (audit S01): `ruv_n_controls` guarded with `is.logical` (the optimiser can return a list when an assay is skipped)** | 595 |
| **K275-nogroup** | **no (audit S03, ruling R-1): a design that lost a whole group by the filter is `no_da`, class `filter_removed_group:<groups>`; the DA block is wrapped in an else** | after 723 |
| **K275-E3** | **no: `checkE3` reduced to the ticket's E3 (the route defaults of 5.2 and the markers). KIT-188's clauses on the ITSD import pattern (KIT-260, not on staging: the pattern is still `w/o MS2:`) and on labelled names matching an ITSD pattern are not carried; ITSD is off in every journey here** | 430-453 |

The header comment (lines 1-40) is replaced by KIT-275's.

## B. Blocks of 0866982f not carried (omitted before the comparison)

| Old lines | Content |
|---|---|
| 86-93 | `BLIND`: KIT-188's blinded estimand columns |
| 809-969 | Part P: `pConfig`, `REQUESTS_P`, `protMatrixOf`, `journeyP` |
| 978-1079 | `scoreFamily`, `ASSAY_COLUMNS`, `ROW_COLUMNS`, `newRow`, `armRows` |
| 1105-1161 | `checkE9L` (the direct-call idea of its LP branch, :1136-1148, is reused in the new `checkE12`) |
| 1180-2283 | `anovaP` to `checkE4`, `ARMS_OF`, the `runCondition*` functions, `k172cSim`, `checkE7Refusal`, `checkE12Patterns`, the planted cases, `runUnit`, the child and the parent |

Kept unchanged from those regions: `groupOfComparison` (:970-977), `choiceOf` (:1080-1089), `checks_env`/`addCheck` (:1091-1098), `positional` (:1099-1104, comment lines 1099-1100 kept), `omicObjectFromMatrices` (:1162-1179).

## C. Added (not copied): the K275 section

`# K275: everything below is added (not copied)` to the end of `sim-k275.R`: `TIER_PLAN` and the tier tables, `armSpec`, `armClass`, `requestK275`, `maskedLog2`, `dropSamples`, `makeVersion`, `plantedShares`, `bestCorrelation`, `failsAt`, `failList`, `pseudoReplicateRatio`, `kappaControls`, `unmoderatedReading`, `ROW_COLUMNS`, `newRow`, `rowFrame`, `daFamily`, `daNumbers`, `journeyRow`, `lxRow`, `appendTable`, `flushChecks`, `checkE12`, `runArm`, `runEmulation`, `checkVersionsE1`, `runUnit`, `runChildMain`, `runParentMain`, and the library switch (`options(k275.library = TRUE)` loads the functions only, for `fixture-k275.R`).

## D. The hunks, line by line

### Hunk 1: header; old 1-3

Old (0866982f):
```r
1: # KIT-188 (VAL-72, VAL-81, VAL-84): calibration of the default route with legacy RUV-III Automatic, lipidomics,
2: # metabolomics (Part L) and proteomics (Part P). Evidence code under docs/qa: no package code. Runs on a GitHub-hosted
3: # runner with MultiScholaR installed from the bound commit (run-tests `scripts`), from the package root.
```
New (sim-k275.R):
```r
1: # KIT-275 (VAL-330, VAL-344, VAL-347): the sample-correlation filter after RUV and one anomalous sample (lipidomics, with a
2: # metabolomics identity check). A measurement; evidence code under docs/qa; no package code. Runs on a GitHub-hosted runner with
3: # MultiScholaR installed from the bound commit (run-tests `scripts`), from the package root. The generator and the journeys are a
4: # copy of KIT-188's sim-k188.R at 0866982f with the marked changes of generator-diff.md (marks K275-*); nothing else of it is carried.
```

### Hunk 2: header; old 6-10

Old (0866982f):
```r
6: #   Rscript sim-k188.R <part> <cells> <first> <last> <out dir>       the job (parent)
7: #   Rscript sim-k188.R --unit <part> <cell> <r> <job> <scenario> <mask> <rds>   one unit (child process)
8: #   part   L (lipidomics and metabolomics: L12, L64, M12, M168) or P (proteomics: P0, P1, P16, P2, P3, P4b)
9: #   cells  one cell or a comma list
10: #   first, last  data sets (S: 1 to 400; lifecycle: 401 to 403)
```
New (sim-k275.R):
```r
7: #   Rscript sim-k275.R <tiers> <design> <first> <last> <out dir>       the job (parent)
8: #   Rscript sim-k275.R --unit <tier> <design> <r> <rds>                one data set of one tier (child process)
9: #   tiers   one tier or a comma list of D, K, F, S, V, M, T, E, ID (tier V waits for KIT-242, section 13 of the ticket)
10: #   design  D3, D4 or D6 (ID: D4); first, last  data sets (counted tiers 1 to 200 or 1 to 100; ID 1 to 20; L1 and the smoke 9101 to 9103)
```

### Hunk 3: K275-guard; old 12-40

Old (0866982f):
```r
12: #   K188_JOBS      comma list of job types (Part L default deciding,reported,munit,ba; Part P deciding,reported,ba):
13: #                  deciding (LA, GA, NO, OR), reported (LA, OR, LM, FP, SC, LCm in Part L; LM, OR in Part P),
14: #                  munit (Part M unit check, L12 and M12), ba (mask BA, L64, M168, P4b), lp (LP: the LA request on the
15: #                  bound commit with KIT-190's markers at their earlier tokens, 5.3; r1.21), las (LAs: the LA request with
16: #                  KIT-224's marker at its earlier token, 34), l1c (the planted cases of 35.6, L12, one data set, K188_PLANT).
17: #                  L24 (38.5, reported) runs the deciding job type only. A job type that does not apply to a cell is not run.
18: #   K188_PLANT     L1c: comma list of the planted cases (itsd, itsdpos, emptysample, ruvskip, ruvdrop, daerror, normplant,
19: #                  exitplant; default all), one child process each; needs K188_JOBS=l1c and K188_LIFECYCLE=1.
20: #   K188_PREDICTIONS  L1c and E7 (R-marked predictions of 35.6 and 7.4, read at the bound commit): absolute path of the TSV
21: #                  (columns key = "<plant or check>/<field>", expected) written into lifecycle.md before the run.
22: #   K188_LIFECYCLE 1: blind lifecycle (L1, L1p): every estimand column is written NA; the rows' SHA-256 is written; the
23: #                  checks E1 to E4, E6, E7, E9 (and E8, E2, E5r) are run and written to checks-<job>.tsv.
24: #   K188_FAULT     1: E7 (crash): one deliberate error in the first unit's first journey.
25: #   K188_WORKERS   child processes at a time (default 4).
26: #   K188_P_REPS    Part P lifecycle: comma list cell=last data set (default "P0=403,P2=403"; others run first only).
27: #   K188_JOB       a name for the job's files (default from part, cells and range).
28: #
29: # One data set (Part L: per job type) or one (data set, scenario, mask) (Part P) per Rscript child process, four at a
30: # time, so that no session owner (duckdb) is forked. Every random draw of a data set is made in generate() before any
31: # package function is called; seeds set.seed(188000000 + 100000 c + 10000 s + r, kind = "Mersenne-Twister",
32: # normal.kind = "Inversion", sample.kind = "Rejection"), c the cell (L12 1, L64 2, M12 3, M168 4, P0 5, P1 6, P16 7,
33: # P2 8, P3 9, P4b 10), s the scenario (N0 1, N1 2, E0 3, E1 4, N2 5, E2 6), r the data set.
34: # The design, estimands and reading rules are KIT-188 revision 1.4, sections 5, 7, 18 and the declaration.
35: # Construction choices this file fixes (5.2 leaves them open): each unwanted component is centred with mean square 1
36: # over the runs (population variance); components are made orthogonal to each other in the order run factor or
37: # plate/batch basis, drift, sample factor; the drift's order is, of 2,000 random orders (within plate in L64), the one
38: # whose drift has the smallest share of variance explained by the groups, and is then residualised; labelled controls
39: # draw their levels and variances as features (the variance scaled by the features' median); every draw is made in
40: # every scenario (a scenario without unwanted variation multiplies it by 0).
```
New (sim-k275.R):
```r
12: #   K275_WORKERS   child processes at a time (default 4)
13: #   K275_LIFECYCLE 1: blind lifecycle (L1): every estimand column of the rows is written NA, the rows' SHA-256 is written, the
14: #                  checks E1 to E12 are run and written
15: #   K275_REPEAT    lifecycle: the tiers whose first data set is run a second time for E9 (default all; the smoke send sets one small tier)
16: #   K275_FAULT     1: E7: one deliberate error in the first journey of the first unit, and the refused journey of E7 (GK3 at D3 on VD)
17: #   K275_SMOKE     1: the smoke send's flag; written to the provenance and to every row; the scorer refuses any input carrying it for an
18: #                  outcome line; a counted request that sets it (without K275_LIFECYCLE=1) is refused here (SMOKE_counted_refusal)
19: #   K275_ONLY_UNITS counted mode: the unit tags (<tier>-<design>-<r>) a continuation runs; refused in the lifecycle; the job name gains -cont
20: #   K275_JOB       a name for the job's files (default from the tiers, design and range)
21: # One data set per child process, all its versions and arms; no session owner (duckdb) is forked. Every random draw of a data set is
22: # made in generate() before any package function is called; seeds in section 5.3 of the ticket.
23: 
24: # Partial-match guard (lesson of KIT-188's smoke checks 1 to 3): `$` on a list or data frame matches a name by its prefix when none is
25: # exact; that warning is an error here, in every process and inside `quiet()` (which would otherwise muffle it). A stop from it is an
26: # instrument defect by definition (VAL-54). Every read of a field that may be absent below uses `[[` with exact matching.
27: options(warnPartialMatchDollar = TRUE) # K275-guard
28: # A `stop()` raised from a GLOBAL calling handler runs with every local handler disabled (R's rule for the handler stack), so it cannot reach a
29: # `tryCatch`: it ends the process with nothing written (smoke 2's planted guard unit). The guard is therefore also installed as a LOCAL calling
30: # handler (`guardWarning`) around each journey and around each unit, where its stop reaches the enclosing `tryCatch`; the global handler stays as
31: # a last net, with an error writer that leaves a trace before the process ends.
32: guardWarning <- function(w) { # K275-guard
33:     if (grepl("partial match of", tolower(conditionMessage(w)), fixed = TRUE)) {
34:         options(k275.guard_trace = paste(vapply(sys.calls(), function(cl) substr(paste(deparse(cl)[[1L]], collapse = ""), 1L, 80L), ""), collapse = " <- "))
35:         stop(paste("partial-match guard:", conditionMessage(w)), call. = FALSE)
36:     }
37: }
38: globalCallingHandlers(warning = guardWarning) # K275-guard
39: if (identical(commandArgs(TRUE)[1L], "--unit")) options(error = function() { # K275-guard   # an error that escaped every handler: write what is known, then end non-zero
40:     dir <- dirname(Sys.getenv("K275_CHECKS_PATH"))
41:     tag <- sub("\\.tsv$", "", sub("^checks-", "", basename(Sys.getenv("K275_CHECKS_PATH"))))
42:     try(writeLines(c(geterrmessage(), getOption("k275.guard_trace", "")), file.path(dir, sprintf("unit-trace-%s.txt", tag))), silent = TRUE)
43:     quit(save = "no", status = 1L)
44: })
```

### Hunk 4: K275-cells; old 62-62

Old (0866982f):
```r
62: cellDef <- function(idx, part, omic, assays, per_assay, groups, per, level, sd, unwanted, blocks = 0L, n2 = FALSE) {
```
New (sim-k275.R):
```r
66: cellDef <- function(idx, part, omic, assays, per_assay, groups, per, level, sd, unwanted, blocks = 0L, n2 = FALSE, no_is = FALSE) { # K275-cells
```

### Hunk 5: K275-cells; old 64-64

Old (0866982f):
```r
64:         sd = sd, unwanted = unwanted, blocks = blocks, n2 = n2)
```
New (sim-k275.R):
```r
68:         sd = sd, unwanted = unwanted, blocks = blocks, n2 = n2, no_is = no_is) # K275-cells
```

### Hunk 6: K275-cells; old 80-80

Old (0866982f):
```r
80:     K172c = cellDef(0L, "L", "lipid", "LCMS_Pos", 515L, 2L, 6L, NA, NA, "none"))
```
New (sim-k275.R):
```r
84:     K172c = cellDef(0L, "L", "lipid", "LCMS_Pos", 515L, 2L, 6L, NA, NA, "none"),
85:     # K275-cells: the designs of KIT-275 (2 groups of 3, 4, 6; two assays of 250 features; one unwanted component; no internal standards), the
86:     # metabolomics identity cell D4M (the same seed index as D4, so the draws are the same), and the exactness fixtures of G1 (one assay, 2,000 features)
87:     D3 = cellDef(1L, "L", "lipid", c("LCMS_Neg", "LCMS_Pos"), 250L, 2L, 3L, 20, 2, "drift_one", no_is = TRUE),
88:     D4 = cellDef(2L, "L", "lipid", c("LCMS_Neg", "LCMS_Pos"), 250L, 2L, 4L, 20, 2, "drift_one", no_is = TRUE),
89:     D6 = cellDef(3L, "L", "lipid", c("LCMS_Neg", "LCMS_Pos"), 250L, 2L, 6L, 20, 2, "drift_one", no_is = TRUE),
90:     D4M = cellDef(2L, "L", "metab", c("LCMS_Neg", "LCMS_Pos"), 250L, 2L, 4L, 20, 2, "drift_one", no_is = TRUE),
91:     X3 = cellDef(5L, "L", "lipid", "LCMS_Pos", 2000L, 2L, 3L, 20, 2, "drift_one", no_is = TRUE),
92:     X4 = cellDef(6L, "L", "lipid", "LCMS_Pos", 2000L, 2L, 4L, 20, 2, "drift_one", no_is = TRUE),
93:     X6 = cellDef(7L, "L", "lipid", "LCMS_Pos", 2000L, 2L, 6L, 20, 2, "drift_one", no_is = TRUE))
```

### Hunk 7: K275-REFUSALS; old 83-83

Old (0866982f):
```r
83: REFUSALS <- c("omics_ruv_insufficient_controls", "prot_ruv_insufficient_controls", "prot_ruv_empty_result",
```
New (sim-k275.R):
```r
96: REFUSALS <- c("multischolar_refusal", "multischolar_ruv_grouping_unusable", "omics_ruv_insufficient_controls", "prot_ruv_insufficient_controls", "prot_ruv_empty_result", # K275-REFUSALS
```

### Hunk 8: K275-guard; old (insert after 108)

Old (0866982f):
```r
(nothing)
```
New (sim-k275.R):
```r
114:         guardWarning(w) # K275-guard
```

### Hunk 9: K275-markers; old 139-139

Old (0866982f):
```r
139:     LAs = list(MultiScholaR.ruv_feature_levels = "shifted"))
```
New (sim-k275.R):
```r
145:     LAs = list(MultiScholaR.ruv_feature_levels = "shifted"),
146:     # K275-markers: the control selection before KIT-190 (positional controls; group-blind), for LG1 and LG2
147:     LG1 = list(MultiScholaR.ruv_control_grouping = "dropped"),
148:     LG2 = list(MultiScholaR.ruv_control_grouping = "dropped"))
```

### Hunk 10: K275-sig; old 164-164

Old (0866982f):
```r
164: generate <- function(cellname, scenario, r) {
```
New (sim-k275.R):
```r
173: generate <- function(cellname, scenario, r, anomaly = NULL) { # K275-sig
```

### Hunk 11: K275-seed; old 166-166

Old (0866982f):
```r
166:     set.seed(188000000 + 100000 * p$idx + 10000 * match(scenario, SCENARIOS) + r, kind = "Mersenne-Twister",
```
New (sim-k275.R):
```r
175:     set.seed(275000000 + 100000 * p$idx + 10000 * match(scenario, SCENARIOS) + r, kind = "Mersenne-Twister", # K275-seed
```

### Hunk 12: K275-noIS; old 183-183

Old (0866982f):
```r
183:     n_is <- if (identical(omic, "lipid")) 5L * n_assays else 0L
```
New (sim-k275.R):
```r
192:     n_is <- if (identical(omic, "lipid") && !isTRUE(p$no_is)) 5L * n_assays else 0L # K275-noIS
```

### Hunk 13: K275-C; old 231-231

Old (0866982f):
```r
231:     C <- switch(p$unwanted, drift_factor = 2L, plate_drift = 2L, batch_drift = p$blocks)
```
New (sim-k275.R):
```r
240:     C <- switch(p$unwanted, drift_one = 1L, drift_factor = 2L, plate_drift = 2L, batch_drift = p$blocks) # K275-C
```

### Hunk 14: K275-plant; old (insert after 285)

Old (0866982f):
```r
(nothing)
```
New (sim-k275.R):
```r
295:     # K275-plant (section 5.3 of the ticket): the planting of one anomalous sample (two for V2) in assay LCMS_Pos, after every draw of
296:     # generate(); its own stream, the session's random-number state restored, so no draw of any version differs from V0 except the planted entries
297:     planted_delta <- NULL
298:     if (!is.null(anomaly)) {   # anomaly = list(form = "g" | "p" | "s" | "t", samples = c(n_runs) or c(n_runs, per))
299:         old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
300:         set.seed(276000000 + 100000 * p$idx + r, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
301:         z_feat <- list(stats::rnorm(n_feat), stats::rnorm(n_feat)); z_lbl <- list(stats::rnorm(n_lbl), stats::rnorm(n_lbl))  # drawn in this order, all, always
302:         u_feat <- stats::runif(n_feat)   # drawn last, all, always: the sparse form's choice of features
303:         cm <- if (identical(anomaly$form, "p")) 0.5 else 0
304:         planted_delta <- list()
305:         for (s in seq_along(anomaly$samples)) {
306:             col <- anomaly$samples[[s]]
307:             rf <- which(feat_assay == "LCMS_Pos"); rl <- which(lbl_assay == "LCMS_Pos")
308:             d_feat <- (p$level + cm * (mu[rf] - p$level) + sqrt(1 - cm^2) * p$sd * z_feat[[s]][rf]) - mu[rf]
309:             if (identical(anomaly$form, "t")) d_feat[!(seq_along(rf) %in% order(u_feat[rf])[seq_len(round(0.1 * length(rf)))])] <- 0
310:             y_feat[rf, col] <- y_feat[rf, col] + d_feat
311:             d_lbl <- if (!anomaly$form %in% c("s", "t")) (p$level + cm * (mu_lbl[rl] - p$level) + sqrt(1 - cm^2) * p$sd * z_lbl[[s]][rl]) - mu_lbl[rl] else numeric(length(rl))
312:             y_lbl[rl, col] <- y_lbl[rl, col] + d_lbl
313:             planted_delta[[s]] <- c(stats::setNames(d_feat, feat_ids[rf]), stats::setNames(d_lbl, lbl_ids[rl]))
314:         }
315:         if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()) else rm(".Random.seed", envir = globalenv())
316:     }
```

### Hunk 15: K275-ret; old 308-308

Old (0866982f):
```r
308:         U = U, order_of = order_of, ba_u = ba_u, ba_block = ba_block, Z = Z, C = C)
```
New (sim-k275.R):
```r
339:         U = U, order_of = order_of, ba_u = ba_u, ba_block = ba_block, Z = Z, C = C, planted_delta = planted_delta) # K275-ret
```

### Hunk 16: K275-mask; old 315-315

Old (0866982f):
```r
315: maskOf <- function(sim, mask) {
```
New (sim-k275.R):
```r
346: maskOf <- function(sim, mask, lod_values = sim$y[sim$feat_ids, , drop = FALSE]) { # K275-mask: the limit of detection is solved on `lod_values` (V0's values)
```

### Hunk 17: K275-mask; old 319-319

Old (0866982f):
```r
319:     out <- sim$U < stats::plogis(solveLod(y, rate) - y)
```
New (sim-k275.R):
```r
350:     out <- sim$U < stats::plogis(solveLod(lod_values, rate) - y) # K275-mask
```

### Hunk 18: K275-E3; old 430-432

Old (0866982f):
```r
430: # 5.2 (r1.20; VAL-241 OD-1 (c)), the import's default pattern as KIT-260 ships it, read at the bound commit by E3 (a difference
431: # stops: the route under test changed)
432: ITSD_DEFAULT_PATTERN <- "^IS_|_d[0-9]+$|ISTD|^IS-|\\(d[0-9]+\\)"
```
New (sim-k275.R):
```r
461: # K275-E3: the ticket's E3 is the route defaults of 5.2 only (ITSD is off in every journey here, so KIT-188's clauses on the ITSD
462: # import pattern and on labelled names matching an ITSD pattern are not carried: at the bound commit the import pattern is still the one
463: # before KIT-260, "w/o MS2:"); the replay markers hold their current token (a leak stops the unit).
```

### Hunk 19: K275-E3; old 438-444

Old (0866982f):
```r
438:     patterns <- c("(?i)itsd|istd|internal.?standard", "(?i)^rt_|_rt$|retention.?time", "(?i)^is_|_is$|^isd_|_isd$",
439:         "(?i)_d[0-9]+$|_13c|_15n|deuterated|labeled")
440:     names_checked <- c(sprintf("LBL_%03d", 1:10), sprintf("CTL_%03d", 1:5))
441:     matched <- names_checked[vapply(names_checked, function(x) any(vapply(patterns, grepl, logical(1L), x = x, perl = TRUE)), logical(1L))]
442:     # r1.21 (R20-P06): the cyclic-loess minimum's state is recorded and either state passes; every marker holds its current
443:     # token (40.3 (5)); 5.2 (r1.20, VAL-241 OD-1 (c)): the import's default ITSD pattern as KIT-260 ships it
444:     loess_minimum <- tryCatch(omicsCyclicLoessMinimumActive(), error = function(e) NA)
```
New (sim-k275.R):
```r
(nothing)
```

### Hunk 20: K275-E3; old 446-453

Old (0866982f):
```r
446:     imported <- tryCatch(if (identical(omic, "metab")) getMetabolomicsColumnDefaults("msdial")$is_pattern else
447:         getLipidomicsColumnDefaults("msdial")$is_pattern, error = function(e) NA_character_)
448:     pattern_ok <- identical(omic, "prot") || identical(imported, ITSD_DEFAULT_PATTERN)
449:     got$cyclic_loess_minimum <- loess_minimum
450:     got$itsd_import_pattern <- imported
451:     list(pass = !length(differ) && !length(matched) && !length(leaked) && pattern_ok, got = got,
452:         note = sprintf("differ: %s; labelled names matching an ITSD pattern: %s; markers off their current token: %s; import pattern as shipped %s",
453:             paste(differ, collapse = ",") %||% "", paste(matched, collapse = ","), paste(leaked, collapse = ","), pattern_ok))
```
New (sim-k275.R):
```r
470:     got$cyclic_loess_minimum <- tryCatch(omicsCyclicLoessMinimumActive(), error = function(e) NA) # recorded, either state passes (as KIT-188)
471:     need <- c("ruv_control_grouping", "ruv_k_bound", "ruv_feature_levels", "correlation_exclusion", strsplit(Sys.getenv("K275_REQUIRED_MARKERS"), ",", fixed = TRUE)[[1L]]) # K275-E3
472:     absent <- setdiff(need[nzchar(need)], names(workflowReplayRuleMarkers)) # K275-E3
473:     ce_ok <- identical(workflowReplayRuleMarkers[["correlation_exclusion"]][["current"]], "within_assay") # K275-E3
474:     got$correlation_exclusion <- workflowReplayRuleMarkers[["correlation_exclusion"]][["current"]] # K275-E3
475:     list(pass = !length(differ) && !length(leaked) && !length(absent) && ce_ok, got = got, # K275-E3
476:         note = sprintf("differ: %s; markers off their current token: %s; markers absent: %s; correlation_exclusion within_assay: %s", paste(differ, collapse = ","),
477:             paste(leaked, collapse = ","), paste(absent, collapse = ","), ce_ok))
```

### Hunk 21: K275-prefix; old 558-559

Old (0866982f):
```r
558:         e12 <- if (!is.null(out$e12)) out$e12[out$e12$assay == a, , drop = FALSE]
559:         e13 <- if (!is.null(out$e13)) out$e13[out$e13$assay == a, , drop = FALSE]
```
New (sim-k275.R):
```r
582:         e12 <- if (!is.null(out$itsd_check)) out$itsd_check[out$itsd_check$assay == a, , drop = FALSE] # K275-prefix
583:         e13 <- if (!is.null(out$norm_check)) out$norm_check[out$norm_check$assay == a, , drop = FALSE] # K275-prefix
```

### Hunk 22: K275-assaycols; old 576-576

Old (0866982f):
```r
576:             if (any(res$status == "skipped_k_exceeds_complete_controls")) reason <- paste0(reason, " | ", oneLine(trace$error))
```
New (sim-k275.R):
```r
600:             if (any(res$status %in% c("skipped_k_exceeds_complete_controls", "skipped_no_controls_selected"))) reason <- paste0(reason, " | ", oneLine(trace$error)) # K275-assaycols: KIT-242's status (RG6)
```

### Hunk 23: K275-prefix; old 582-582

Old (0866982f):
```r
582:         da_status <- if (da_failed) "failed" else if (!isTRUE(out$da_reached)) "not_reached" else
```
New (sim-k275.R):
```r
606:         da_status <- if (da_failed) "failed" else if (!isTRUE(out$reached_da)) "not_reached" else # K275-prefix
```

### Hunk 24: K275-assaycols; old (insert after 587)

Old (0866982f):
```r
(nothing)
```
New (sim-k275.R):
```r
612:             ruv_not_applied_reason = if (requested && !declared) as.character(pickAssay(product[["ruv_not_applied_reason"]], a)) else NA_character_, # K275-assaycols
```

### Hunk 25: K275-assaycols; old 591-591

Old (0866982f):
```r
591:             ruv_n_controls = if (requested && !declared && !is.null(trace)) sum(trace$control_genes_index, na.rm = TRUE) else NA_integer_,
```
New (sim-k275.R):
```r
616:             ruv_n_controls = if (requested && !declared && !is.null(trace) && is.logical(trace[["control_genes_index"]])) sum(trace[["control_genes_index"]], na.rm = TRUE) else NA_integer_, # K275-assaycols
```

### Hunk 26: K275-guard; old 632-632

Old (0866982f):
```r
632:     result <- tryCatch({
```
New (sim-k275.R):
```r
657:     result <- tryCatch(withCallingHandlers({ # K275-guard
```

### Hunk 27: K275-guard; old (insert after 638)

Old (0866982f):
```r
(nothing)
```
New (sim-k275.R):
```r
664:         if (isTRUE(getOption("k275.guard_now"))) { # K275-guard: the planted partial `$` match on a data frame (the self-check of the guard)
665:             options(k275.guard_now = FALSE)
666:             out$guard_planted <- TRUE
667:             invisible(data.frame(abcdef = 1)$abc)
668:         }
```

### Hunk 28: K275-prefix; old 659-661

Old (0866982f):
```r
659:             out$e12 <- checkItsdActed(m_qc, out$m_itsd, itsd_record)
660:             if (!all(out$e12$itsd_ok)) {
661:                 bad <- out$e12[!out$e12$itsd_ok, , drop = FALSE]
```
New (sim-k275.R):
```r
689:             out$itsd_check <- checkItsdActed(m_qc, out$m_itsd, itsd_record) # K275-prefix
690:             if (!all(out$itsd_check$itsd_ok)) {
691:                 bad <- out$itsd_check[!out$itsd_check$itsd_ok, , drop = FALSE]
```

### Hunk 29: K275-pairs; old (insert after 671)

Old (0866982f):
```r
(nothing)
```
New (sim-k275.R):
```r
702:         step <- "pairs_pre" # K275-pairs
703:         out$pairs_pre <- quiet(pearsonCorForSamplePairs(theObject = owner$getState(paste0(prefix, "_normalized")), correlation_group = "group")) # K275-pairs
704:         step <- "between_sample"
```

### Hunk 30: K275-prefix; old 677-679

Old (0866982f):
```r
677:             out$e13 <- checkNormActed(compared, out$m_log2)
678:             if (!all(out$e13$norm_ok)) {
679:                 bad <- out$e13[!out$e13$norm_ok, , drop = FALSE]
```
New (sim-k275.R):
```r
710:             out$norm_check <- checkNormActed(compared, out$m_log2) # K275-prefix
711:             if (!all(out$norm_check$norm_ok)) {
712:                 bad <- out$norm_check[!out$norm_check$norm_ok, , drop = FALSE]
```

### Hunk 31: K275-W; old 691-691

Old (0866982f):
```r
691:                 "control_count", "ruv_k_estimate"))])
```
New (sim-k275.R):
```r
724:                 "control_count", "ruv_k_estimate", "W"))]) # K275-W
```

### Hunk 32: K275-pairs; old (insert after 694)

Old (0866982f):
```r
(nothing)
```
New (sim-k275.R):
```r
728:             step <- "pairs_post" # K275-pairs
729:             out$pairs_post <- quiet(pearsonCorForSamplePairs(theObject = owner$getState(paste0(prefix, "_ruv_corrected")), correlation_group = "group")) # K275-pairs
730:             step <- "ruv"
```

### Hunk 33: K275-guard; old 696-696

Old (0866982f):
```r
696:             record <- tryCatch(quiet(buildWorkflowAnalysisRecord(owner)), error = function(e) NULL)
```
New (sim-k275.R):
```r
732:             record <- tryCatch(quiet(buildWorkflowAnalysisRecord(owner)), error = function(e) if (grepl("partial-match guard", conditionMessage(e), fixed = TRUE)) stop(e) else NULL) # K275-guard
```

### Hunk 34: K275-pairs; old (insert after 703)

Old (0866982f):
```r
(nothing)
```
New (sim-k275.R):
```r
740:             # K275-pairs: per assay, the samples that held a value before the filter and none after it (the package filter's own removed set, E5)
741:             after_filter <- stateMatrices(owner, binding, owner$current_state)
742:             out$corr_removed <- lapply(stats::setNames(names(after_filter), names(after_filter)), function(a) {
743:                 before <- out$m_out[[a]] %||% out$m_in[[a]] # K275-pairs
744:                 had <- colnames(before)[colSums(is.finite(before)) > 0L]
745:                 kept <- colnames(after_filter[[a]])[colSums(is.finite(after_filter[[a]])) > 0L]
746:                 setdiff(had, kept)
747:             }) # K275-pairs
```

### Hunk 35: K275-nogroup; old (insert after 706)

Old (0866982f):
```r
(nothing)
```
New (sim-k275.R):
```r
751:         left <- unique(as.character(final@design_matrix[[final@group_id]]))
752:         if (!all(sim$groups %in% left)) { # K275-nogroup: the filter removed a whole group from the design: the DA has nothing to test (a result, `no_da`)
753:             out$class <- paste0("filter_removed_group:", paste(setdiff(sim$groups, left), collapse = ","))
754:         } else {
```

### Hunk 36: K275-prefix; old 708-708

Old (0866982f):
```r
708:         out$da_reached <- TRUE
```
New (sim-k275.R):
```r
756:         out$reached_da <- TRUE # K275-prefix
```

### Hunk 37: K275-nogroup; old (insert after 734)

Old (0866982f):
```r
(nothing)
```
New (sim-k275.R):
```r
783:         } # K275-nogroup
```

### Hunk 38: K275-guard; old 752-752

Old (0866982f):
```r
752:     }, error = function(e) {
```
New (sim-k275.R):
```r
801:     }, warning = guardWarning), error = function(e) { # K275-guard
```

### Hunk 39: K275-guard; old (insert after 753)

Old (0866982f):
```r
(nothing)
```
New (sim-k275.R):
```r
803:         guard <- grepl("partial-match guard", conditionMessage(e), fixed = TRUE) # K275-guard
```

### Hunk 40: K275-guard; old 756-756

Old (0866982f):
```r
756:         if (identical(out$status, "instrument")) {
```
New (sim-k275.R):
```r
806:         if (guard) { # K275-guard: a guard stop is an instrument defect, never a crash
807:             out$status <<- "instrument"
808:             out$class <<- oneLine(paste("partial_match_guard:", conditionMessage(e)))
809:         } else if (identical(out$status, "instrument")) {
```

### Hunk 41: K275-prefix; old 797-797

Old (0866982f):
```r
797:         k_used <- suppressWarnings(as.integer(unlist(journey$ruv_value$best_k[[a]])[1L]))
```
New (sim-k275.R):
```r
850:         k_used <- suppressWarnings(as.integer(unlist(journey[["ruv_value"]][["best_k"]][[a]])[1L])) # K275-prefix
```
