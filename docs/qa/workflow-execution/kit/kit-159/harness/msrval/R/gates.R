# G1(a)'s re-freeze chain reader, its live-difference classifier and G1(r)'s adoption of earlier
# reference rows (KIT-159 section 2).

# `limpa_da_last_digits` is not a class of G1(a). D19 was answered by VAL-109 (5) and the KIT-199
# plan audit's ruling C1 (a): KIT-199 re-freezes the DIA-with-limpa reference, and its link is read
# as `declared_change (KIT-199)` like any re-freezing ticket's; the harness keeps no separate
# admission of the class (the stored value is a digest of masked text, so the kind and size of the
# laptop-runner difference are not known).
G1A_CLASSES <- c("integration_artifacts", "kit128_technical_replicates", "metab_fixture_4",
                 "platform_rounding")

#' The class a re-freezing ticket's record gives a difference it declared
#' @noRd
valDeclaredClass <- function(ticket) sprintf("declared_change (%s)", ticket)

#' The lane-range classes of G1(a): lanes and the step each class starts at
#' @noRd
G1A_LANE_RANGES <- list(
    kit128_technical_replicates = list(lanes = "trx_screen", from = "design"),
    metab_fixture_4 = list(lanes = c("metab_lc", "metab_gc", "metab_combined"), from = "norm")
)

#' Class G1(a)'s live differences by the bounds of its pass column (audit I15)
#'
#' One row per difference of a lane's live capture from its oracle: `lane`, `step`, `field`
#' (`artifacts`, `states`, `outputs`, `store_entries`, ...), `key`, `base` and `live` (digests or
#' counts; `NA` where the key is on one side only). A difference a re-freezing ticket declared
#' (`declared`: `lane`, `step`, `key`, `ticket`) is `declared_change (<ticket>)` (the KIT-199 plan
#' audit's ruling C1 (a)); the step's other differences are classed among themselves, so a step may
#' carry a declared change and a class. A lane-range class takes the rest of the step
#' (`kit128_technical_replicates` on `trx_screen` from `design` on, `metab_fixture_4` on the three
#' metabolomics lanes from `norm` on, even where integration names also differ, C08);
#' `integration_artifacts` when every other difference of the step is an `artifacts` key
#' `integration/<omic>/...` present on one side only or a `store_entries` count;
#' `platform_rounding` only for a key with a committed reference fingerprint whose comparison found
#' rounding only (`fingerprints`: `lane`, `step`, `key`, `rounding_only`). Anything else is `NA`:
#' outside the classes, G1(a) stops on it. There is no `limpa_da_last_digits` class.
#'
#' @param diffs The differences.
#' @param lane_steps Named list: each lane's steps in order.
#' @param declared Declared changes of re-freezing tickets, or `NULL`.
#' @param fingerprints Committed reference fingerprints, or `NULL` (an oracle has none).
#' @return `diffs` with `class` (and `digest`, the live value, for the chain reader).
#' @export
valClassifyLiveDiffs <- function(diffs, lane_steps, declared = NULL, fingerprints = NULL) {
    unknown <- setdiff(unique(diffs$lane), names(lane_steps))
    if (length(unknown)) {
        valAbort("lanes without a step list: {paste(unknown, collapse = ', ')}",
                 "val_spec_invalid")
    }
    cls <- rep(NA_character_, nrow(diffs))
    if (!is.null(declared) && nrow(declared)) {
        k <- match(paste(diffs$lane, diffs$step, diffs$key, sep = "\r"),
                   paste(declared$lane, declared$step, declared$key, sep = "\r"))
        cls[!is.na(k)] <- valDeclaredClass(declared$ticket[k[!is.na(k)]])
    }
    pos <- function(lane, step) match(step, lane_steps[[lane]])
    open <- which(is.na(cls))
    for (g in split(open, paste(diffs$lane[open], diffs$step[open], sep = "\r"))) {
        d <- diffs[g, , drop = FALSE]
        lane <- d$lane[[1L]]
        at <- pos(lane, d$step[[1L]])
        if (is.na(at)) next
        range_class <- NA_character_
        for (rc in names(G1A_LANE_RANGES)) {
            r <- G1A_LANE_RANGES[[rc]]
            if (lane %in% r$lanes && at >= pos(lane, r$from)) range_class <- rc
        }
        if (!is.na(range_class)) {
            cls[g] <- range_class
            next
        }
        integ <- (d$field == "artifacts" & startsWith(d$key, "integration/") &
                      (is.na(d$base) | is.na(d$live))) | d$field == "store_entries"
        if (all(integ)) {
            cls[g] <- "integration_artifacts"
            next
        }
        if (!is.null(fingerprints)) {
            fp <- fingerprints[fingerprints$lane == lane & fingerprints$step == d$step[[1L]] &
                                   fingerprints$rounding_only %in% TRUE, , drop = FALSE]
            if (all(d$key %in% fp$key)) cls[g] <- "platform_rounding"
        }
    }
    diffs$class <- cls
    diffs$digest <- diffs$live
    diffs
}

#' Counts per class
#' @noRd
valClassCounts <- function(cls) {
    if (length(cls) == 0L) return(stats::setNames(integer(), character()))
    t <- table(cls)
    stats::setNames(as.integer(t), names(t))
}

#' Counts per class of (step, class) units: G1(a)'s reference table has one row per lane and
#' differing step with its class, so a step with several differing keys of one class is one unit
#' (KIT-199 re-check 1, Q01); rows without a step count one each
#' @noRd
valStepClassCounts <- function(d) {
    if (is.null(d) || nrow(d) == 0L) return(valClassCounts(character()))
    if (is.null(d$step)) return(valClassCounts(d$class))
    u <- unique(data.frame(step = d$step, class = d$class, stringsAsFactors = FALSE))
    valClassCounts(u$class)
}

#' Equal counts per class
#' @noRd
valSameCounts <- function(a, b) {
    identical(sort(as.character(names(a))), sort(as.character(names(b)))) &&
        all(a[names(b)] == b)
}

#' Read one lane's re-freeze chain (G1(a), round-2 C01; r1.5 N1, N2)
#'
#' With no re-freeze since 12d68ca7 the live differences must fall in G1(a)'s classes with the
#' reference table's counts (else `stop`). Each link is a ticket's control record: (i) the base's
#' differences against the file at that base: for the first link in G1(a)'s classes with the
#' reference counts, for later links none or `platform_rounding` only; (ii) before the changed step
#' only keys classed `declared_change (<ticket>)`, each listed with its expected value and equal to
#' it; (iii) from the changed step on every difference classed by the ticket's own change. A missing
#' record, or a link failing (i), (ii) or (iii), makes the lane `not_certified`.
#'
#' @param lane Lane id.
#' @param links Control records in merge order (`NULL` for a re-frozen file without one): lists
#'   with `ticket`, `changed_step`, `base_diffs`, `before_diffs`, `after_diffs`, `declared`,
#'   `change_classes`.
#' @param reference `g1a-reference.tsv` rows (lane, step, class, fields, keys).
#' @param live_diffs For a lane with no link: the live run's differences from the oracle.
#' @return A list: `status` (`certified`, `not_certified`, `stop`), `failed_link`, `reason`.
#' @export
valReadChain <- function(lane, links, reference, live_diffs = NULL) {
    classes <- G1A_CLASSES
    ref_counts <- valStepClassCounts(reference[reference$lane == lane, , drop = FALSE])
    if (length(links) == 0L) {
        ok <- !is.null(live_diffs) && all(live_diffs$class %in% classes) &&
            valSameCounts(valStepClassCounts(live_diffs), ref_counts)
        return(list(status = if (ok) "certified" else "stop", failed_link = NA_integer_,
                    reason = if (ok) "" else "differences outside the classes or other counts"))
    }
    fail <- function(i, why) {
        list(status = "not_certified", failed_link = as.integer(i), reason = why)
    }
    for (i in seq_along(links)) {
        link <- links[[i]]
        if (is.null(link)) return(fail(i, "a re-frozen oracle without its control record"))
        bd <- link$base_diffs
        # The KIT-199 plan audit's ruling C1 (a): link (i) accepts a base difference the link's own
        # record classes `declared_change (<ticket>)` and lists (step and key) in `declared`; the
        # reference counts compare the rest, so a step may carry a declared change and a class.
        own <- bd$class %in% valDeclaredClass(link$ticket)
        listed <- paste(bd$step, bd$key, sep = "\r") %in%
            paste(link$declared$step, link$declared$key, sep = "\r")
        rest <- bd[!own, , drop = FALSE]
        base_counts <- ref_counts[!startsWith(names(ref_counts), "declared_change")]
        # KIT-199 12.3: each declared row's stored and live values equal `declared-keys.tsv`'s
        # `frozen` and `expected` (the KIT-199 implementation audit's owed item)
        valued <- valDeclaredValuesOk(bd[own, , drop = FALSE], link$declared)
        ok_base <- valued && all(listed[own]) && if (i == 1L) {
            all(rest$class %in% classes) &&
                valSameCounts(valStepClassCounts(rest), base_counts)
        } else {
            all(rest$class %in% "platform_rounding")
        }
        if (!ok_base) return(fail(i, "(i) the base's differences"))
        bf <- link$before_diffs
        declared_class <- sprintf("declared_change (%s)", link$ticket)
        ok_before <- all(vapply(seq_len(nrow(bf)), function(r) {
            d <- link$declared[link$declared$key == bf$key[[r]] &
                                   link$declared$step == bf$step[[r]], , drop = FALSE]
            identical(bf$class[[r]], declared_class) && nrow(d) == 1L &&
                identical(as.character(d$expected), as.character(bf$head_value[[r]]))
        }, logical(1)))
        if (!ok_before) return(fail(i, "(ii) a difference before the changed step"))
        af <- link$after_diffs
        if (!all(!is.na(af$class) & af$class %in% link$change_classes)) {
            return(fail(i, "(iii) a difference without a class of the ticket's change"))
        }
    }
    list(status = "certified", failed_link = NA_integer_, reason = "")
}

#' Link (i)'s value check of a link's own declared rows (KIT-199 12.3)
#' @noRd
valDeclaredValuesOk <- function(rows, declared) {
    if (nrow(rows) == 0L) return(TRUE)
    if (is.null(rows$base) || is.null(rows$live) || is.null(declared$frozen) ||
            is.null(declared$expected)) {
        return(FALSE)
    }
    k <- match(paste(rows$step, rows$key, sep = "\r"),
               paste(declared$step, declared$key, sep = "\r"))
    !anyNA(k) && !anyNA(rows$base) &&
        identical(as.character(rows$base), as.character(declared$frozen[k])) &&
        identical(as.character(rows$live), as.character(declared$expected[k]))
}

#' Adopt reference rows a re-freezing ticket committed earlier, when they agree with G1(r)'s own
#' rows for its lanes (r1.5; N2 of the KIT-172 re-check; KIT-199 12.3)
#'
#' Compared by (lane, step, class) and by key sets (`field/key`), in either committed column form
#' (`fields`/`keys` or `differing_fields_and_keys`), never by the free text or the column names. A
#' row whose source listed only part of the names (`names_list_truncated`, KIT-172) must have its
#' keys among G1(r)'s own.
#'
#' A committed key erratum (`inst/extdata/g1r-key-errata.tsv`; coordinator's ruling on re-audit 3,
#' R3.B3) replaces an earlier row's key by G1(r)'s when lane, step, class and the earlier key all
#' match a listed pair exactly; nothing else is corrected.
#'
#' @param own G1(r)'s own rows (lane, step, class, fields, keys).
#' @param earlier The ticket's committed rows (read by [valReadReferenceRows()], or raw).
#' @param errata Key errata (`lane`, `step`, `class`, `earlier_key`, `key`), or `NULL`.
#' @return `own`, or `val_spec_invalid` on any disagreement (neither file is edited to agree).
#' @export
valAdoptReferenceRows <- function(own, earlier, errata = valG1rKeyErrata()) {
    if (nrow(valReferenceDisagreements(own, earlier, errata))) {
        valAbort("a ticket's committed reference rows disagree with G1(r)'s own rows",
                 "val_spec_invalid")
    }
    own
}

#' The (lane, step, class) groups on which a ticket's committed reference rows and G1(r)'s own
#' rows for its lanes disagree (the rule of [valAdoptReferenceRows()], for the gate's report)
#'
#' @param own,earlier,errata As [valAdoptReferenceRows()].
#' @return A data frame, one row per disagreeing group: `lane`, `step`, `class`, `in_own`,
#'   `in_earlier`, `own_only` and `earlier_only` (`field/key`, semicolon-separated), `truncated`.
#' @export
valReferenceDisagreements <- function(own, earlier, errata = valG1rKeyErrata()) {
    if (!is.null(errata) && nrow(errata)) {
        if (is.null(earlier$keyset)) {
            earlier$keyset <- lapply(seq_len(nrow(earlier)), function(i) {
                valReferenceKeys(earlier[i, , drop = FALSE])
            })
        }
        for (i in seq_len(nrow(earlier))) {
            e <- errata[errata$lane == earlier$lane[[i]] & errata$step == earlier$step[[i]] &
                            errata$class == earlier$class[[i]], , drop = FALSE]
            k <- earlier$keyset[[i]]
            hit <- match(k, e$earlier_key)
            k[!is.na(hit)] <- e$key[hit[!is.na(hit)]]
            earlier$keyset[[i]] <- k
        }
    }
    sets <- function(x) {
        if (is.null(x$keyset)) {
            x$keyset <- lapply(seq_len(nrow(x)), function(i) valReferenceKeys(x[i, , drop = FALSE]))
        }
        if (is.null(x$truncated)) {
            x$truncated <- (x$names_list_truncated %||% rep("", nrow(x))) %in% "yes"
        }
        g <- paste(x$lane, x$step, x$class, sep = "\r")
        list(keys = lapply(split(x$keyset, g), function(k) sort(unique(unlist(k)))),
             truncated = vapply(split(x$truncated, g), any, logical(1)))
    }
    a <- sets(own[own$lane %in% earlier$lane, , drop = FALSE])
    b <- sets(earlier)
    groups <- sort(union(names(a$keys), names(b$keys)))
    rows <- lapply(groups, function(g) {
        ka <- a$keys[[g]]
        kb <- b$keys[[g]]
        truncated <- isTRUE(b$truncated[g])
        agree <- !is.null(ka) && !is.null(kb) &&
            (if (truncated) all(kb %in% ka) else identical(ka, kb))
        if (agree) return(NULL)
        p <- strsplit(g, "\r", fixed = TRUE)[[1L]]
        data.frame(lane = p[[1L]], step = p[[2L]], class = p[[3L]], in_own = !is.null(ka),
                   in_earlier = !is.null(kb), own_only = paste(setdiff(ka, kb), collapse = ";"),
                   earlier_only = paste(setdiff(kb, ka), collapse = ";"), truncated = truncated,
                   stringsAsFactors = FALSE)
    })
    empty <- data.frame(lane = character(), step = character(), class = character(),
                        in_own = logical(), in_earlier = logical(), own_only = character(),
                        earlier_only = character(), truncated = logical(),
                        stringsAsFactors = FALSE)
    do.call(rbind, c(list(empty), rows))
}

#' G1(r)'s committed key errata: KIT-172's 18 keys that kept waldo's line window, each with the
#' output id it names (coordinator's ruling on re-audit 3, R3.B3; KIT-172's file is not edited)
#' @noRd
valG1rKeyErrata <- function() {
    valReadTextTable(system.file("extdata", "g1r-key-errata.tsv", package = "msrval",
                                 mustWork = TRUE))
}
