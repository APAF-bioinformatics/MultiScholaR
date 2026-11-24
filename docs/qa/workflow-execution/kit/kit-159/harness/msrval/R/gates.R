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
    ref_counts <- valClassCounts(reference$class[reference$lane == lane])
    if (length(links) == 0L) {
        ok <- !is.null(live_diffs) && all(live_diffs$class %in% classes) &&
            valSameCounts(valClassCounts(live_diffs$class), ref_counts)
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
        ok_base <- all(listed[own]) && if (i == 1L) {
            all(rest$class %in% classes) &&
                valSameCounts(valClassCounts(rest$class), base_counts)
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

#' Adopt reference rows a re-freezing ticket committed earlier, when they agree with G1(r)'s own
#' rows for its lanes (r1.5; N2 of the KIT-172 re-check)
#'
#' @param own G1(r)'s own rows (lane, step, class, fields, keys).
#' @param earlier The ticket's committed rows.
#' @return `own`, or `val_spec_invalid` on any disagreement (neither file is edited to agree).
#' @export
valAdoptReferenceRows <- function(own, earlier) {
    cols <- c("lane", "step", "class", "fields", "keys")
    key <- function(x) sort(do.call(paste, c(x[cols], sep = "\r")))
    mine <- own[own$lane %in% earlier$lane, , drop = FALSE]
    if (!identical(key(mine), key(earlier))) {
        valAbort("a ticket's committed reference rows disagree with G1(r)'s own rows",
                 "val_spec_invalid")
    }
    own
}
