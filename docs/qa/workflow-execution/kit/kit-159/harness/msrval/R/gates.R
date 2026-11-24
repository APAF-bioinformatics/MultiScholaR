# G1(a)'s re-freeze chain reader and G1(r)'s adoption of earlier reference rows (KIT-159 section 2).

# `limpa_da_last_digits` is not a default class: G1D did not admit it (its base reproduction
# failed) and D19 was ruled in VAL-109: a separate ticket (KIT-199) re-freezes the DIA-with-limpa
# reference, and until it merges G1(a) does not admit the class (audit I15).
G1A_CLASSES <- c("integration_artifacts", "kit128_technical_replicates", "metab_fixture_4",
                 "platform_rounding")
G1D_LIMPA_DIGEST_PREFIX <- "8a4c2dd5"

#' G1(a)'s classes, with `limpa_da_last_digits` only through an explicit admission
#'
#' @param admit `NULL`, or `list(ruling, digest)`: the ruling that admits the class (KIT-199's
#'   merge) and the measured masked digest, which must begin with G1D's `8a4c2dd5`.
#' @return The admitted classes.
#' @noRd
valG1aClasses <- function(admit = NULL) {
    if (is.null(admit)) return(G1A_CLASSES)
    if (!nzchar(admit[["ruling"]] %||% "") ||
            !startsWith(admit[["digest"]] %||% "", G1D_LIMPA_DIGEST_PREFIX)) {
        valAbort("limpa_da_last_digits is admitted only with its ruling and G1D's digest",
                 "val_spec_invalid")
    }
    c(G1A_CLASSES, "limpa_da_last_digits")
}

#' Every difference of the admitted limpa class carries G1D's digest
#' @noRd
valLimpaDigestsOk <- function(diffs, admit) {
    i <- diffs$class %in% "limpa_da_last_digits"
    !any(i) || (!is.null(diffs$digest) && all(diffs$digest[i] %in% admit[["digest"]]))
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
#' @param admit `NULL`, or the admission of `limpa_da_last_digits` ([valG1aClasses()]).
#' @return A list: `status` (`certified`, `not_certified`, `stop`), `failed_link`, `reason`.
#' @export
valReadChain <- function(lane, links, reference, live_diffs = NULL, admit = NULL) {
    classes <- valG1aClasses(admit)
    ref_counts <- valClassCounts(reference$class[reference$lane == lane])
    if (length(links) == 0L) {
        ok <- !is.null(live_diffs) && all(live_diffs$class %in% classes) &&
            valLimpaDigestsOk(live_diffs, admit) &&
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
        ok_base <- if (i == 1L) {
            all(bd$class %in% classes) && valLimpaDigestsOk(bd, admit) &&
                valSameCounts(valClassCounts(bd$class), ref_counts)
        } else {
            all(bd$class %in% "platform_rounding")
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
