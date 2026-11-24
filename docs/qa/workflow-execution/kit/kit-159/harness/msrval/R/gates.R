# G1(a)'s re-freeze chain reader and G1(r)'s adoption of earlier reference rows (KIT-159 section 2).

# `limpa_da_last_digits` is not a default class: G1D did not admit it (its base reproduction
# failed) and D19 was ruled in VAL-109: a separate ticket (KIT-199) re-freezes the DIA-with-limpa
# reference, and until it merges G1(a) does not admit the class (audit I15).
G1A_CLASSES <- c("integration_artifacts", "kit128_technical_replicates", "metab_fixture_4",
                 "platform_rounding")
G1D_LIMPA_DIGEST_PREFIX <- "8a4c2dd5"

#' G1(a)'s classes, with `limpa_da_last_digits` only through an explicit admission
#'
#' D19 (VAL-109): the class is admitted only once KIT-199 (the re-freeze of the DIA-with-limpa
#' reference) has merged into the commit under test. The admission names that merge commit; it
#' must be an ancestor of (or equal to) `HEAD` of `pkg_root` and its subject must begin
#' `Merge KIT-199`, so a ruling's text alone admits nothing.
#'
#' @param admit `NULL`, or `list(ruling, digest, merge_commit, pkg_root)`: the ruling, the measured
#'   masked digest (beginning with G1D's `8a4c2dd5`), KIT-199's merge commit and the checkout.
#' @return The admitted classes.
#' @noRd
valG1aClasses <- function(admit = NULL) {
    if (is.null(admit)) return(G1A_CLASSES)
    refuse <- function(why) {
        valAbort(paste("limpa_da_last_digits is not admitted:", why, "(D19, VAL-109)"),
                 "val_spec_invalid")
    }
    if (!nzchar(admit[["ruling"]] %||% "")) refuse("no ruling")
    if (!startsWith(admit[["digest"]] %||% "", G1D_LIMPA_DIGEST_PREFIX)) {
        refuse("the digest is not G1D's 8a4c2dd5")
    }
    commit <- admit[["merge_commit"]] %||% ""
    root <- admit[["pkg_root"]] %||% ""
    if (!nzchar(commit) || !nzchar(root) || !dir.exists(root)) {
        refuse("KIT-199's merge commit and the checkout are not named")
    }
    git <- function(...) {
        suppressWarnings(system2("git", c("-C", shQuote(root), ...), stdout = TRUE,
                                 stderr = FALSE))
    }
    subject <- git("log", "-1", "--format=%s", shQuote(commit))
    if (!identical(attr(subject, "status"), NULL) || length(subject) != 1L ||
            !startsWith(subject, "Merge KIT-199")) {
        refuse("the named commit is not KIT-199's merge")
    }
    anc <- suppressWarnings(system2("git", c("-C", shQuote(root), "merge-base", "--is-ancestor",
                                             shQuote(commit), "HEAD"), stdout = FALSE,
                                    stderr = FALSE))
    if (!identical(as.integer(anc), 0L)) refuse("KIT-199's merge is not in the commit under test")
    c(G1A_CLASSES, "limpa_da_last_digits")
}

G1A_LIMPA_ARTIFACT <- "results/proteomics/da_proteins/da_proteins_groupKO-groupWT_long_annot.tsv"

#' The lane-range classes of G1(a): lanes and the step each class starts at
#' @noRd
G1A_LANE_RANGES <- list(
    kit128_technical_replicates = list(lanes = "trx_screen", from = "design"),
    metab_fixture_4 = list(lanes = c("metab_lc", "metab_gc", "metab_combined"), from = "norm"),
    limpa_da_last_digits = list(lanes = "prot_dia_limpa", from = "da")
)

#' Class G1(a)'s live differences by the bounds of its pass column (audit I15)
#'
#' One row per difference of a lane's live capture from its oracle: `lane`, `step`, `field`
#' (`artifacts`, `states`, `outputs`, `store_entries`, ...), `key`, `base` and `live` (digests or
#' counts; `NA` where the key is on one side only). Per step: a lane-range class takes the whole
#' step (`kit128_technical_replicates` on `trx_screen` from `design` on, `metab_fixture_4` on the three
#' metabolomics lanes from `norm` on, even where integration names also differ, C08);
#' `limpa_da_last_digits` only when admitted ([valG1aClasses()]), on `prot_dia_limpa` from `da` on,
#' for the DA artifact alone with a live digest beginning with the admitted digest and no `states`
#' or `outputs` difference in the step; `integration_artifacts` when every difference of the step is
#' an `artifacts` key `integration/<omic>/...` present on one side only or a `store_entries` count;
#' `platform_rounding` only for a key with a committed reference fingerprint whose comparison found
#' rounding only (`fingerprints`: `lane`, `step`, `key`, `rounding_only`). Anything else is `NA`:
#' outside the classes, G1(a) stops on it.
#'
#' @param diffs The differences.
#' @param lane_steps Named list: each lane's steps in order.
#' @param admit The admission of `limpa_da_last_digits`, or `NULL`.
#' @param fingerprints Committed reference fingerprints, or `NULL` (an oracle has none).
#' @return `diffs` with `class` (and `digest`, the live value, for the chain reader).
#' @export
valClassifyLiveDiffs <- function(diffs, lane_steps, admit = NULL, fingerprints = NULL) {
    classes <- valG1aClasses(admit)
    unknown <- setdiff(unique(diffs$lane), names(lane_steps))
    if (length(unknown)) {
        valAbort("lanes without a step list: {paste(unknown, collapse = ', ')}",
                 "val_spec_invalid")
    }
    cls <- rep(NA_character_, nrow(diffs))
    pos <- function(lane, step) match(step, lane_steps[[lane]])
    for (g in split(seq_len(nrow(diffs)), paste(diffs$lane, diffs$step, sep = "\r"))) {
        d <- diffs[g, , drop = FALSE]
        lane <- d$lane[[1L]]
        at <- pos(lane, d$step[[1L]])
        if (is.na(at)) next
        range_class <- NA_character_
        for (rc in names(G1A_LANE_RANGES)) {
            r <- G1A_LANE_RANGES[[rc]]
            if (rc %in% classes && lane %in% r$lanes && at >= pos(lane, r$from)) range_class <- rc
        }
        if (identical(range_class, "limpa_da_last_digits")) {
            live <- ifelse(is.na(d$live), "", d$live)
            ok <- all(d$field == "artifacts" & d$key == G1A_LIMPA_ARTIFACT) &&
                all(startsWith(live, admit[["digest"]]))
            if (ok) cls[g] <- range_class
            next
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
        if (!is.null(fingerprints) && "platform_rounding" %in% classes) {
            fp <- fingerprints[fingerprints$lane == lane & fingerprints$step == d$step[[1L]] &
                                   fingerprints$rounding_only %in% TRUE, , drop = FALSE]
            if (all(d$key %in% fp$key)) cls[g] <- "platform_rounding"
        }
    }
    diffs$class <- cls
    diffs$digest <- diffs$live
    diffs
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
