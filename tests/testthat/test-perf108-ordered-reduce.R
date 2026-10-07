test_that("ordered reduction preserves exact sums across transport boundaries", {
    skip_if(Sys.getenv("MULTISCHOLAR_PERF108_STRUCTURAL") != "S")
    reduce <- function(raw, normalised, forms = rep("m", length(raw)), sample = rep("a", length(raw))) {
        state <- .Call(C_workflow_dia_ordered_new)
        on.exit(.Call(C_workflow_dia_ordered_close, state))
        input <- data.frame(sample = sample, protein = rep("P", length(raw)), peptide = rep("SEQ", length(raw)),
            raw = as.double(raw), normalised = as.double(normalised), form = forms,
            source_order = as.double(seq_along(raw)))
        pieces <- list()
        if (nrow(input)) for (start in seq.int(1L, nrow(input), by = 10000L)) {
            batch <- input[start:min(nrow(input), start + 9999L), ]
            pieces[[length(pieces) + 1L]] <- .Call(C_workflow_dia_ordered_push, state, batch, FALSE)
            batch$sample[] <- "mutated after delivery"
            batch$raw[] <- -1
            invisible(gc())
        }
        final <- .Call(C_workflow_dia_ordered_push, state, input[0, ], TRUE)
        statistics <- attr(final, "ordered_reducer_statistics")
        attr(final, "ordered_reducer_statistics") <- NULL
        pieces[[length(pieces) + 1L]] <- final
        list(data = do.call(rbind, pieces), statistics = statistics)
    }
    values <- list(c(1e16, 1, 1), c(0, NA_real_, 1), rep(NA_real_, 3),
        c(.Machine$double.xmin / 2, .Machine$double.xmin / 2, 0),
        c(.Machine$double.xmax, .Machine$double.xmax), rep(c(1e16, 1, 1), 10001L))
    expected <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
    for (raw in values) {
        norm <- rev(raw)
        forms <- ifelse(seq_along(raw) %% 11 == 0, NA_character_, paste0("m", seq_along(raw) %% 211))
        result <- reduce(raw, norm, forms)
        expect_identical(result$data$Peptide.RawQuantity, expected(raw))
        expect_identical(result$data$Peptide.Normalised, expected(norm))
        expect_identical(result$data$peptidoform_count, length(unique(forms[!is.na(forms)])))
        expect_identical(result$data$source_count, length(raw))
        expect_equal(result$statistics[["input_rows"]], length(raw))
        expect_equal(result$statistics[["max_active_groups"]], 1)
        expect_lte(result$statistics[["max_group_state_bytes"]], 67108864)
    }
    raw <- c(rep(1, 10001L), rep(NA_real_, 20003L), 0, 1e16, 1, 1)
    norm <- c(rep(NA_real_, 10001L), rep(c(1e16, 1, 1), length.out = 20003L), rep(0, 4L))
    sample <- rep(c("a", "b", "\u03b2"), c(10001L, 20003L, 4L))
    result <- reduce(raw, norm, sample = sample)
    expect_identical(result$data$sample, c("a", "b", "\u03b2"))
    expect_identical(result$data$Peptide.RawQuantity,
        unname(vapply(split(raw, sample), expected, numeric(1))))
    expect_identical(result$data$Peptide.Normalised,
        unname(vapply(split(norm, sample), expected, numeric(1))))
    expect_equal(result$statistics[["groups"]], 3)
    empty <- reduce(double(), double())
    expect_equal(nrow(empty$data), 0L)
    expect_equal(empty$statistics[["max_active_groups"]], 0)
})

test_that("ordered reducers fail closed on invalid values or ordering", {
    skip_if(Sys.getenv("MULTISCHOLAR_PERF108_STRUCTURAL") != "S")
    frame <- data.frame(sample = "b", protein = "P", peptide = "SEQ", raw = 1,
        normalised = 1, form = "m", source_order = 1)
    for (kind in c("negative", "infinite", "nan", "key_order", "ordinal", "missing_key")) {
        state <- .Call(C_workflow_dia_ordered_new)
        .Call(C_workflow_dia_ordered_push, state, frame, FALSE)
        bad <- frame
        bad$source_order <- 2
        if (kind %in% c("negative", "infinite", "nan")) bad$raw <- switch(kind, negative = -1, infinite = Inf, nan = NaN)
        if (kind == "key_order") bad$sample <- "a"
        if (kind == "ordinal") bad$source_order <- 1
        if (kind == "missing_key") bad$sample <- NA_character_
        expect_error(.Call(C_workflow_dia_ordered_push, state, bad, FALSE),
            "nonnegative|ascending|strictly increasing|missing ordered")
        expect_error(.Call(C_workflow_dia_ordered_push, state, frame[0, ], TRUE), "failed")
        expect_no_error(.Call(C_workflow_dia_ordered_close, state))
        expect_no_error(.Call(C_workflow_dia_ordered_close, state))
        expect_error(.Call(C_workflow_dia_ordered_push, state, frame, FALSE), "closed")
    }
})
