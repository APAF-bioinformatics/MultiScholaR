publicationScheduleAbort <- function(message) {
    publicationAbort(message, "multischolar_publication_schedule_error")
}

publicationScheduleDigest <- function(schedule) {
    candidate <- schedule
    candidate$generated_at <- NULL
    candidate$schedule_digest <- NULL
    publicationObjectDigest(candidate)
}

publicationScheduleProjectRecords <- function(projects, subject_id) {
    records <- projects[[subject_id]]
    if (!is.list(records) || length(records) < 3L) {
        publicationScheduleAbort(paste(
            "Subject requires at least three project records:",
            subject_id
        ))
    }
    required <- c("project_id", "source_id")
    for (record in records) {
        if (!is.list(record) || !all(required %in% names(record)) ||
            !all(vapply(record[required], publicationScalarString, logical(1)))) {
            publicationScheduleAbort("Project record is malformed")
        }
    }
    ids <- vapply(records, `[[`, character(1), "project_id")
    sources <- vapply(records, `[[`, character(1), "source_id")
    if (anyDuplicated(ids) || anyDuplicated(sources)) {
        publicationScheduleAbort("Project and source ids must be independent")
    }
    records
}

publicationSchedulePairAssignments <- function(
    project_records,
    pair_count,
    session_count,
    seed
) {
    project_count <- length(project_records)
    if (pair_count < 30L || pair_count > 60L || pair_count %% project_count != 0L ||
        (pair_count / project_count) %% 2L != 0L || session_count < 3L) {
        publicationScheduleAbort(
            "Pair count must be 30-60, project-balanced, and arm-balanced"
        )
    }
    assignments <- do.call(rbind, lapply(seq_along(project_records), \(index) {
        count <- pair_count / project_count
        data.frame(
            project_index = index,
            within_project_pair = seq_len(count),
            first_arm = rep(c("A", "B"), length.out = count),
            stringsAsFactors = FALSE
        )
    }))
    assignments <- publicationWithPreservedSeed(
        seed,
        assignments[sample.int(nrow(assignments)), , drop = FALSE]
    )
    assignments$session_index <- rep(seq_len(session_count), length.out = nrow(assignments))
    assignments$pair_index <- seq_len(nrow(assignments))
    assignments
}

publicationScheduleWorkKey <- function(subject_id, project_id, estimand_id) {
    paste(subject_id, project_id, estimand_id, sep = "::")
}

publicationScheduleWorkBinding <- function(
    work_bindings,
    subject_id,
    project_id,
    estimand_id
) {
    key <- publicationScheduleWorkKey(subject_id, project_id, estimand_id)
    binding <- work_bindings[[key]]
    expected <- c("work_unit_id", "work_count", "exact_input_bytes", "input_sha256")
    if (!is.list(binding) || !setequal(names(binding), expected) ||
        !publicationScalarString(binding$work_unit_id) ||
        !publicationScalarNumber(binding$work_count, positive = TRUE) ||
        !publicationScalarNumber(binding$exact_input_bytes, positive = TRUE) ||
        !publicationScalarString(binding$input_sha256) ||
        !grepl("^[0-9a-f]{64}$", binding$input_sha256)) {
        publicationScheduleAbort(paste("Work binding is invalid:", key))
    }
    binding
}

publicationScheduleWorkDigest <- function(slots) {
    slot_keys <- vapply(slots, \(slot) publicationScheduleWorkKey(
        slot$subject_id,
        slot$project_id,
        slot$estimand_id
    ), character(1))
    keys <- sort(unique(slot_keys), method = "radix")
    records <- stats::setNames(lapply(keys, \(key) {
        slot <- slots[[which(slot_keys == key)[[1L]]]]
        slot[c("work_unit_id", "work_count", "exact_input_bytes", "input_sha256")]
    }), keys)
    publicationObjectDigest(records)
}

publicationScheduleWarmupSlots <- function(
    subject_id,
    project_records,
    estimand_id,
    comparison_id,
    host_id,
    work_bindings,
    warmups_per_arm,
    ordinal_start
) {
    slots <- list()
    ordinal <- ordinal_start
    for (project in project_records) {
        work <- publicationScheduleWorkBinding(
            work_bindings,
            subject_id,
            project$project_id,
            estimand_id
        )
        for (arm in c("A", "B")) {
            for (warmup_index in seq_len(warmups_per_arm)) {
                slots[[length(slots) + 1L]] <- list(
                    slot_id = paste(
                        "warmup", subject_id, project$project_id, estimand_id,
                        comparison_id, arm, warmup_index,
                        sep = "::"
                    ),
                    ordinal = ordinal,
                    subject_id = subject_id,
                    project_id = project$project_id,
                    source_id = project$source_id,
                    estimand_id = estimand_id,
                    comparison_id = comparison_id,
                    pair_id = NULL,
                    block_id = paste(
                        "warmup-block", subject_id, project$project_id,
                        estimand_id, comparison_id, arm, warmup_index,
                        sep = "::"
                    ),
                    session_id = "warmup",
                    host_id = host_id,
                    arm = arm,
                    sequence_in_pair = NULL,
                    warmup = TRUE,
                    promotion_authority = FALSE,
                    cache_stratum = "standardized_warm_input",
                    work_count = work$work_count,
                    work_unit_id = work$work_unit_id,
                    exact_input_bytes = work$exact_input_bytes,
                    input_sha256 = work$input_sha256
                )
                ordinal <- ordinal + 1L
            }
        }
    }
    list(slots = slots, next_ordinal = ordinal)
}

publicationScheduleMeasuredSlots <- function(
    subject_id,
    project_records,
    estimand_id,
    comparison_id,
    host_id,
    work_bindings,
    pair_count,
    session_count,
    seed,
    ordinal_start
) {
    assignments <- publicationSchedulePairAssignments(
        project_records,
        pair_count,
        session_count,
        seed
    )
    slots <- list()
    ordinal <- ordinal_start
    for (row_index in seq_len(nrow(assignments))) {
        assignment <- assignments[row_index, ]
        project <- project_records[[assignment$project_index]]
        work <- publicationScheduleWorkBinding(
            work_bindings,
            subject_id,
            project$project_id,
            estimand_id
        )
        pair_id <- paste(
            subject_id,
            project$project_id,
            estimand_id,
            comparison_id,
            sprintf("pair-%03d", assignment$pair_index),
            sep = "::"
        )
        arms <- if (identical(assignment$first_arm, "A")) {
            c("A", "B")
        } else {
            c("B", "A")
        }
        for (sequence_index in seq_along(arms)) {
            slots[[length(slots) + 1L]] <- list(
                slot_id = paste(pair_id, arms[[sequence_index]], sep = "::"),
                ordinal = ordinal,
                subject_id = subject_id,
                project_id = project$project_id,
                source_id = project$source_id,
                estimand_id = estimand_id,
                comparison_id = comparison_id,
                pair_id = pair_id,
                block_id = pair_id,
                session_id = paste0("session-", assignment$session_index),
                host_id = host_id,
                arm = arms[[sequence_index]],
                sequence_in_pair = sequence_index,
                warmup = FALSE,
                promotion_authority = TRUE,
                cache_stratum = "standardized_warm_input",
                work_count = work$work_count,
                work_unit_id = work$work_unit_id,
                exact_input_bytes = work$exact_input_bytes,
                input_sha256 = work$input_sha256
            )
            ordinal <- ordinal + 1L
        }
    }
    list(slots = slots, next_ordinal = ordinal)
}

publicationBuildSchedule <- function(
    projects,
    estimand_ids,
    comparison_ids,
    pair_count,
    host_id,
    work_bindings,
    session_count = 3L,
    warmups_per_arm = 2L,
    seed = 106229L,
    schedule_id = "publication-schedule-v1"
) {
    if (!is.list(projects) || is.null(names(projects)) || any(!nzchar(names(projects)))) {
        publicationScheduleAbort("Projects must be a named subject list")
    }
    if (!length(estimand_ids) || !length(comparison_ids) ||
        anyDuplicated(estimand_ids) || anyDuplicated(comparison_ids) ||
        !publicationScalarString(host_id) || !publicationScalarString(schedule_id)) {
        publicationScheduleAbort("Schedule dimensions are invalid")
    }
    slots <- list()
    ordinal <- 1L
    combination_index <- 0L
    for (subject_id in names(projects)) {
        records <- publicationScheduleProjectRecords(projects, subject_id)
        for (estimand_id in estimand_ids) {
            for (comparison_id in comparison_ids) {
                combination_index <- combination_index + 1L
                warmups <- publicationScheduleWarmupSlots(
                    subject_id,
                    records,
                    estimand_id,
                    comparison_id,
                    host_id,
                    work_bindings,
                    warmups_per_arm,
                    ordinal
                )
                slots <- c(slots, warmups$slots)
                ordinal <- warmups$next_ordinal
                measured <- publicationScheduleMeasuredSlots(
                    subject_id,
                    records,
                    estimand_id,
                    comparison_id,
                    host_id,
                    work_bindings,
                    pair_count,
                    session_count,
                    seed + combination_index,
                    ordinal
                )
                slots <- c(slots, measured$slots)
                ordinal <- measured$next_ordinal
            }
        }
    }
    schedule <- list(
        schema = "multischolar.omics_publication_schedule",
        schema_version = "1.0.0",
        schedule_id = schedule_id,
        generated_at = publicationUtcNow(),
        seed = as.integer(seed),
        pair_count = as.integer(pair_count),
        session_count = as.integer(session_count),
        warmups_per_arm_project_estimand = as.integer(warmups_per_arm),
        host_id = host_id,
        arm_mapping = "opaque_external_receipt",
        single_concurrency = TRUE,
        optional_stopping = FALSE,
        work_binding_digest = publicationScheduleWorkDigest(slots),
        slots = slots,
        schedule_digest = NULL
    )
    schedule$schedule_digest <- publicationScheduleDigest(schedule)
    publicationValidateSchedule(schedule)
    schedule
}

publicationValidateSchedule <- function(schedule) {
    expected <- c(
        "schema", "schema_version", "schedule_id", "generated_at", "seed",
        "pair_count", "session_count", "warmups_per_arm_project_estimand",
        "host_id", "arm_mapping", "single_concurrency", "optional_stopping",
        "work_binding_digest", "slots", "schedule_digest"
    )
    publicationRequireNames(schedule, expected, "Publication schedule")
    if (!identical(schedule$schema, "multischolar.omics_publication_schedule") ||
        !identical(schedule$schema_version, "1.0.0") ||
        !identical(schedule$arm_mapping, "opaque_external_receipt") ||
        !isTRUE(schedule$single_concurrency) ||
        isTRUE(schedule$optional_stopping) ||
        !identical(
            schedule$work_binding_digest,
            publicationScheduleWorkDigest(schedule$slots)
        )) {
        publicationScheduleAbort("Schedule header or execution policy differs")
    }
    slot_ids <- vapply(schedule$slots, `[[`, character(1), "slot_id")
    ordinals <- vapply(schedule$slots, `[[`, integer(1), "ordinal")
    if (anyDuplicated(slot_ids) || !identical(ordinals, seq_along(ordinals))) {
        publicationScheduleAbort("Schedule slots are duplicate or unordered")
    }
    slot_fields <- c(
        "slot_id", "ordinal", "subject_id", "project_id", "source_id",
        "estimand_id", "comparison_id", "pair_id", "block_id", "session_id",
        "host_id", "arm", "sequence_in_pair", "warmup", "promotion_authority",
        "cache_stratum", "work_count", "work_unit_id", "exact_input_bytes",
        "input_sha256"
    )
    for (slot in schedule$slots) {
        identity_fields <- c(
            "slot_id", "subject_id", "project_id", "source_id", "estimand_id",
            "comparison_id", "block_id", "session_id", "host_id", "arm",
            "cache_stratum", "work_unit_id", "input_sha256"
        )
        if (!setequal(names(slot), slot_fields) ||
            any(!vapply(slot[identity_fields], publicationScalarString, logical(1))) ||
            !slot$arm %in% c("A", "B") ||
            !identical(slot$cache_stratum, "standardized_warm_input") ||
            !publicationScalarNumber(slot$work_count, positive = TRUE) ||
            !publicationScalarNumber(slot$exact_input_bytes, positive = TRUE) ||
            !publicationScalarString(slot$work_unit_id) ||
            !publicationScalarString(slot$input_sha256) ||
            !grepl("^[0-9a-f]{64}$", slot$input_sha256)) {
            publicationScheduleAbort("Schedule slot work binding is invalid")
        }
    }
    publicationValidateWarmups(schedule)
    publicationValidateMeasuredPairs(schedule)
    expected_digest <- schedule$schedule_digest
    if (!identical(publicationScheduleDigest(schedule), expected_digest)) {
        publicationScheduleAbort("Schedule digest mismatch")
    }
    invisible(schedule)
}

publicationValidateWarmups <- function(schedule) {
    warmups <- Filter(\(slot) isTRUE(slot$warmup), schedule$slots)
    measured <- Filter(\(slot) !isTRUE(slot$warmup), schedule$slots)
    groups <- split(warmups, vapply(warmups, \(slot) paste(
        slot$subject_id,
        slot$project_id,
        slot$estimand_id,
        slot$comparison_id,
        slot$arm,
        sep = "\x1f"
    ), character(1)))
    expected_groups <- length(unique(vapply(measured, \(slot) paste(
        slot$subject_id,
        slot$project_id,
        slot$estimand_id,
        slot$comparison_id,
        sep = "\x1f"
    ), character(1)))) * 2L
    valid <- length(groups) == expected_groups && all(vapply(groups, \(slots) {
        length(slots) == schedule$warmups_per_arm_project_estimand &&
            all(vapply(slots, \(slot) {
                is.null(slot$pair_id) && is.null(slot$sequence_in_pair) &&
                    identical(slot$session_id, "warmup") &&
                    !isTRUE(slot$promotion_authority)
            }, logical(1)))
    }, logical(1)))
    if (!valid) publicationScheduleAbort("Warm-up schedule differs")
    invisible(TRUE)
}

publicationValidateMeasuredPairs <- function(schedule) {
    measured <- Filter(\(slot) !isTRUE(slot$warmup), schedule$slots)
    pair_ids <- unique(vapply(measured, `[[`, character(1), "pair_id"))
    for (pair_id in pair_ids) {
        pair <- Filter(\(slot) identical(slot$pair_id, pair_id), measured)
        binding_fields <- c(
            "subject_id", "project_id", "source_id", "estimand_id",
            "comparison_id", "block_id", "session_id", "host_id", "cache_stratum",
            "work_count", "work_unit_id", "exact_input_bytes", "input_sha256"
        )
        if (length(pair) != 2L ||
            !setequal(vapply(pair, `[[`, character(1), "arm"), c("A", "B")) ||
            !identical(vapply(pair, `[[`, integer(1), "sequence_in_pair"), 1:2) ||
            !all(vapply(pair, \(slot) identical(slot$block_id, pair_id), logical(1))) ||
            any(!vapply(binding_fields, \(field) {
                identical(pair[[1L]][[field]], pair[[2L]][[field]])
            }, logical(1))) || any(!vapply(pair, \(slot) {
                isTRUE(slot$promotion_authority)
            }, logical(1)))) {
            publicationScheduleAbort("Measured pair is malformed")
        }
    }
    combinations <- split(
        measured,
        vapply(measured, \(slot) paste(
            slot$subject_id,
            slot$estimand_id,
            slot$comparison_id,
            sep = "\x1f"
        ), character(1))
    )
    for (slots in combinations) {
        pairs <- unique(vapply(slots, `[[`, character(1), "pair_id"))
        projects <- unique(vapply(slots, `[[`, character(1), "project_id"))
        first <- Filter(\(slot) slot$sequence_in_pair == 1L, slots)
        first_arms <- vapply(first, `[[`, character(1), "arm")
        project_ids <- vapply(first, `[[`, character(1), "project_id")
        session_ids <- vapply(first, `[[`, character(1), "session_id")
        project_counts <- table(project_ids)
        project_arm_balance <- vapply(projects, \(project_id) {
            arms <- first_arms[project_ids == project_id]
            sum(arms == "A") == sum(arms == "B")
        }, logical(1))
        project_session_coverage <- vapply(projects, \(project_id) {
            length(unique(session_ids[project_ids == project_id])) ==
                schedule$session_count
        }, logical(1))
        if (length(pairs) != schedule$pair_count || length(projects) < 3L ||
            length(unique(session_ids)) != schedule$session_count ||
            length(unique(as.integer(project_counts))) != 1L ||
            !all(project_arm_balance) || !all(project_session_coverage) ||
            sum(first_arms == "A") != sum(first_arms == "B")) {
            publicationScheduleAbort("Measured schedule is not balanced")
        }
    }
    invisible(TRUE)
}
