test_that("schedule is deterministic project-balanced and arm-counterbalanced", {
    set.seed(123L)
    previous_seed <- .Random.seed
    first <- publicationBuildSchedule(
        projects = publicationScheduleFixtureProjects(),
        estimand_ids = "estimand-1",
        comparison_ids = "backend-effect",
        pair_count = 30L,
        host_id = "primary-host",
        work_bindings = publicationScheduleFixtureWorkBindings()
    )
    second <- publicationBuildSchedule(
        projects = publicationScheduleFixtureProjects(),
        estimand_ids = "estimand-1",
        comparison_ids = "backend-effect",
        pair_count = 30L,
        host_id = "primary-host",
        work_bindings = publicationScheduleFixtureWorkBindings()
    )

    expect_silent(publicationValidateSchedule(first))
    expect_identical(first$schedule_digest, second$schedule_digest)
    expect_identical(first$slots, second$slots)
    expect_identical(.Random.seed, previous_seed)
    expect_length(first$slots, 72L)

    measured <- Filter(\(slot) !slot$warmup, first$slots)
    expect_length(measured, 60L)
    expect_setequal(vapply(measured, `[[`, character(1), "project_id"),
        c("project-1", "project-2", "project-3"))
    expect_setequal(vapply(measured, `[[`, character(1), "session_id"),
        paste0("session-", 1:3))
    expect_false(any(grepl(
        "memory|artifact",
        vapply(measured, `[[`, character(1), "arm"),
        ignore.case = TRUE
    )))
    expect_true(all(vapply(measured, \(slot) {
        publicationScalarNumber(slot$work_count, positive = TRUE) &&
            publicationScalarNumber(slot$exact_input_bytes, positive = TRUE) &&
            grepl("^[0-9a-f]{64}$", slot$input_sha256)
    }, logical(1))))
})

test_that("schedule rejects undercoverage imbalance duplicates and digest drift", {
    projects <- publicationScheduleFixtureProjects()

    expect_error(
        publicationBuildSchedule(
            projects = projects,
            estimand_ids = "estimand-1",
            comparison_ids = "backend-effect",
            pair_count = 32L,
            host_id = "primary-host",
            work_bindings = publicationScheduleFixtureWorkBindings()
        ),
        class = "multischolar_publication_schedule_error"
    )

    duplicate <- projects
    duplicate[[1L]][[3L]]$source_id <- "source-1"
    expect_error(
        publicationBuildSchedule(
            projects = duplicate,
            estimand_ids = "estimand-1",
            comparison_ids = "backend-effect",
            pair_count = 30L,
            host_id = "primary-host",
            work_bindings = publicationScheduleFixtureWorkBindings()
        ),
        class = "multischolar_publication_schedule_error"
    )

    schedule <- publicationBuildSchedule(
        projects = projects,
        estimand_ids = "estimand-1",
        comparison_ids = "backend-effect",
        pair_count = 30L,
        host_id = "primary-host",
        work_bindings = publicationScheduleFixtureWorkBindings()
    )
    valid_schedule <- publicationGovernanceCopy(schedule)
    schedule$slots[[1L]]$arm <- "artifact"
    expect_error(
        publicationValidateSchedule(schedule),
        class = "multischolar_publication_schedule_error"
    )

    missing_work <- publicationScheduleFixtureWorkBindings()
    missing_work[[1L]] <- NULL
    expect_error(
        publicationBuildSchedule(
            projects = projects,
            estimand_ids = "estimand-1",
            comparison_ids = "backend-effect",
            pair_count = 30L,
            host_id = "primary-host",
            work_bindings = missing_work
        ),
        class = "multischolar_publication_schedule_error"
    )

    no_warmup <- valid_schedule
    no_warmup$slots <- Filter(\(slot) !isTRUE(slot$warmup), no_warmup$slots)
    for (index in seq_along(no_warmup$slots)) {
        no_warmup$slots[[index]]$ordinal <- as.integer(index)
    }
    no_warmup$work_binding_digest <- publicationScheduleWorkDigest(
        no_warmup$slots
    )
    no_warmup$schedule_digest <- publicationScheduleDigest(no_warmup)
    expect_error(
        publicationValidateSchedule(no_warmup),
        class = "multischolar_publication_schedule_error"
    )
})

test_that("schedule supports the frozen 60-pair upper bound", {
    schedule <- publicationBuildSchedule(
        projects = publicationScheduleFixtureProjects(),
        estimand_ids = "estimand-1",
        comparison_ids = "backend-effect",
        pair_count = 60L,
        host_id = "primary-host",
        work_bindings = publicationScheduleFixtureWorkBindings()
    )

    expect_silent(publicationValidateSchedule(schedule))
    measured <- Filter(\(slot) !slot$warmup, schedule$slots)
    expect_length(measured, 120L)
})

test_that("precision-selected 36-pair schedule is project and arm balanced", {
    schedule <- publicationBuildSchedule(
        projects = publicationScheduleFixtureProjects(),
        estimand_ids = "estimand-1",
        comparison_ids = "backend-effect",
        pair_count = 36L,
        host_id = "primary-host",
        work_bindings = publicationScheduleFixtureWorkBindings()
    )
    measured <- Filter(\(slot) !slot$warmup, schedule$slots)
    first <- Filter(\(slot) slot$sequence_in_pair == 1L, measured)

    expect_silent(publicationValidateSchedule(schedule))
    expect_length(measured, 72L)
    project_counts <- table(vapply(first, `[[`, character(1), "project_id"))
    arm_counts <- table(vapply(first, `[[`, character(1), "arm"))
    expect_identical(names(project_counts), paste0("project-", 1:3))
    expect_identical(as.integer(project_counts), rep(12L, 3L))
    expect_identical(names(arm_counts), c("A", "B"))
    expect_identical(as.integer(arm_counts), c(18L, 18L))
})
