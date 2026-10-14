test_that("cohort-size experiments retain the 16 GiB guard and short runtime", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.R"), api)
    path <- withr::local_tempfile()
    for (budget in c("256", "768", "768-threads8")) {
        fields <- c("dia-controller-request-v1", "multischolar-dia003-browser-smoke.ABC123.service",
            "reader-check", "/fixture/library", "17179869184", "0", "900", "Rscript",
            "--internal-reader-check", "/fixture/library", "/fixture/output", "/fixture/repo",
            paste0("gui-cohort-size:", budget))
        writeBin(fields, path)
        request <- api$readCheckpointControllerRequest(path)
        expect_equal(request$memory_max_bytes, 17179869184)
        expect_equal(request$runtime_max_seconds, 900)
        fields[[7L]] <- "5400"
        writeBin(fields, path)
        expect_error(api$readCheckpointControllerRequest(path), "invalid controller request")
    }
})

test_that("the 300 checkpoint replay has an exact separate 90-minute guard", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.R"), api)
    path <- withr::local_tempfile()
    fields <- c("dia-controller-request-v1", "multischolar-dia003-browser-smoke.ABC123.service",
        "reader-check", "/fixture/library", "17179869184", "0", "5400", "Rscript",
        "--internal-reader-check", "/fixture/library", "/fixture/output", "/fixture/repo",
        "gui-raw-replay:300")
    writeBin(fields, path)
    request <- api$readCheckpointControllerRequest(path)
    expect_equal(request$runtime_max_seconds, 5400)
    expect_equal(request$memory_max_bytes, 17179869184)
    for (change in list(c(5L, "8589934592"), c(6L, "1"), c(7L, "7200"),
        c(length(fields), "gui-raw-replay:small"))) {
        wrong <- fields
        wrong[[as.integer(change[[1L]])]] <- change[[2L]]
        writeBin(wrong, path)
        expect_error(api$readCheckpointControllerRequest(path), "invalid controller request")
    }
})

test_that("the import design diagnostic stops before Save without changing normal routing", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), api)
    design <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))@design_matrix
    path <- tempfile(fileext = ".tsv")
    withr::defer(unlink(path))
    utils::write.table(data.frame(sample = design$Run, group = design$group), path,
        sep = "\t", row.names = FALSE, quote = FALSE)
    clicks <- character()
    driver <- list(get_js = function(script) { clicks <<- c(clicks, script); TRUE })
    api$checkpointSetInput <- function(...) NULL
    api$checkpointWaitForJs <- function(...) NULL
    api$checkpointWaitForIdle <- function(...) NULL
    api$checkpointChooseExisting <- function(driver, id, values, multiple = FALSE) values
    finished <- 0L
    api$waitCheckpointStepComplete <- function(...) { finished <<- finished + 1L }
    api$runCheckpointDesignActions(driver, path, stop_before_save = TRUE)
    expect_false(any(grepl("save_results", clicks, fixed = TRUE)))
    expect_identical(finished, 0L)
    api$runCheckpointDesignActions(driver, path)
    expect_true(any(grepl("save_results", clicks, fixed = TRUE)))
    expect_identical(finished, 1L)
})

test_that("diagnostic memory acknowledgements require an idle owner and deduplicate requests", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_gui_diagnostic.R"), api)
    root <- withr::local_tempdir()
    events <- character()
    api$checkpointGuiDiagnosticEvent <- function(...) NULL
    api$checkpointAllocatorBoundary <- function(context, label, collect = FALSE, trim = FALSE) {
        events <<- c(events, label)
    }
    shiny::testServer(function(input, output, session) {
        workflow_data <- shiny::reactiveValues()
        manager <- new.env(parent = emptyenv())
        class(manager) <- "DIAStateAdapter"
        manager$.context <- new.env(parent = emptyenv())
        manager$.context$scope_depth <- 1L
        manager$.context$store <- list(open_resources = 0L, table_action = NULL)
        manager$.context$snapshot <- list(current_state = "initial")
        manager$.context$working <- new.env(parent = emptyenv())
        session$userData$dia_capacity_manager <- manager
        context <- new.env(parent = emptyenv())
        context$directory <- root
        context$design_memory <- TRUE
        module_frame <- new.env(parent = environment())
        module_frame$session <- session$makeScope("import")
        module_frame$workflow_data <- workflow_data
        api$checkpointGuiObserveOwners(context, module_frame)
    }, {
        session$setInputs(.dia_memory_probe = "design_setup")
        expect_length(events, 0L)
        manager$.context$scope_depth <- 0L
        session$flushReact()
        expect_identical(events, c("design_setup:observed", "design_setup:gc", "design_setup:trim"))
        expect_identical(jsonlite::read_json(file.path(root, "memory-boundary-design_setup.json"))$status,
            "complete")
        session$flushReact()
        expect_length(events, 3L)
    })
})

test_that("explicit native launch captures reclamation without changing the default", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), api)
    library <- dirname(find.package("MultiScholaR"))
    root <- withr::local_tempdir()
    default <- api$checkpointInstalledApp(library, root, "qs2_duckdb", "native_v1")
    opted <- api$checkpointInstalledApp(library, root, "qs2_duckdb", "native_v1", TRUE)
    expect_false(environment(default)$fixed_reclaim)
    expect_true(environment(opted)$fixed_reclaim)
    expect_error(api$checkpointInstalledApp(library, root, "qs2_duckdb", "legacy", TRUE),
        "invalid reclamation")
    expect_error(api$checkpointInstalledApp(library, root, "qs2_duckdb", "native_v1", NA),
        "invalid reclamation")
})

test_that("timed waits observe application completion and fail promptly on session loss", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), api)
    withr::local_options(multischolar.dia_memory_audit_directory = NULL)
    launch <- function() shiny::shinyApp(
        shiny::fluidPage(shiny::actionButton("finish", "Finish"),
            shiny::actionButton("disconnect", "Disconnect"), shiny::textOutput("result")),
        function(input, output, session) {
            done <- shiny::reactiveVal("waiting")
            shiny::observeEvent(input$finish, done("complete"))
            shiny::observeEvent(input$disconnect, session$close())
            output$result <- shiny::renderText(done())
        })
    environment(launch) <- baseenv()
    browser <- shinytest2::AppDriver$new(app_dir = launch, name = "reference-wait-lifecycle",
        timeout = 20000, load_timeout = 20000, check_names = FALSE)
    withr::defer(browser$stop())
    browser$get_js("document.getElementById('finish').click();true")
    expect_no_error(api$checkpointWaitForJs(browser,
        "document.getElementById('result').textContent==='complete'", timeout = 10000))
    browser$get_js("document.getElementById('disconnect').click();true")
    start <- unname(proc.time()[["elapsed"]])
    expect_error(api$checkpointWaitForJs(browser, "false", timeout = 10000), "disconnected")
    expect_lt(unname(proc.time()[["elapsed"]]) - start, 10)
})

test_that("design assignment flushes buffered values before clicking the production button", {
    skipIfNoBrowserRuntime()
    shell <- readLines(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.sh"))
    browser_selector <- shell[grepl('if [[ "$proof_memory_guard" != true', shell, fixed = TRUE)]
    expect_length(browser_selector, 1L)
    expect_match(browser_selector, '"$proof_sequence_selector" == gui-reference-first-qc', fixed = TRUE)
    expect_match(browser_selector, '"$proof_sequence_selector" == gui-native-session', fixed = TRUE)
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), api)
    browser <- chromote::ChromoteSession$new()
    withr::defer(browser$close())
    setup <- paste0("document.body.innerHTML=",
        jsonlite::toJSON(paste0('<select id="design-tech_rep_samples"></select>',
            '<input id="design-manual_replicate_number" value="2">',
            '<div id="design-tech_rep_assignment_mode"><input type="radio" value="manual" checked></div>',
            '<button id="design-assign_tech_reps">Assign</button>'), auto_unbox = TRUE), ";",
        "document.getElementById('design-tech_rep_samples').selectize={getValue:()=>['s3','s4']};",
        "window.sent=[];window.committed={};window.clicks=[];",
        "window.Shiny={setInputValue:(id,value,options)=>{",
        "sent.push({id,value,priority:options.priority});",
        "if(options.priority==='event')committed[id]=value;}};",
        "document.getElementById('design-assign_tech_reps').onclick=()=>clicks.push(JSON.parse(JSON.stringify(committed)));",
        "document.getElementById('design-assign_tech_reps').click();")
    expect_null(browser$Runtime$evaluate(setup)$exceptionDetails)
    script <- api$checkpointDesignAssignmentScript("design-", c("s3", "s4"), 2)
    result <- browser$Runtime$evaluate(script, returnByValue = TRUE)
    expect_null(result$exceptionDetails)
    expect_true(result$result$value)
    observed <- jsonlite::fromJSON(browser$Runtime$evaluate("JSON.stringify({sent,clicks})",
        returnByValue = TRUE)$result$value, simplifyVector = FALSE)
    expect_length(observed$clicks[[1L]], 0L) # An ordinary click can precede input delivery.
    expect_identical(unlist(observed$clicks[[2L]]$`design-tech_rep_samples`), c("s3", "s4"))
    expect_identical(observed$clicks[[2L]]$`design-manual_replicate_number`, 2L)
    expect_identical(observed$clicks[[2L]]$`design-tech_rep_assignment_mode`, "manual")
    expect_true(all(vapply(observed$sent, \(entry) entry$priority == "event", logical(1))))
    browser$Runtime$evaluate("document.getElementById('design-manual_replicate_number').value='3'")
    expect_match(browser$Runtime$evaluate(script)$exceptionDetails$exception$description,
        "controls differ from frozen input")
    expect_equal(browser$Runtime$evaluate("clicks.length", returnByValue = TRUE)$result$value, 2)
    browser$Runtime$evaluate(paste0("window.selection=[];",
        "document.getElementById('design-tech_rep_samples').selectize={",
        "options:{s3:{},s4:{}},getValue:()=>selection,setValue:(x,silent)=>{selection=x;}};"))
    prepared <- api$checkpointDesignAssignmentScript("design-", c("s3", "s4"), 2, prepare = TRUE)
    expect_null(browser$Runtime$evaluate(prepared)$exceptionDetails)
    observed <- jsonlite::fromJSON(browser$Runtime$evaluate("JSON.stringify(clicks[2])",
        returnByValue = TRUE)$result$value)
    expect_identical(observed$`design-tech_rep_samples`, c("s3", "s4"))
    expect_equal(observed$`design-manual_replicate_number`, 2)
    absent <- api$checkpointDesignAssignmentScript("design-", c("s3", "missing"), 2, prepare = TRUE)
    expect_match(browser$Runtime$evaluate(absent)$exceptionDetails$exception$description, "option absent")
    expect_equal(browser$Runtime$evaluate("clicks.length", returnByValue = TRUE)$result$value, 3)
})

test_that("saved-design validation rejects the retained assignment failure before QC", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), api)
    root <- withr::local_tempdir()
    expected <- data.frame(sample = paste0("s", 1:6), group = rep(c("control", "treated"), c(4, 2)),
        replicates = c(1, 1, 2, 2, 1, 1), tech_rep_group = c("control_1", "control_1",
            "control_2", "control_2", "treated_1", "treated_1"))
    observed <- expected
    names(observed)[1L] <- "Run"
    observed$tech_reps <- rep(1:2, 3)
    input <- file.path(root, "input.tsv")
    saved <- file.path(root, "design_matrix.tab")
    utils::write.table(expected, input, sep = "\t", quote = FALSE, row.names = FALSE)
    write_saved <- function(value) utils::write.table(value, saved, sep = "\t", quote = FALSE, row.names = FALSE)
    write_saved(observed)
    expect_identical(api$checkCheckpointSavedDesign(input, saved)$status, "passed_before_qc")
    for (field in c("replicates", "tech_reps", "tech_rep_group")) {
        changed <- observed
        changed[[field]][1L] <- NA
        write_saved(changed)
        expect_error(api$checkCheckpointSavedDesign(input, saved), paste0("input: ", field))
    }
    write_saved(observed[6:1, ])
    expect_error(api$checkCheckpointSavedDesign(input, saved), "identities/order")
    receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    if (nzchar(receipt)) {
        proof <- dirname(dirname(receipt))
        declared <- file.path(proof, "browser-smoke.ykIZwA", "design.tsv")
        for (member in c("vOSqsa", "SfeXyz")) {
            path <- file.path(proof, paste0("browser-smoke.", member), "dia_baseline_smoke",
                "scripts", "proteomics", "design_matrix.tab")
            if (file.exists(declared) && file.exists(path)) {
                if (member == "vOSqsa") expect_identical(api$checkCheckpointSavedDesign(declared, path)$status,
                    "passed_before_qc") else expect_error(api$checkCheckpointSavedDesign(declared, path), "replicates")
            }
        }
    }
})

test_that("sample-scaled diagnostics retain their separate recipes and exact 18 GiB limits", {
    driver <- new.env(parent = baseenv())
    runner <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), driver)
    sys.source(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.R"), runner)
    spec <- driver$readCheckpointRuvPearsonSpec(normalizePath(test_path("..", "..")))
    expect_identical(vapply(spec$recipes, function(x) as.numeric(x$n_samples), numeric(1)), c(24, 96, 192))
    expect_identical(vapply(spec$recipes, function(x) as.numeric(x$precursor_rows), numeric(1)),
        c(192000, 7680000, 15360000))
    expect_identical(spec$settings$ruv, "automatic")
    expect_identical(spec$settings$correlation, "apply")
    expect_identical(spec$stress_96$recipe, spec$recipes[[2L]])
    expect_identical(spec$final_launch_headroom_bytes, 16 * 1024^3)
    expect_equal(spec$stress_96$limits, spec$stress_300$limits)
    expect_identical(spec$integrated_candidate_rdb_sha256,
        "060455e405f7b7dbd5ee723b2f3e3ba8eff9dc4a4f4833f575ff463696e240d0")
    expect_equal(spec$stress_500$recipe$n_samples, 500)
    expect_equal(spec$stress_500$recipe$precursor_rows, 40000000)
    expect_equal(spec$stress_500$limits$memory_max_bytes, 16 * 1024^3)
    expect_equal(spec$stress_500$limits$early_stop_bytes, 15.5 * 1024^3)
    chunk <- makeWorkflowPrecursorChunk(spec$stress_500$recipe, 1:50)
    expect_equal(nrow(chunk), 200000)
    expect_equal(length(unique(chunk$Run)), 500)
    expect_equal(nrow(makeWorkflowCheckpointDesign(500)), 500)
    expect_error(makeWorkflowPrecursorChunk(spec$stress_500$recipe, 1:51), "200000-row")
    expect_equal(spec$stress_300$recipe$n_samples, 300)
    expect_equal(spec$stress_300$recipe$precursor_rows, 24000000)
    expect_equal(spec$stress_300$limits, spec$stress_500$limits)
    chunk <- makeWorkflowPrecursorChunk(spec$stress_300$recipe, 1:83)
    expect_equal(nrow(chunk), 199200)
    expect_equal(length(unique(chunk$Run)), 300)
    expect_equal(nrow(makeWorkflowCheckpointDesign(300)), 300)
    expect_error(makeWorkflowPrecursorChunk(spec$stress_300$recipe, 1:84), "200000-row")
    rm(chunk)
    path <- tempfile()
    fields <- c("dia-controller-request-v1", "multischolar-dia003-browser-smoke.ABC123.service",
        "browser-baseline", "library", "19327352832", "0", "3600", "Rscript", "--vanilla",
        "runner.R", "--internal-browser-baseline", "library", "run", "report", "fasta", "design",
        "capacity:qs2_duckdb:ruv-pearson")
    writeBin(fields, path)
    limits <- runner$readCheckpointControllerRequest(path)
    expect_equal(limits$memory_max_bytes, 18 * 1024^3)
    expect_equal(limits$runtime_max_seconds, 3600)
    verifier <- fields[seq_len(14L)]
    verifier[[3L]] <- "verify-da"
    verifier[[11L]] <- "--internal-verify-da"
    verifier[[14L]] <- "member"
    verifier <- c(verifier, "diagnostic:ruv-pearson")
    writeBin(verifier, path)
    expect_equal(runner$readCheckpointControllerRequest(path)$memory_max_bytes, 18 * 1024^3)
    pilot <- fields
    pilot[c(5L, 7L)] <- c("8589934592", "900")
    pilot[[length(pilot)]] <- "capacity:qs2_duckdb:ruv-pearson-small"
    writeBin(pilot, path)
    expect_equal(runner$readCheckpointControllerRequest(path)$memory_max_bytes, 8 * 1024^3)
    expect_equal(runner$readCheckpointControllerRequest(path)$runtime_max_seconds, 900)
    for (suffix in c(":native", ":timed", ":timed:native")) {
        modified <- fields
        modified[[length(modified)]] <- paste0("capacity:qs2_duckdb:ruv-pearson", suffix)
        writeBin(modified, path)
        expect_equal(runner$readCheckpointControllerRequest(path)$memory_max_bytes, 18 * 1024^3)
        modified <- verifier
        modified[[length(modified)]] <- paste0("diagnostic:ruv-pearson", suffix)
        writeBin(modified, path)
        expect_equal(runner$readCheckpointControllerRequest(path)$runtime_max_seconds, 3600)
    }
    for (bytes in c("25769803776", "30064771072")) {
        changed <- fields
        changed[[5L]] <- bytes
        writeBin(changed, path)
        expect_error(runner$readCheckpointControllerRequest(path), "identity or limits")
    }
    for (samples in c(96, 300, 500)) for (case in list(fields, verifier)) {
        case[[5L]] <- "17179869184"
        case[[length(case)]] <- if (case[[3L]] == "browser-baseline")
            paste0("capacity:qs2_duckdb:ruv-pearson-", samples, ":timed:native") else
                paste0("diagnostic:ruv-pearson-", samples, ":timed:native")
        writeBin(case, path)
        expect_equal(runner$readCheckpointControllerRequest(path)$memory_max_bytes, 16 * 1024^3)
        expect_equal(runner$readCheckpointControllerRequest(path)$runtime_max_seconds, 3600)
        case[[5L]] <- "19327352832"
        writeBin(case, path)
        expect_error(runner$readCheckpointControllerRequest(path), "identity or limits")
    }
    fields[[length(fields)]] <- "capacity:qs2_duckdb"
    writeBin(fields, path)
    expect_error(runner$readCheckpointControllerRequest(path), "identity or limits")
})

test_that("large-sample shell routing uses 16 GiB for launch and verification and stops at 15.5", {
    skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))
    root <- withr::local_tempdir()
    bin <- file.path(root, "bin")
    library_path <- file.path(root, "library")
    dir.create(bin)
    dir.create(library_path)
    inputs <- file.path(root, c("report.tsv", "proteins.fasta", "design.tsv"))
    for (path in inputs) writeLines("routing fixture only", path)
    launch_log <- file.path(root, "launch-arguments")
    stop_log <- file.path(root, "stopped-unit")
    writeLines(c("#!/usr/bin/env bash", "exit 0"), file.path(bin, "Rscript"))
    writeLines(c("#!/usr/bin/env bash", 'printf "%s\\n" "$@" > "$DIA_500_LAUNCH_LOG"'),
        file.path(bin, "systemd-run"))
    writeLines(c("#!/usr/bin/env bash", 'echo "${DIA_500_HEADROOM:-21474836480}"'),
        file.path(bin, "awk"))
    writeLines(c("#!/usr/bin/env bash", 'if [[ "$2" == show ]]; then',
        'if [[ "${DIA_500_GUARD_TRIP:-0}" == 1 ]]; then',
        'printf "Id=%s\\nSubState=running\\nMemoryCurrent=12000000000\\nMemoryPeak=16642998272\\n" "$3"',
        'else printf "Id=%s\\nSubState=exited\\nMainPID=0\\nMemoryCurrent=0\\nMemoryPeak=1048576\\nResult=success\\nExecMainStatus=0\\n" "$3"; fi',
        'elif [[ "$2" == stop ]]; then printf "%s\\n" "$3" >> "$DIA_500_STOP_LOG"; else exit 3; fi'),
        file.path(bin, "systemctl"))
    Sys.chmod(file.path(bin, c("Rscript", "systemd-run", "systemctl", "awk")), "0755")
    withr::local_envvar(c(PATH = paste(bin, Sys.getenv("PATH"), sep = .Platform$path.sep),
        CHROMOTE_CHROME = Sys.which("true"), DIA_500_LAUNCH_LOG = launch_log,
        DIA_500_STOP_LOG = stop_log))
    runner <- normalizePath(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.sh"))
    api <- new.env(parent = baseenv())
    sys.source(sub("[.]sh$", ".R", runner), api)
    for (samples in c(96, 300, 500)) for (operation in c("browser-baseline", "verify-da"))
        for (first_qc in if (samples != 500 && operation == "browser-baseline") c(FALSE, TRUE) else FALSE)
        for (extended in if (samples == 300 && operation == "browser-baseline" && !first_qc) c(FALSE, TRUE) else FALSE)
        for (trip in c(FALSE, TRUE)) {
        before <- list.dirs(root, recursive = FALSE)
        result <- withr::with_envvar(c(DIA_500_GUARD_TRIP = if (trip) "1" else "0",
            DIA_500_HEADROOM = if (samples == 500) "21474836480" else "18253611008"),
            processx::run("bash", c(runner, "--service", operation, library_path, root,
                if (operation == "browser-baseline") c(inputs, "--capacity-qs2_duckdb") else library_path,
                paste0("--stress-", samples), "--native-session", if (first_qc) "--first-qc",
                if (extended) c("--reclaim-memory", "--runtime-120-minutes")), error_on_status = FALSE))
        directory <- setdiff(list.dirs(root, recursive = FALSE), before)
        expect_length(directory, 1L)
        expect_identical(result$status, if (trip) 76L else 0L)
        request <- api$readCheckpointControllerRequest(file.path(directory, "controller-request.nul"))
        expect_equal(request$memory_max_bytes, 16 * 1024^3)
        expect_equal(request$swap_max_bytes, 0)
        expect_equal(request$runtime_max_seconds, if (extended) 7200 else 3600)
        expect_identical(endsWith(utils::tail(request$command, 1L), ":first-qc"), first_qc)
        launched <- readLines(launch_log)
        expect_true(all(c("--property=MemoryMax=17179869184", "--property=MemorySwapMax=0",
            paste0("--property=RuntimeMaxSec=", if (extended) 7200 else 3600)) %in% launched))
        cap <- if (operation == "browser-baseline") "GUI" else "VERIFY"
        expect_true(paste0("--setenv=MULTISCHOLAR_DIA_SMALL_", cap, "_CAP=17179869184") %in% launched)
        expect_identical(utils::tail(readLines(stop_log), 1L), request$unit)
        if (trip) expect_match(paste(readLines(file.path(directory, "guard-stop.txt")), collapse = "\n"),
            "reason=early_memory_stop")
    }
    for (samples in c(96, 300, 500)) {
        denied <- withr::with_envvar(c(DIA_500_HEADROOM =
            if (samples == 500) "18253611008" else "16642998272"), processx::run("bash",
            c(runner, "--service", "browser-baseline", library_path, root, inputs,
                "--capacity-qs2_duckdb", paste0("--stress-", samples), "--native-session"), error_on_status = FALSE))
        expect_identical(denied$status, 76L)
        expect_match(denied$stderr, paste("headroom is below", if (samples == 500) 18 else 16, "GiB"))
    }
})

test_that("final 96 and 300 routing binds retained inputs to the qualified builds", {
    receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    skip_if(!nzchar(receipt), "guarded retained-input routing preflight only")
    proof <- dirname(dirname(receipt))
    # Historical complete-journey approval remains bound to the frozen build.
    # The new first-QC candidate does not inherit that approval.
    candidate <- file.path(proof, "browser-smoke.x3Dkg2", "library")
    control <- file.path(proof, "browser-smoke.xz9cb6", "library")
    inputs96 <- file.path(proof, "browser-smoke.nDqUl6", "inputs.json")
    inputs300 <- file.path(proof, "browser-smoke.OQADra", "inputs.json")
    expect_true(all(file.exists(c(inputs96, inputs300))))
    runner <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.R"), runner)
    runner$checkpointRunnerDirectory <- function()
        normalizePath(test_path("..", "..", "tools", "profiling"))
    runner$checkpointServiceMemory <- function(memory_max) {
        expect_identical(memory_max, "17179869184")
    }
    runner$runCheckpointBrowserSmoke <- function(library_path, run_dir, inputs) {
        expect_identical(inputs$profile_application, isTRUE(inputs$design_memory))
        expect_identical(inputs$session_format, "native_v1")
        expect_identical(inputs$capacity$treatment, "qs2_duckdb")
        inputs
    }
    run_dir <- withr::local_tempdir()
    arguments <- function(library, input, samples) {
        metadata <- jsonlite::read_json(input)
        c("--internal-browser-baseline", library, run_dir,
            unlist(metadata$files, use.names = FALSE),
            paste0("capacity:qs2_duckdb:ruv-pearson-", samples, ":timed:native"))
    }
    for (library in c(control, candidate)) {
        result <- runner$checkpointProofMain(arguments(library, inputs96, 96))
        expect_identical(result$expected_sha256, jsonlite::read_json(inputs96)$sha256)
    }
    current <- dirname(find.package("MultiScholaR"))
    spec <- jsonlite::read_json(test_path("..", "testdata", "workflow-checkpoints", "gui-ruv-pearson-diagnostic.json"))
    if (digest::digest(file = file.path(current, "MultiScholaR", "R", "MultiScholaR.rdb"),
        algo = "sha256") %in% unlist(spec$native_first_300_rdb_sha256)) {
        request <- arguments(current, inputs300, 300)
        request[[7L]] <- paste0(request[[7L]], ":reclaim")
        opted <- runner$checkpointProofMain(request)
        expect_true(opted$reclaim_memory)
        expect_false(opted$profile_application)
        expect_identical(opted$expected_sha256, jsonlite::read_json(inputs300)$sha256)
        request[[7L]] <- sub(":reclaim$", ":design-memory", request[[7L]])
        diagnostic <- runner$checkpointProofMain(request)
        expect_true(diagnostic$design_memory)
        expect_true(diagnostic$profile_application)
        expect_false(diagnostic$peptide_diagnostic)
        expect_false(diagnostic$summary_diagnostic)
        expect_false(diagnostic$reclaim_memory)
    }
    current_hash <- digest::digest(file = file.path(current, "MultiScholaR", "R", "MultiScholaR.rdb"), algo = "sha256")
    if (current_hash %in% unlist(spec$extended_120_rdb_sha256)) {
        request <- arguments(current, inputs300, 300)
        request[[7L]] <- paste0(request[[7L]], ":reclaim:extended120")
        extended <- withr::with_envvar(c(MULTISCHOLAR_DIA_DIAGNOSTIC_SECONDS = "7200"),
            runner$checkpointProofMain(request))
        expect_true(extended$reclaim_memory)
        expect_true(extended$summary_diagnostic)
        expect_false(extended$first_qc)
        expect_false(extended$profile_application)
        wrong <- arguments(candidate, inputs300, 300)
        wrong[[7L]] <- paste0(wrong[[7L]], ":reclaim:extended120")
        expect_error(withr::with_envvar(c(MULTISCHOLAR_DIA_DIAGNOSTIC_SECONDS = "7200"),
            runner$checkpointProofMain(wrong)), "120-minute route is not approved")
    }
    result <- runner$checkpointProofMain(arguments(candidate, inputs300, 300))
    expect_identical(result$expected_sha256, jsonlite::read_json(inputs300)$sha256)
    native_control <- file.path(proof, "browser-smoke.jfR7Ma", "library")
    reference_candidate <- file.path(proof, "browser-smoke.Qb1Btf", "library")
    for (library in c(native_control, reference_candidate)) {
        result <- runner$checkpointProofMain(arguments(library, inputs96, 96))
        expect_identical(result$expected_sha256, jsonlite::read_json(inputs96)$sha256)
        expect_error(runner$checkpointProofMain(arguments(library, inputs300, 300)),
            "sample stress build, recipe or route differs")
    }
    first_qc_candidate <- Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY")
    first_qc_hash <- digest::digest(file = file.path(first_qc_candidate, "MultiScholaR", "R", "MultiScholaR.rdb"), algo = "sha256")
    spec <- jsonlite::read_json(test_path("..", "testdata", "workflow-checkpoints", "gui-ruv-pearson-diagnostic.json"))
    if (first_qc_hash %in% unlist(spec$reference_first_qc_rdb_sha256)) {
        for (samples in c(96, 300)) {
            args <- arguments(first_qc_candidate, if (samples == 96) inputs96 else inputs300, samples)
            args[[length(args)]] <- paste0(args[[length(args)]], ":first-qc")
            result <- runner$checkpointProofMain(args)
            expect_true(result$first_qc)
            expect_false(result$protein_diagnostic)
            expect_false(result$normalization_diagnostic)
            expect_false(result$da_diagnostic)
            args[[length(args)]] <- paste0(args[[length(args)]], ":reclaim")
            reclaimed <- runner$checkpointProofMain(args)
            expect_true(reclaimed$first_qc)
            expect_true(reclaimed$reclaim_memory)
            expect_false(reclaimed$profile_application)
            expect_false(reclaimed$summary_diagnostic)
            expect_error(runner$checkpointProofMain(arguments(first_qc_candidate,
                if (samples == 96) inputs96 else inputs300, samples)), "runtime or input binding differs")
        }
    }
    input500 <- file.path(proof, "browser-smoke.570cC2", "inputs.json")
    if (first_qc_hash %in% unlist(spec$native_session_rdb_sha256) &&
        !first_qc_hash %in% unlist(spec$native_first_96_rdb_sha256) &&
        !first_qc_hash %in% unlist(spec$reference_first_qc_rdb_sha256)) {
        for (samples in c(96, 300)) {
            allowed <- unlist(spec[[paste0("native_first_", samples, "_rdb_sha256")]])
            if (!first_qc_hash %in% allowed) expect_error(runner$checkpointProofMain(
                arguments(first_qc_candidate, if (samples == 96) inputs96 else inputs300, samples)),
                "sample stress build, recipe or route differs")
        }
    }
    expect_error(runner$checkpointProofMain(arguments(candidate, input500, 500)),
        "sample stress build, recipe or route differs")
    expect_error(runner$checkpointProofMain(arguments(candidate, inputs300, 96)),
        "sample stress build, recipe or route differs")
})

test_that("the independent verifier requires the real RUV history and preserves the old route", {
    runner <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.R"), runner)
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    states <- c("raw_data_s4", "qvalue_filtered", "protein_peptide_filtered", "precursor_rollup",
        "intensity_filtered", "sample_filtered", "replicate_filtered", "imputed", "protein_s4_created",
        "protein_accession_cleaned", "protein_intensity_filtered", "duplicates_removed",
        "protein_replicate_filtered", "ruv_corrected", "correlation_filtered")
    manager <- WorkflowState$new(audit_enabled = FALSE)
    for (name in states) manager$saveState(name, object, list(), name)
    session <- list(r6_complete_states = manager$states, r6_state_history = manager$state_history,
        r6_current_state_name = manager$current_state, current_s4_object = object,
        correlation_filtered_s4 = object, final_protein_count = nrow(object@protein_quant_table),
        final_sample_count = nrow(object@design_matrix), ruv_mode = "automatic", ruv_applied = TRUE, ruv_k = 1L)
    table <- object@protein_quant_table
    table[[object@protein_id_column]] <- as.character(table[[object@protein_id_column]])
    expect_equal(runner$checkCheckpointSessionContents(session, object, table, ruv_diagnostic = TRUE)$history_entries, 16)
    expect_error(runner$checkCheckpointSessionContents(session, object, table))
    session$ruv_applied <- FALSE
    expect_error(runner$checkCheckpointSessionContents(session, object, table, ruv_diagnostic = TRUE))
})

test_that("the diagnostic drives automatic RUV and Pearson Apply while the old route still skips", {
    driver <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), driver)
    inputs <- list()
    clicks <- character()
    driver$checkpointSetInput <- function(driver, name, value) inputs[[name]] <<- value
    driver$checkpointClickAndWaitText <- function(driver, input_id, selector, text, ...) {
        clicks <<- c(clicks, input_id)
        invisible(TRUE)
    }
    result <- driver$runCheckpointNormalizationActions(NULL, "ruv_pearson")
    expect_identical(result$ruv, "automatic")
    expect_identical(result$correlation, "apply")
    expect_equal(result$correlation_threshold, 0.5)
    expect_true(any(grepl("apply_correlation_filter$", clicks)))
    expect_false(any(grepl("skip_correlation_filter$", clicks)))
    expect_identical(inputs$`proteomics_workflow-normalization-auto_percentage_max`, 20)
    clicks <- character()
    result <- driver$runCheckpointNormalizationActions(NULL)
    expect_identical(result, list(method = "cyclicloess", ruv = "skip", correlation = "skip",
        actions = c("run_normalization", "skip_correlation_filter", "export_filtered_session")))
    expect_true(any(grepl("skip_correlation_filter$", clicks)))
})

test_that("real RUV and Pearson traces observe matrices without changing numerical results", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_gui_diagnostic.R"), api)
    ns <- asNamespace("MultiScholaR")
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    samples <- as.character(object@design_matrix[[object@sample_id]])
    set.seed(902023)
    for (sample in samples) object@protein_quant_table[[sample]] <- stats::rnorm(nrow(object@protein_quant_table), 20, 1)
    object@protein_quant_table[[samples[[1L]]]][[1L]] <- NA_real_
    ids <- object@protein_quant_table[[object@protein_id_column]]
    controls <- stats::setNames(seq_along(ids) <= length(ids) / 2, ids)
    calculate <- function() {
        set.seed(902023)
        list(correlation = get("pearsonCorForSamplePairs", ns)(object, correlation_group = "group"),
            ruv = get("ruvIII_C_Varying", ns)(object, ruv_grouping_variable = "group", ruv_number_k = 1, ctrl = controls),
            cancor = get("ruvCancor", ns)(object, ctrl = controls, num_components_to_impute = 2, ruv_grouping_variable = "group")$data)
    }
    expected <- calculate()
    root <- withr::local_tempdir()
    context <- api$startCheckpointGuiDiagnostic(root)
    stopped <- FALSE
    withr::defer(if (!stopped) api$stopCheckpointGuiDiagnostic(context))
    api$installCheckpointRuvPearsonTraces(context, ns)
    installed <- lapply(readLines(file.path(root, "application-phases.jsonl")), jsonlite::fromJSON)
    expect_false(any(startsWith(vapply(installed, `[[`, character(1), "stage"), "enter:")))
    actual <- calculate()
    expect_equal(actual$correlation, expected$correlation, tolerance = 1e-12)
    expect_equal(actual$ruv, expected$ruv, tolerance = 1e-12)
    expect_equal(actual$cancor, expected$cancor, tolerance = 1e-12)
    expect_error(get("cor", parent.env(ns))(matrix(letters[1:12], 3)), "numeric")
    api$stopCheckpointGuiDiagnostic(context)
    stopped <- TRUE
    expect_null(getOption("multischolar.ruv_pearson_observer"))
    records <- lapply(readLines(file.path(root, "application-phases.jsonl")), jsonlite::fromJSON)
    receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    if (nzchar(receipt)) expect_true(file.copy(file.path(root, "application-phases.jsonl"),
        file.path(dirname(receipt), "ruv-pearson-small-traces.jsonl"), overwrite = FALSE))
    stages <- vapply(records, `[[`, character(1), "stage")
    expect_true(all(c("enter:import:RUVIII_C_Varying", "exit:import:RUVIII_C_Varying",
        "enter:import:ruv_cancorplot", "exit:import:ruv_cancorplot",
        "enter:import:cor", "exit:import:cor",
        "enter:pearsonCorForSamplePairs[ProteinQuantitativeData]",
        "exit:pearsonCorForSamplePairs[ProteinQuantitativeData]") %in% stages))
    event <- records[[which(stages == "enter:import:RUVIII_C_Varying")[[1L]]]]
    expect_equal(unlist(event$detail$Y$shape), c(length(samples), length(ids)))
    expect_gt(event$memory$VmRSS, 0)
    expect_lt(file.size(file.path(root, "application-phases.jsonl")), 65536)
    expect_false(any(c("object", "frame", "Y", "M") %in% names(as.list(context))))
    expect_false(file.exists(file.path(root, "application.Rprofmem")))
    evidence <- api$readCheckpointRuvPearsonEvidence(root)
    expect_true(evidence$coverage_complete)
    expect_equal(evidence$calls[["import:RUVIII_C_Varying"]], c(enter = 1L, exit = 1L))
    expect_equal(unlist(evidence$first_call_shapes[["import:cor"]]$x$shape), c(length(ids), length(samples)))
})

test_that("the retained automatic RUV plot is portable through export validation", {
    receipt_path <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    skip_if(!nzchar(receipt_path), "requires the guarded DIA-025 evidence runner")
    proof_root <- dirname(dirname(receipt_path))
    result_path <- file.path(proof_root, "browser-smoke.Y9atKE", "dia_baseline_smoke",
        "scripts", "proteomics", "ruv_optimization_results.RDS")
    skip_if_not(file.exists(result_path), "retained Y9atKE pilot is unavailable")
    expect_lt(file.size(result_path), 1024^2)
    validator <- get("checkSavedWorkflowPayload", asNamespace("MultiScholaR"))
    expect_false(inherits(validator, "functionWithTrace"))
    expect_null(getOption("multischolar.ruv_pearson_observer"))
    result <- readRDS(result_path)
    # Capture only scalar failure locations, never call frames or payloads.
    inspectPayload <- function(value) {
        location <- list()
        error <- tryCatch(withCallingHandlers({
            validator(value)
            NULL
        }, error = function(condition) {
            calls <- sys.calls()
            frames <- sys.frames()
            indexes <- which(vapply(seq_along(calls), \(i) {
                identical(sys.function(i), validator)
            }, logical(1)))
            location <<- lapply(utils::tail(indexes, 8L), \(i) {
                frame <- frames[[i]]
                list(call = paste(deparse(calls[[i]]), collapse = " "),
                    field = if (exists("field", frame, inherits = FALSE))
                        get("field", frame, inherits = FALSE) else NULL,
                    type = tryCatch(typeof(get("object", frame, inherits = FALSE)),
                        error = \(condition) "unavailable"))
            })
        }), error = \(condition) conditionMessage(condition))
        list(error = error, location = location)
    }
    fields <- lapply(result, inspectPayload)
    plain_plot <- ggplot2::ggplot(data.frame(x = 1:3, y = 3:1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()
    plain <- inspectPayload(plain_plot)
    evidence <- list(source_member = "Y9atKE", source_bytes = file.size(result_path),
        source_sha256 = digest::digest(file = result_path, algo = "sha256"),
        traces_active = FALSE, retained_fields = fields, ordinary_plot = plain,
        promotion_authorized = FALSE)
    jsonlite::write_json(evidence, file.path(dirname(receipt_path), "ruv-export-validation.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null")
    expect_null(fields$best_cancor_plot$error)
    expect_null(plain$error)
    if (is.null(fields$best_cancor_plot$error)) {
        root <- withr::local_tempdir()
        store <- newWorkflowEntityStore(root, "ruv-export", "DIA")
        withr::defer(closeWorkflowEntityStore(store))
        identity <- list(input = evidence$source_sha256, design = "pilot-24-samples",
            parameters = "automatic", selectors = "complete", software = "repair-test")
        metadata <- list(class = "list", shape = integer(), factors = list())
        original_plot <- result$best_cancor_plot
        original_data <- ggplot2::ggplot_build(original_plot)$data
        original_png <- file.path(root, "original.png")
        ggplot2::ggsave(original_png, original_plot, width = 6, height = 4, dpi = 100)
        for (codec in c("rds", "qs2")) {
            reference <- saveWorkflowCheckpoint(store, result, codec, identity, metadata,
                codec = codec)
            restored <- hydrateSavedWorkflowCheckpoint(store, reference)
            fields_without_plot <- setdiff(names(result), "best_cancor_plot")
            expect_identical(restored[fields_without_plot], result[fields_without_plot])
            expect_identical(class(restored$best_cancor_plot), class(original_plot))
            expect_identical(ggplot2::ggplot_build(restored$best_cancor_plot)$data, original_data)
            rendered <- file.path(root, paste0(codec, ".png"))
            ggplot2::ggsave(rendered, restored$best_cancor_plot, width = 6, height = 4, dpi = 100)
            expect_identical(digest::digest(file = rendered, algo = "sha256"),
                digest::digest(file = original_png, algo = "sha256"))
        }
    }
})

test_that("missing plot bindings do not bypass payload resource checks", {
    validator <- get("checkSavedWorkflowPayload", asNamespace("MultiScholaR"))
    payload <- new.env(parent = emptyenv())
    payload$missing_argument <- rlang::missing_arg()
    payload$`...` <- rlang::missing_arg()
    payload$cycle <- payload
    expect_true(validator(payload))
    expect_true(validator(list(rlang::missing_arg(), formals(function(x, ...) NULL))))
    expect_true(rlang::is_missing(.subset2(payload, "...")))
    connection <- file(withr::local_tempfile(), "w")
    withr::defer(close(connection))
    for (resource in list(connection, methods::new("externalptr"),
        rlang::new_weakref(new.env(parent = emptyenv())))) {
        payload$resource <- resource
        expect_error(validator(payload), class = "saved_entity_unsupported")
    }
    rm("resource", envir = payload)
    touched <- FALSE
    makeActiveBinding("active", function() { touched <<- TRUE; 1L }, payload)
    expect_error(validator(payload), class = "saved_entity_unsupported")
    expect_false(touched)
    rm("active", envir = payload)
    delayedAssign("lazy", { touched <<- TRUE; connection }, assign.env = payload)
    expect_error(validator(payload), class = "saved_entity_unsupported")
    expect_false(touched)
    rm("lazy", envir = payload)
    # Class methods must not hide ordinary environment fields from inspection.
    class(payload) <- "payload_validation_fixture"
    payload$resource <- connection
    expect_error(validator(payload), class = "saved_entity_unsupported")
})

test_that("the real RUV inventory requires its actual filenames and extra figures", {
    receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    skip_if(!nzchar(receipt), "requires guarded retained GUI evidence")
    member <- file.path(dirname(dirname(receipt)), "browser-smoke.vOSqsa")
    skip_if_not(file.exists(file.path(member, "smoke-receipt.json")))
    runner <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.R"), runner)
    inventory <- runner$checkCheckpointIqQcOutputs(member, ruv_diagnostic = TRUE)
    expect_identical(inventory$status, "passed_iq_qc_presence_only")
    expect_true(all(paste0("ruv_normalised_results_cln_with_replicates", c(".RDS", ".tsv")) %in%
        basename(inventory$files)))
    expect_true(all(paste0("ruv_corrected_", c("pca", "density", "rle", "correlation", "cancor"),
        ".png") %in% basename(inventory$files)))
    expect_error(runner$checkCheckpointIqQcOutputs(member), "normalised_results_cln_with_replicates")
    runner$file.info <- function(paths) {
        info <- base::file.info(paths)
        info$size[basename(paths) == "ruv_corrected_cancor.png"] <- 0
        info
    }
    expect_error(runner$checkCheckpointIqQcOutputs(member, TRUE), "ruv_corrected_cancor.png")
})

test_that("Summary verification selects exact RUV or Skip publication names", {
    receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    skip_if(!nzchar(receipt), "requires guarded retained GUI evidence")
    proof_root <- dirname(dirname(receipt))
    members <- file.path(proof_root, c("browser-smoke.vOSqsa", "browser-smoke.tCCZKy"))
    skip_if_not(all(file.exists(file.path(members, "smoke-receipt.json"))))
    for (i in seq_along(members)) {
        runner <- new.env(parent = baseenv())
        sys.source(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.R"), runner)
        ruv <- i == 1L
        expected <- if (ruv) c("RUV_normalised_results.tsv", "ruv_normalised_results.RDS") else
            c("normalised_results.tsv", "normalised_results.RDS")
        result <- runner$checkCheckpointSummaryFiles(members[[i]], "DIA Experiment", ruv)
        expect_identical(result$status, "passed_small_summary_files")
        expect_identical(basename(unlist(result$files[c("table_tsv", "table_rds")],
            use.names = FALSE)), expected)
        expect_error(runner$checkCheckpointSummaryFiles(members[[i]], "DIA Experiment", !ruv),
            "table_tsv,table_rds")
        runner$file.info <- function(paths) {
            info <- base::file.info(paths)
            info$size[basename(paths) == expected[[1L]]] <- 0
            info
        }
        expect_error(runner$checkCheckpointSummaryFiles(members[[i]], "DIA Experiment", ruv),
            "table_tsv")
    }
})


test_that("120-minute wait budgets require the explicit route option", {
    driver <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), driver)
    withr::local_options(multischolar.dia_memory_audit_seconds = 7200,
        multischolar.dia_extended120 = FALSE)
    expect_error(driver$checkpointWaitBudget(elapsed_seconds = 1), "invalid diagnostic time ceiling")
    options(multischolar.dia_extended120 = TRUE)
    expect_identical(driver$checkpointWaitBudget(elapsed_seconds = 7199), 1000L)
    expect_identical(driver$checkpointWaitBudget(control = TRUE, elapsed_seconds = 1000), 30000L)
    expect_error(driver$checkpointWaitBudget(elapsed_seconds = 7200), "exhausted")
})
