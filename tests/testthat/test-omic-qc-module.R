# KIT-104 W1: the QC tab through its public entry points, written before the shared QC shell
# (R/mod_omic_qc.R) and its sub-tabs (R/mod_omic_qc_tabs.R) exist. Each omic's tab is driven as the
# app shell drives it (`omicGuiComponentUi/Server("<prefix>_qc", "qc", ...)`, R/mod_omic_workflow.R),
# and the ITSD sub-tab through its exported servers, so the blocks hold through phases B to D.
#   (a) Parity: the finalise refuses a design that does not align with the assays, with no commit.
#   (b) Identity: the QC UI's HTML digests and its namespaced input and output ids, frozen today.
#   (c) F-24: a Detect duplicates or Analyze internal standards click made before the sub-modules
#       start is not run when they start.
#   (d) The ITSD capture: analyze_is's summary, visualisation tabs, CV-plot and intensity-plot data.
# (b) and (d) pass on the start commit (b423b764), where their values were frozen. The parity blocks
# for metab and trx and the F-24 blocks fail there (docs/qa/workflow-execution/kit/kit-104/w1/
# test-first.txt); W5 landed lipid's F-24 block and W6 the metab and trx parity blocks and metab's
# F-24 block, so none is skipped.

#' The rendered HTML's sha256, with shiny's random tabset ids masked (as the entry suite does).
omicQcHtml <- function(x) {
    html <- if (inherits(x, c("shiny.tag", "shiny.tag.list"))) htmltools::renderTags(x)$html else x$html
    html <- gsub("data-tabsetid=\"[0-9]+\"", "data-tabsetid=\"<ID>\"", as.character(html))
    gsub("tab-[0-9]+-([0-9]+)", "tab-<ID>-\\1", html)
}

omicQcDigest <- function(x) digest::digest(omicQcHtml(x), algo = "sha256", serialize = FALSE)

#' A value's sha256 (version-3 serialisation).
omicQcSha <- function(x) digest::digest(serialize(x, NULL, version = 3L), algo = "sha256", serialize = FALSE)

#' The element ids (`id` attributes) and test ids (`data-testid`, as `testid:<id>`) in rendered HTML,
#' sorted, without shiny's random tab ids.
omicQcIds <- function(x) {
    html <- omicQcHtml(x)
    found <- regmatches(html, gregexpr("[[:space:]](data-test)?id=\"[^\"]+\"", html))[[1L]]
    ids <- sub("^[[:space:]]data-testid=\"(.*)\"$", "testid:\\1", sub("^[[:space:]]id=\"(.*)\"$", "\\1", found))
    sort(unique(ids[!grepl("^tab-<ID>-", ids)]))
}

#' Capture shiny's notifications as "[type] text" in the caller's frame.
omicQcNotes <- function(env = parent.frame()) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    testthat::local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
        seen$notes <- c(seen$notes, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }, .package = "shiny", .env = env)
    seen
}

#' Drive one omic's QC tab as the app shell does, under module id `qc`. `trigger` is the design's QC
#' trigger; `body` runs in the test session (its first argument is the session).
omicQcTab <- function(workflow, component, omic_type, trigger = NULL, paths = list(), body) {
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(function(input, output, session) {
        session$userData$trigger <- shiny::reactiveVal(trigger)
        context <- list(workflow_data = workflow, experiment_paths = paths, omic_type = omic_type,
            experiment_label = "kit104", qc_trigger = session$userData$trigger, selected_tab = NULL)
        omicGuiComponentServer(component, "qc", context)
    }, body(session)))))
    invisible(NULL)
}

# ------------------------------------------------------------------------------------------------
# (a) Parity: the QC-finalise design-alignment check for every omic.
# ------------------------------------------------------------------------------------------------

#' A state whose design does not align with its assays: a design run absent from every assay
#' (`run`), or a run with a blank group (`group`). The object is built valid and then given the
#' design, as a state committed before a design check existed would hold it.
omicQcMisaligned <- function(object, how) {
    design <- object@design_matrix
    if (identical(how, "run")) {
        extra <- design[1L, , drop = FALSE]
        extra[[object@sample_id]] <- "S_absent"
        design <- rbind(design, extra)
    } else {
        design[[object@group_id]][[1L]] <- NA_character_
    }
    object@design_matrix <- design
    object
}

#' Start the QC tab and click Finalize QC once; returns the history and QC status before and after,
#' and the notifications.
omicQcFinalize <- function(workflow, component, omic_type, paths = list()) {
    notes <- omicQcNotes()
    out <- new.env(parent = emptyenv())
    out$before <- workflow$state_manager$getHistory()
    omicQcTab(workflow, component, omic_type, trigger = TRUE, paths = paths, body = function(session) {
        session$flushReact()
        session$setInputs(`qc-s4_finalize-finalize_qc` = 1L)
    })
    out$after <- workflow$state_manager$getHistory()
    out$qc <- workflow$tab_status$quality_control
    out$notes <- notes$notes
    as.list(out)
}

#' The finalise refused: nothing committed, QC not complete, an error notification holding `text`.
expectOmicQcRefused <- function(out, text, info) {
    expect_identical(out$after, out$before, info = info)
    expect_false(identical(out$qc, "complete"), info = info)
    errors <- out$notes[startsWith(out$notes, "[error] Error finalizing QC: ")]
    expect_length(errors, 1L)
    expect_true(any(grepl(text, errors, fixed = TRUE)), info = paste(info, paste(out$notes, collapse = " | ")))
}

#' A trx session owner on the trx_screen corpus whose design state holds `how`'s misaligned design.
omicQcTrxMisaligned <- function(how, env = parent.frame()) {
    fixture <- trxDeclarationsOwner(withr::local_tempdir(.local_envir = env), module_ci_trx_corpus())
    withr::defer(closeDIAStateAdapter(fixture$owner), envir = env)
    # The design state replaced in place, under its name and description.
    fixture$owner$saveState("trx_raw_data_s4", omicQcMisaligned(fixture$object, how), list(),
        "Initial TranscriptomicsAssayData S4 object created after design matrix")
    fixture
}

test_that("(a) lipid refuses a misaligned design at finalise, with no commit (the control)", {
    out <- omicQcFinalize(module_ci_lipid_qc_workflow(omicQcMisaligned(module_ci_lipid_qc_object(), "run")),
        "lipid_qc", "lipidomics")
    expectOmicQcRefused(out, "S_absent", "lipid, a design run absent from the assays")
    out <- omicQcFinalize(module_ci_lipid_qc_workflow(omicQcMisaligned(module_ci_lipid_qc_object(), "group")),
        "lipid_qc", "lipidomics")
    expectOmicQcRefused(out, paste("Current LipidomicsAssayData object has assay/design mismatch: Design",
        "matrix group column 'group' contains blank assignments."), "lipid, a blank group")
})

test_that("(a) metab refuses a design run absent from its assays today (its class's validity on the shell read)", {
    out <- omicQcFinalize(module_ci_metab_qc_workflow(omicQcMisaligned(module_ci_metab_qc_object(), "run")),
        "metab_qc", "metabolomics")
    expectOmicQcRefused(out, "S_absent", "metab, a design run absent from the assays")
})

test_that("(a) parity: metab refuses a design the alignment check refuses, with no commit (fails before)", {
    out <- omicQcFinalize(module_ci_metab_qc_workflow(omicQcMisaligned(module_ci_metab_qc_object(), "group")),
        "metab_qc", "metabolomics")
    # D11: the shared alignment with lipid's arguments and lipid's text for the class.
    expectOmicQcRefused(out, paste("Current MetaboliteAssayData object has assay/design mismatch: Design",
        "matrix group column 'group' contains blank assignments."), "metab, a blank group")
})

test_that("(a) parity: trx refuses a misaligned design at finalise, with no commit (fails before)", {
    skip_if_not(workflowOwnerCapable())
    skip_if_not_installed("edgeR")
    # D10: trx's own design check (validateTrxDesign), as every other trx step applies it.
    for (case in list(c(how = "run", text = "S_absent"), c(how = "group", text = "outside the declared levels"))) {
        fixture <- omicQcTrxMisaligned(case[["how"]])
        out <- omicQcFinalize(fixture$workflow, "trx_qc", "transcriptomics", fixture$paths)
        expectOmicQcRefused(out, case[["text"]], paste("trx", case[["how"]]))
    }
})

# ------------------------------------------------------------------------------------------------
# (b) Identity: the QC UI's digests and namespaced ids, frozen at b423b764.
# ------------------------------------------------------------------------------------------------

omicQcFrozen <- list(
    # The tab's own UI (the well panel, title and `dynamic_qc_tabs` slot).
    static = c(lipid_qc = "8961d9dce8d1c7b7cc0e567d3178fcf5ce894824e37c7e3975ab0d10de268878",
        metab_qc = "9801a70d58115260973dd9941aec919c4d32a9885bf3416127b823c27acc354d",
        trx_qc = "f7ee1eae6db0e81466ff1cfd1aefe1742af4385c90b5bc867d65d9b36d6e2d3d",
        prot_qc = "d63c209d7d09e7c861a4fb8d8ce551dba62601f4a438acf37cd6026a71c38e28"),
    # The tabset the server renders into `dynamic_qc_tabs`.
    tabs = c(lipid_qc = "d5724b7ded8cf9851f8f262f7f81803fbd6f672ff57008e8c7b1489fb2ab7ba5",
        metab_qc = "b63cab30bb386095e8b2e2a7dee290eded2495f5eb7a277032f560c9ccdf2a6d",
        trx_qc = "9d6924956f5e70f8076fdebdb31e9ea728df900fe36cfd527760ead918d7bfc1",
        prot_DIA = "002b2539c64a30b9dba67cc8d603a598206fcfad00c42642039bf1fa7044cb2f",
        prot_DIA_LIMPA = "ef0c045f70768b03f435c3c8ca55dd350b66a31e67b254fc933cf3cf173758b5",
        prot_TMT = "fcc59ebc749699d96e0392876691d3da5fe783c39b0914e9ac287bc73186f7d3",
        prot_LFQ = "fcc59ebc749699d96e0392876691d3da5fe783c39b0914e9ac287bc73186f7d3"),
    # The exported sub-tab and group UIs.
    subtabs = c(lipid_intensity = "cdd5302311c548ea1674f5f3df604e435870d99aacdf0396751d47b9412c5419",
        lipid_duplicates = "ab0f3fc1010376677b0084178f055b14438c37af1ddd73b645fc200965d0b825",
        lipid_itsd = "f214992241d651fb815bef303869e5c20e5fc2c0c7286a7f92f4302d66dafcbe",
        lipid_finalize = "e21323fc516e940694e91f0bcf69329625e0b6d90a405ae41826bc8cf61a1054",
        metab_intensity = "8af27b71c011f13ff33e8492d58ef8507ed7e76e5554c1ab487ba1c226fb2c98",
        metab_duplicates = "f766e81a7074d73bd21c1c8db63b4ed2d0e36abd40a8695c45e2eea19dc575cf",
        metab_itsd = "786265276d71b4fdbcf44f1e6ec9369d661bceeb481a2dcc107223af71aaf9a1",
        metab_finalize = "e21323fc516e940694e91f0bcf69329625e0b6d90a405ae41826bc8cf61a1054",
        trx_filter = "64458a3f45faf1d137e7266f450cfb17cc235ae8bc04f6f5a97bbd0a897ba162",
        trx_finalize = "bb5978191704f746d821032c894afaee695f9b6f0dc17319ba9ead806e0ae0b8",
        prot_peptide_DIA = "0754b131dad654c76c6f02cf822f907b42f3cdf209e37938257f9cb236220b67",
        prot_peptide_DIA_LIMPA = "0e907c496b7023350eacd8c7ff8c08010d70cc53ce332a3fee123117e65c6528",
        prot_protein_DIA = "7b635143802b8c07789b657c54b5a537f6e21730b9984e10a801191b931d9f52",
        prot_protein_DIA_LIMPA = "617f4c9e7a720ca96b12a5ee0332771fb0c2ee31d61885016c49f9e6c648f5eb",
        prot_protein_TMT = "0619511dd8d0da9ee2f8b023f49bc7e9ce4dfb6b7082e9eed5305a397a8ba6c9",
        prot_protein_LFQ = "0619511dd8d0da9ee2f8b023f49bc7e9ce4dfb6b7082e9eed5305a397a8ba6c9"),
    # The namespaced ids and test ids in the rendered tabset.
    ids = c(lipid_qc = "52668f1002eb15eb88a3d5117909127e96f1ea5e2647c550b96c5eedb69cd8d7",
        metab_qc = "a65948c8fe4cf9e5972b8338a3cb0b019526938bb73d502a11f6ad047973c3ff",
        trx_qc = "93a1d79e1e11b35de297b8f2267709e871709bd1cacb5081c727744a719208a2",
        prot_DIA = "718c024732fed8dc11f1cbd664da539ef3a8a82ad922103dda44001888f46aa1",
        prot_DIA_LIMPA = "b78791ba561d4441ad8cefc91881f4e1734c55b9a29bc95fb233c0b4a938a252",
        prot_TMT = "4952ee03c41953bfa964e2d2f9510234b55749561dc0fa41fe7d7facf584737c",
        prot_LFQ = "4952ee03c41953bfa964e2d2f9510234b55749561dc0fa41fe7d7facf584737c"),
    id_counts = c(lipid_qc = 29L, metab_qc = 29L, trx_qc = 14L, prot_DIA = 94L, prot_DIA_LIMPA = 47L,
        prot_TMT = 39L, prot_LFQ = 39L),
    # The outputs the tab's server and its started sub-modules register (proteomics' start on the
    # design's trigger only).
    outputs = list(
        lipid_qc = c("qc-duplicates-duplicate_summary", "qc-duplicates-duplicate_tables",
            "qc-duplicates-filter_plot", "qc-dynamic_qc_tabs", "qc-intensity-assay_results_tabs",
            "qc-intensity-filter_plot", "qc-itsd-cv_plot", "qc-itsd-intensity_plot", "qc-itsd-is_summary",
            "qc-itsd-is_viz_tabs", "qc-s4_finalize-assay_stats_table", "qc-s4_finalize-data_summary",
            "qc-s4_finalize-filter_plot", "qc-s4_finalize-state_history"),
        metab_qc = c("qc-duplicates-duplicate_summary", "qc-duplicates-duplicate_tables",
            "qc-duplicates-filter_plot", "qc-dynamic_qc_tabs", "qc-intensity-assay_results_tabs",
            "qc-intensity-filter_plot", "qc-itsd-cv_plot", "qc-itsd-intensity_plot", "qc-itsd-is_summary",
            "qc-itsd-is_viz_tabs", "qc-s4_finalize-assay_stats_table", "qc-s4_finalize-data_summary",
            "qc-s4_finalize-filter_plot", "qc-s4_finalize-state_history"),
        trx_qc = c("qc-dynamic_qc_tabs", "qc-filter-filter_plot", "qc-filter-raw_pca", "qc-filter-raw_rle",
            "qc-filter-sample_metrics", "qc-s4_finalize-data_summary", "qc-s4_finalize-filter_plot",
            "qc-s4_finalize-state_history"),
        prot_DIA = "qc-dynamic_qc_tabs", prot_DIA_LIMPA = "qc-dynamic_qc_tabs", prot_TMT = "qc-dynamic_qc_tabs",
        prot_LFQ = "qc-dynamic_qc_tabs"))

#' The QC tab of one omic as the app shell renders it: the tab's own UI, the tabset its server
#' renders, its ids, and the outputs its server and the started sub-modules register.
omicQcRendered <- function(component, omic_type, workflow, paths = list()) {
    out <- new.env(parent = emptyenv())
    omicQcTab(workflow, component, omic_type, paths = paths, body = function(session) {
        session$flushReact()
        out$tabs <- session$getOutput("qc-dynamic_qc_tabs")
        out$outputs <- sort(names(session$.__enclos_env__$private$outs))
    })
    list(tabs = omicQcDigest(out$tabs), ids = omicQcIds(out$tabs), outputs = out$outputs)
}

#' Every identity value, by the names `omicQcFrozen` uses.
omicQcIdentity <- function(env = parent.frame()) {
    components <- c(lipid_qc = "lipid_qc", metab_qc = "metab_qc", trx_qc = "trx_qc", prot_qc = "prot_qc")
    static <- vapply(components, function(component) omicQcDigest(omicGuiComponentUi(component, "qc")), "")
    trx <- trxDeclarationsOwner(withr::local_tempdir(.local_envir = env), module_ci_trx_corpus())
    withr::defer(closeDIAStateAdapter(trx$owner), envir = env)
    rendered <- list(
        lipid_qc = omicQcRendered("lipid_qc", "lipidomics", module_ci_lipid_qc_workflow()),
        metab_qc = omicQcRendered("metab_qc", "metabolomics", module_ci_metab_qc_workflow()),
        trx_qc = omicQcRendered("trx_qc", "transcriptomics", trx$workflow, trx$paths))
    for (route in c("DIA", "DIA_LIMPA", "TMT", "LFQ")) {
        # The route as the design leaves it on the owner (test-perf107-audit.R's construction).
        workflow <- shiny::reactiveValues(state_manager = list(workflow_type = route))
        rendered[[paste0("prot_", route)]] <- omicQcRendered("prot_qc", "proteomics", workflow)
    }
    trx_profile <- trxOmicProfile()
    subtabs <- c(
        lipid_intensity = omicQcDigest(mod_lipid_qc_intensity_ui("x")),
        lipid_duplicates = omicQcDigest(mod_lipid_qc_duplicates_ui("x")),
        lipid_itsd = omicQcDigest(mod_lipid_qc_itsd_ui("x")),
        lipid_finalize = omicQcDigest(mod_lipid_qc_s4_ui("x")),
        metab_intensity = omicQcDigest(mod_metab_qc_intensity_ui("x")),
        metab_duplicates = omicQcDigest(mod_metab_qc_duplicates_ui("x")),
        metab_itsd = omicQcDigest(mod_metab_qc_itsd_ui("x")),
        metab_finalize = omicQcDigest(mod_metab_qc_s4_ui("x")),
        trx_filter = omicQcDigest(mod_trx_qc_filter_ui("x")),
        trx_finalize = omicQcDigest(mod_lipid_qc_s4_ui("x", profile = trx_profile)),
        prot_peptide_DIA = omicQcDigest(mod_prot_qc_peptide_ui("x", workflowType = "DIA")),
        prot_peptide_DIA_LIMPA = omicQcDigest(mod_prot_qc_peptide_ui("x", workflowType = "DIA_LIMPA")),
        prot_protein_DIA = omicQcDigest(mod_prot_qc_protein_ui("x", workflow_type = "DIA")),
        prot_protein_DIA_LIMPA = omicQcDigest(mod_prot_qc_protein_ui("x", workflow_type = "DIA_LIMPA")),
        prot_protein_TMT = omicQcDigest(mod_prot_qc_protein_ui("x", workflow_type = "TMT")),
        prot_protein_LFQ = omicQcDigest(mod_prot_qc_protein_ui("x", workflow_type = "LFQ")))
    # The ids as their sha256 and count (the lists are in docs/qa/workflow-execution/kit/kit-104/w1/
    # identity-ids.txt); the outputs as their names.
    list(static = static, tabs = vapply(rendered, `[[`, "", "tabs"), subtabs = subtabs,
        ids = vapply(rendered, function(r) omicQcSha(r$ids), ""),
        id_counts = vapply(rendered, function(r) length(r$ids), 0L),
        outputs = lapply(rendered, `[[`, "outputs"))
}

test_that("(b) each omic's QC UI, sub-tab UIs and namespaced ids are today's (lipid, metab, trx, prot by route)", {
    skip_if_not(workflowOwnerCapable())
    skip_if_not_installed("edgeR")
    actual <- omicQcIdentity()
    for (part in names(omicQcFrozen)) {
        expect_identical(actual[[part]], omicQcFrozen[[part]], info = part)
    }
})

# ------------------------------------------------------------------------------------------------
# (c) F-24: the two display-only actions do not run a click made before the sub-modules start.
# ------------------------------------------------------------------------------------------------

#' KIT-101 F-17's construction (test-module-ci-{lipid,metab}-qc.R): the QC tab as the design leaves
#' it; `button` is clicked (1) before the design's trigger starts the sub-modules (before the design
#' the tab has no state), then clicked again (2) after the start. Returns the notifications at the
#' start and after the second click.
omicQcClickBeforeStart <- function(workflow, component, omic_type, state, button) {
    owner <- workflow$state_manager
    selectWorkflowStateReference(owner, state)
    workflow$state_manager <- NULL
    notes <- omicQcNotes()
    out <- new.env(parent = emptyenv())
    omicQcTab(workflow, component, omic_type, body = function(session) {
        click <- function(value) do.call(session$setInputs, stats::setNames(list(value), paste0("qc-", button)))
        # The pattern box as the browser sends it, empty: the state's own pattern is used.
        session$setInputs(`qc-itsd-is_pattern` = "")
        click(1L)
        workflow$state_manager <- owner
        session$userData$trigger(TRUE)
        session$flushReact()
        out$at_start <- notes$notes
        click(2L)
        out$after_click <- notes$notes
    })
    as.list(out)
}

# What one click after the start shows (lipid and metab alike on the MCI fixtures). Today the click
# made before the start also runs when the sub-modules start, so these show twice.
omicQcDisplayOnly <- list(
    `duplicates-detect_duplicates` = "[warning] Detection complete: 4 duplicate IDs found",
    `itsd-analyze_is` = c("[default] Analyzing internal standards...", "[message] Found 4 internal standards"))

test_that("(c) F-24: a lipid Detect duplicates or Analyze internal standards click before the start is not run", {
    for (button in names(omicQcDisplayOnly)) {
        out <- omicQcClickBeforeStart(module_ci_lipid_qc_workflow(), "lipid_qc", "lipidomics", "lipid_raw_data_s4",
            button)
        expect_identical(out$at_start, character(), info = button)
        # A click after the start runs as before.
        expect_identical(out$after_click, omicQcDisplayOnly[[button]], info = button)
    }
})

test_that("(c) F-24: a metab Detect duplicates or Analyze internal standards click before the start is not run", {
    for (button in names(omicQcDisplayOnly)) {
        out <- omicQcClickBeforeStart(module_ci_metab_qc_workflow(), "metab_qc", "metabolomics", "metab_raw_data_s4",
            button)
        expect_identical(out$at_start, character(), info = button)
        expect_identical(out$after_click, omicQcDisplayOnly[[button]], info = button)
    }
})

# ------------------------------------------------------------------------------------------------
# (d) The ITSD capture: the success path, which no lane reaches.
# ------------------------------------------------------------------------------------------------

#' Analyze internal standards through an omic's exported ITSD server (pattern empty: the state's
#' own) and capture what it shows: the summary and visualisation-tab HTML, the results text, the
#' data the CV and intensity plots draw, and the notifications.
omicQcItsdCapture <- function(server, workflow, omic_type) {
    notes <- omicQcNotes()
    drawn <- new.env(parent = emptyenv())
    drawn$data <- list()
    ggplot <- ggplot2::ggplot
    testthat::local_mocked_bindings(ggplot = function(data = NULL, ...) {
        drawn$data <- c(drawn$data, list(data))
        ggplot(data, ...)
    }, .package = "ggplot2")
    out <- new.env(parent = emptyenv())
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(function(input, output, session) {
        server("itsd", workflow, omic_type, "kit104")
    }, {
        session$setInputs(`itsd-is_pattern` = "", `itsd-analyze_is` = 1L)
        out$summary <- omicQcDigest(output[["itsd-is_summary"]])
        out$tabs <- omicQcDigest(output[["itsd-is_viz_tabs"]])
        out$results <- omicQcSha(output[["itsd-is_results"]])
        # The plots are drawn in the flush; each one's data is told apart by its columns.
        expect_true(all(c("src", "coordmap") %in% names(output[["itsd-cv_plot"]])))
        expect_true(all(c("src", "coordmap") %in% names(output[["itsd-intensity_plot"]])))
    }))))
    drawn_by <- function(column) Filter(function(data) column %in% names(data), drawn$data)
    out$cv_data <- omicQcSha(drawn_by("cv_status"))
    out$intensity_data <- omicQcSha(drawn_by("IS_ID"))
    out$plots_drawn <- length(drawn$data)
    out$notes <- notes$notes
    as.list(out)[c("summary", "tabs", "results", "cv_data", "intensity_data", "plots_drawn", "notes")]
}

omicQcItsdFrozen <- list(
    lipid = list(summary = "18820161286938d7a3f935aa2cd889f8806b0913894b8584bec1f4cded4fef86",
        tabs = "ef19f6e589787c704250c9e62bb4dc2e66a1844d6b5b04e472500379acecfcdb",
        results = "6af22df84db1bd93d64e2547be3e6249b3406a72452a8e476ffea1b7057d1a2a",
        cv_data = "cbe87ff4ebeeab996f60ba013589bb95a370d3f1b9245325a45f8ed581f7532e",
        intensity_data = "941d0808c1d3b70597f5a57cfee080e15a32d1d996f243273d4e8328481e8d42",
        plots_drawn = 2L,
        notes = c("[default] Analyzing internal standards...", "[message] Found 4 internal standards")),
    # The metab fixture mirrors the lipid one, so its summary, tabs and text are lipid's; its plot
    # data differ.
    metab = list(summary = "18820161286938d7a3f935aa2cd889f8806b0913894b8584bec1f4cded4fef86",
        tabs = "ef19f6e589787c704250c9e62bb4dc2e66a1844d6b5b04e472500379acecfcdb",
        results = "6af22df84db1bd93d64e2547be3e6249b3406a72452a8e476ffea1b7057d1a2a",
        cv_data = "69f605c90da9573afedaae8b6e2eed02d6c544cc75f0798fe923a628bd56cbff",
        intensity_data = "a42db1e2454694f6b7eb453f91a5ec2533f740f7da872905194b3b083dc727bb",
        plots_drawn = 2L,
        notes = c("[default] Analyzing internal standards...", "[message] Found 4 internal standards")))

test_that("(d) analyze_is through the lipid and metab ITSD servers shows today's summary, tabs, text and plot data", {
    lipid <- omicQcItsdCapture(mod_lipid_qc_itsd_server, module_ci_lipid_qc_workflow(), "lipidomics")
    expect_identical(lipid, omicQcItsdFrozen$lipid)
    metab <- omicQcItsdCapture(mod_metab_qc_itsd_server, module_ci_metab_qc_workflow(), "metabolomics")
    expect_identical(metab, omicQcItsdFrozen$metab)
})

# ------------------------------------------------------------------------------------------------
# KIT-104 W5: the shared shell (R/mod_omic_qc.R) and sub-tabs (R/mod_omic_qc_tabs.R), lipid first.
# Lipid runs on them through its wrappers ((a) to (d) above, MCI-020). These blocks hold what the
# later consumers need of them, against today's modules: proteomics' tabsets, start and completion
# by route (mod_prot_qc_server), and trx's one-table summary and text (registerTrxQcSummaryOutput,
# trxQcFinalizeText).
# ------------------------------------------------------------------------------------------------

#' The tabset digest, ids and registered outputs of `mod_omic_qc_server("qc", ...)` under `config`.
omicQcShellRendered <- function(config, workflow, omic_type) {
    held <- new.env(parent = emptyenv())
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(function(input, output, session) {
        session$userData$trigger <- shiny::reactiveVal(NULL)
        mod_omic_qc_server("qc", workflow, list(), omic_type, "kit104", session$userData$trigger, config = config)
    }, {
        session$flushReact()
        held$tabs <- session$getOutput("qc-dynamic_qc_tabs")
        held$outputs <- sort(names(session$.__enclos_env__$private$outs))
    }))))
    list(tabs = omicQcDigest(held$tabs), ids = omicQcSha(omicQcIds(held$tabs)), outputs = held$outputs)
}

test_that("W5 the shell renders proteomics' tabsets by route as mod_prot_qc_server does, before the trigger", {
    config <- omicGuiConfig("proteomics")
    expect_identical(omicQcDigest(mod_omic_qc_ui("qc", config)), omicQcFrozen$static[["prot_qc"]])
    for (route in c("DIA", "DIA_LIMPA", "TMT", "LFQ")) {
        key <- paste0("prot_", route)
        rendered <- omicQcShellRendered(config, shiny::reactiveValues(state_manager = list(workflow_type = route)),
            "proteomics")
        expect_identical(rendered, list(tabs = omicQcFrozen$tabs[[key]], ids = omicQcFrozen$ids[[key]],
            outputs = omicQcFrozen$outputs[[key]]), info = route)
    }
})

test_that("W5 the shell starts proteomics' groups on the trigger by route, keeps `qc_workflow_type` and completes", {
    started <- character()
    testthat::local_mocked_bindings(
        mod_prot_qc_peptide_server = function(id, ...) started <<- c(started, id),
        mod_prot_qc_protein_server = function(id, ...) started <<- c(started, id),
        .package = "MultiScholaR")
    config <- omicGuiConfig("proteomics")
    for (case in list(list(route = "DIA", groups = c("peptide_qc", "protein_qc")),
        list(route = "TMT", groups = "protein_qc"))) {
        started <- character()
        bus <- shiny::reactiveValues(state_manager = list(workflow_type = case$route, states = list()),
            tab_status = list(quality_control = "pending", normalization = "disabled"))
        suppressMessages(shiny::testServer(function(input, output, session) {
            # NULL until the design saves (a FALSE at the start ends the once-only start observer, as in
            # mod_prot_qc.R today).
            session$userData$go <- shiny::reactiveVal(NULL)
            mod_omic_qc_server("qc", bus, list(), "proteomics", "kit104", session$userData$go, config = config)
        }, {
            session$flushReact()
            expect_identical(started, character(), info = case$route)
            session$userData$go(TRUE)
            session$flushReact()
            expect_identical(started, case$groups, info = case$route)
            # The final protein state completes QC and opens normalisation (completeProtQcWorkflowStatus).
            bus$state_manager <- list(workflow_type = case$route, states = list(protein_replicate_filtered = TRUE))
            session$flushReact()
            expect_identical(shiny::isolate(bus$tab_status), list(quality_control = "complete",
                normalization = "pending"), info = case$route)
        }))
    }
    # test-perf107-audit.R's reading of the shell's local route, under an owner that changes its route in place.
    native <- WorkflowState$new(audit_enabled = FALSE)
    bus <- shiny::reactiveValues(state_manager = native)
    go <- shiny::reactiveVal(FALSE)
    suppressMessages(shiny::testServer(mod_omic_qc_server, args = list(workflow_data = bus, experiment_paths = list(),
        omic_type = "proteomics", experiment_label = "kit104", qc_trigger = go, config = config), {
        session$flushReact()
        expect_identical(qc_workflow_type(), "LFQ")
        native$setWorkflowType("DIA")
        go(TRUE)
        session$flushReact()
        expect_identical(qc_workflow_type(), "DIA")
    }))
})

test_that("W5 the shared finalise's one-table summary and text are trx's today (registerTrxQcSummaryOutput, trxQcFinalizeText)", {
    skip_if_not(workflowOwnerCapable())
    skip_if_not_installed("edgeR")
    fixture <- trxDeclarationsOwner(withr::local_tempdir(), module_ci_trx_corpus())
    withr::defer(closeDIAStateAdapter(fixture$owner))
    config <- omicGuiConfig("transcriptomics")
    held <- new.env(parent = emptyenv())
    suppressWarnings(suppressMessages(shiny::testServer(function(input, output, session) {
        registerOmicQcFinalizeSummary(output, fixture$workflow, config)
    }, {
        held$shared <- output$data_summary
    })))
    expect_false(is.null(held$shared$html))
    # KIT-104 W6 deleted trx's copies: the summary's HTML and the text are the sha256 of what
    # registerTrxQcSummaryOutput and trxQcFinalizeText gave on this fixture at 827f2328 (w6/test-first.txt).
    sha <- function(text) digest::digest(as.character(text), algo = "sha256", serialize = FALSE)
    expect_identical(sha(held$shared$html), "01b96be2ab75dabe19023041fc713e37f9e4b979118f607310fd5cdbd8045201")
    meta <- readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, fixture$owner)
    history <- fixture$owner$getHistory()
    expect_identical(sha(omicQcFinalizeText(sum(meta$assay_rows), NULL, history, config)),
        "068f1028210d199b0ee99687b5ca511234bdf96ead9ad8d3dc6553ce4c2cb179")
})

test_that("W5 lipid's QC copies are deleted; its exports are wrappers over the shared shell and sub-tabs", {
    ns <- asNamespace("MultiScholaR")
    deleted <- c("registerLipidIntensityApplyFilterObserver", "registerLipidIntensityAssayResultsOutput",
        "registerLipidIntensityFilterPlotOutput", "handleLipidDuplicateDetection", "handleLipidDuplicateResolution",
        "handleLipidDuplicateRevert", "registerLipidDuplicateServerBindings", "runLipidDuplicateModuleServerShell",
        "registerLipidItsdSummaryOutput", "registerLipidItsdCvPlotOutput", "buildLipidQcS4DataSummaryUi",
        "buildLipidQcS4AssayStatsTable", "buildLipidQcS4FinalizeResultsText", "validateLipidQcS4FinalizeState",
        "updateLipidQcS4TrackingPlot", "getLipidQcS4FinalizeState", "saveLipidQcS4CompletedState")
    expect_identical(deleted[vapply(deleted, exists, logical(1), envir = ns, inherits = FALSE)], character())
    # qcClickAfterStart moved into the shell's file with its name.
    expect_match(utils::getSrcFilename(qcClickAfterStart), "mod_omic_qc.R", fixed = TRUE)
    for (name in c("mod_omic_qc_ui", "mod_omic_qc_server", "omicQcComponents", "registerOmicQcRevertObserver")) {
        expect_match(utils::getSrcFilename(get(name, envir = ns)), "mod_omic_qc.R", fixed = TRUE, info = name)
    }
    for (part in c("intensity", "duplicates", "itsd", "finalize")) {
        for (side in c("ui", "server")) {
            name <- sprintf("mod_omic_qc_%s_%s", part, side)
            expect_match(utils::getSrcFilename(get(name, envir = ns)), "mod_omic_qc_tabs.R", fixed = TRUE, info = name)
        }
    }
    # The registry names a module pair for every component the four omics' qc groups name.
    components <- omicQcComponents()
    for (omic in c("lipidomics", "metabolomics", "transcriptomics", "proteomics")) {
        for (group in omicGuiConfig(omic)$qc$groups) {
            ids <- c(group$component, vapply(group$tabs, `[[`, "", "component"))
            expect_true(all(ids %in% names(components)), info = omic)
            for (id in intersect(ids, names(components))) {
                expect_true(exists(components[[id]]$ui, envir = ns) && exists(components[[id]]$server, envir = ns),
                    info = id)
            }
        }
    }
})

test_that("W5 the shell's and sub-tabs' module frames hold no name a test's own variables would be masked by", {
    # The variables the QC tests (this file, MCI-020, MCI-028, test-perf107-audit.R) use inside testServer.
    test_names <- c("out", "seen", "notes", "workflow", "fixture", "owner", "click", "button", "renders", "first",
        "native", "trigger", "go", "bus", "started", "case", "counts", "captured", "probe", "design", "after_click",
        "config", "route", "drawn", "ggplot", "held", "legacy", "manager", "lipid", "metab")
    held <- new.env(parent = emptyenv())
    held$frames <- list()
    capture <- function(server, args) {
        suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(server, args = args, {
            held$frames <- c(held$frames, list(ls(session$env, all.names = TRUE)))
        }))))
    }
    lipid <- omicGuiConfig("lipidomics")
    capture(mod_omic_qc_server, list(workflow_data = module_ci_lipid_qc_workflow(), experiment_paths = list(),
        omic_type = "lipidomics", experiment_label = "kit104", qc_trigger = function() NULL, config = lipid))
    # The shell keeps the local a test reads by name (test-perf107-audit.R).
    expect_true("qc_workflow_type" %in% held$frames[[1L]])
    for (part in c("intensity", "duplicates", "itsd", "finalize")) {
        capture(get(sprintf("mod_omic_qc_%s_server", part)), list(workflow_data = module_ci_lipid_qc_workflow(),
            omic_type = "lipidomics", experiment_label = "kit104", config = lipid))
    }
    expect_length(held$frames, 5L)
    for (frame in held$frames) expect_identical(intersect(frame, test_names), character())
})

# ------------------------------------------------------------------------------------------------
# KIT-104 W7: proteomics' QC tab on the shared shell. mod_prot_qc_ui/server are wrappers over
# mod_omic_qc_*(..., omicGuiConfig("proteomics")); the peptide and protein orchestrators stay as the
# group plugins; getProtQc*ModuleSpecs return today's lists, built from the configuration (D9).
# P-16 (the lead's decision, 2026-11-08): DIA_LIMPA's QC completes at peptide_count_filtered.
# ------------------------------------------------------------------------------------------------

test_that("W7 mod_prot_qc_ui/server are wrappers over the shared shell, without the DEBUG66 messages", {
    expect_true("mod_omic_qc_ui" %in% all.names(body(mod_prot_qc_ui)))
    expect_true("mod_omic_qc_server" %in% all.names(body(mod_prot_qc_server)))
    expect_false(any(grepl("DEBUG66", deparse(mod_prot_qc_server), fixed = TRUE)))
    expect_identical(names(formals(mod_prot_qc_ui)), "id")
    expect_identical(names(formals(mod_prot_qc_server)),
        c("id", "workflow_data", "experiment_paths", "omic_type", "experiment_label", "qc_trigger"))
    # The static UI is today's (identity (b)).
    expect_identical(omicQcDigest(mod_prot_qc_ui("qc")), omicQcFrozen$static[["prot_qc"]])
})

test_that("W7 P-16: QC completes at the route's final state through mod_prot_qc_server, DIA_LIMPA at peptide_count_filtered", {
    testthat::local_mocked_bindings(
        mod_prot_qc_peptide_server = function(...) NULL,
        mod_prot_qc_protein_server = function(...) NULL,
        .package = "MultiScholaR")
    completes <- function(route, state) {
        bus <- shiny::reactiveValues(state_manager = list(workflow_type = route, states = list()),
            tab_status = list(quality_control = "pending", normalization = "disabled"))
        suppressMessages(shiny::testServer(mod_prot_qc_server, args = list(workflow_data = bus,
            experiment_paths = list(), omic_type = "proteomics", experiment_label = "kit104",
            qc_trigger = function() TRUE), {
            session$flushReact()
            bus$state_manager <- list(workflow_type = route, states = stats::setNames(list(TRUE), state))
            session$flushReact()
        }))
        shiny::isolate(bus$tab_status)
    }
    complete <- list(quality_control = "complete", normalization = "pending")
    pending <- list(quality_control = "pending", normalization = "disabled")
    # P-16: DIA_LIMPA's last QC state (fails before: its QC status stays pending).
    expect_identical(completes("DIA_LIMPA", "peptide_count_filtered"), complete)
    # The other routes complete as today, and not on another route's final state.
    for (route in c("DIA", "TMT", "LFQ")) {
        expect_identical(completes(route, "protein_replicate_filtered"), complete, info = route)
        expect_identical(completes(route, "peptide_count_filtered"), pending, info = route)
    }
})

test_that("W7 getProtQc*ModuleSpecs return today's lists for every route, built from the configuration", {
    # Each list's sha256 (version-3 serialisation) at the start commit 827f2328.
    peptide <- c(default = omicQcSha(getProtQcPeptideModuleSpecs()),
        DIA = omicQcSha(getProtQcPeptideModuleSpecs("DIA")),
        DIA_LIMPA = omicQcSha(getProtQcPeptideModuleSpecs("DIA_LIMPA")),
        TMT = omicQcSha(getProtQcPeptideModuleSpecs("TMT")), LFQ = omicQcSha(getProtQcPeptideModuleSpecs("LFQ")),
        null = omicQcSha(getProtQcPeptideModuleSpecs(NULL)))
    dia_peptide <- "f12eb2444c894fab4390fe01098d8abc03f7722853e9ee2af06ff474393c0751"
    expect_identical(peptide, c(default = dia_peptide, DIA = dia_peptide,
        DIA_LIMPA = "ce13cbbfe8d0c1375ed775321f18b32b4e6528855c2e4512b9fb7edae3e7eb30", TMT = dia_peptide,
        LFQ = dia_peptide, null = dia_peptide))
    protein <- c(default = omicQcSha(getProtQcProteinModuleSpecs()),
        DIA = omicQcSha(getProtQcProteinModuleSpecs("DIA")),
        DIA_LIMPA = omicQcSha(getProtQcProteinModuleSpecs("DIA_LIMPA")),
        TMT = omicQcSha(getProtQcProteinModuleSpecs("TMT")), LFQ = omicQcSha(getProtQcProteinModuleSpecs("LFQ")),
        rollup = omicQcSha(getProtQcProteinRollupModuleSpec()), common = omicQcSha(getProtQcProteinCommonModuleSpecs()))
    dia_protein <- "c835f2cf6cfe081ed71e515bdb1d302506d1eeb57c2e77a9ce25b428f4bab530"
    common <- "e8163cffe36c5926145d1d189152a6dbf06a414c5d5def6742d17d790f12cef1"
    expect_identical(protein, c(default = dia_protein, DIA = dia_protein,
        DIA_LIMPA = "519e3103e07568bfdc6f5c51c1dc2497e49eb41c0600826bdcba067cb3b0968d", TMT = common, LFQ = common,
        rollup = "c4d065800d06aaee39fcc99c54c468f2122cd080c8a5cedc3d6c048678f05cd5", common = common))
    # Built from the configuration: a title changed there is the spec's (fails before: the tables were literals).
    original <- omicGuiProteomicsQc
    testthat::local_mocked_bindings(omicGuiProteomicsQc = function() {
        qc <- original()
        qc$groups[[1L]]$tabs[[1L]]$title <- "Q-Value Filter (configured)"
        qc$groups[[2L]]$tabs[[2L]]$title <- "Accession Cleanup (configured)"
        qc
    }, .package = "MultiScholaR")
    expect_identical(getProtQcPeptideModuleSpecs("DIA_LIMPA")[[1L]]$title, "Q-Value Filter (configured)")
    expect_identical(getProtQcProteinCommonModuleSpecs()[[1L]]$title, "Accession Cleanup (configured)")
    expect_identical(getProtQcProteinModuleSpecs("TMT")[[1L]]$title, "Accession Cleanup (configured)")
})

# ------------------------------------------------------------------------------------------------
# KIT-104 W6: metab and trx on the shared shell and sub-tabs; trx's finalise on a structural shell
# (F104-04); the metab and trx QC copies deleted.
# ------------------------------------------------------------------------------------------------

test_that("W6 F104-04: trx's shared finalise reads a structural shell (no validity), checks it and commits trx_qc_complete", {
    skip_if_not(workflowOwnerCapable())
    skip_if_not_installed("edgeR")
    fixture <- trxDeclarationsOwner(withr::local_tempdir(), module_ci_trx_corpus())
    withr::defer(closeDIAStateAdapter(fixture$owner))
    owner <- fixture$owner
    config <- omicGuiConfig("transcriptomics")
    binding <- transcriptomicsDeclarationBinding
    meta <- readWorkflowOmicsMetadata(binding, owner)
    # The generic shell read refuses every trx state: the class's validity on the zero-row prototype.
    expect_error(readWorkflowOmicsShell(binding, owner, measures = "features"), class = "trx_counts_invalid")
    # The structural shell: the stored slots and design, the count table as its zero-row prototype, the stored rows.
    expect_true(is.function(config$qc$finalize_shell))
    shell <- config$qc$finalize_shell(binding, owner)
    expect_s4_class(shell, "TranscriptomicsAssayData")
    expect_identical(attr(shell, "assay_rows"), meta$assay_rows)
    expect_identical(names(shell@count_data), meta$assay_names)
    expect_identical(nrow(shell@count_data[[1L]]), 0L)
    expect_identical(names(shell@count_data[[1L]]), names(fixture$object@count_data[[1L]]))
    expect_identical(shell@design_matrix, meta$design_matrix)
    expect_identical(shell@args, fixture$object@args)
    expect_identical(validateWorkflowOmicsQcFinalize(shell, config$qc, binding), shell)
    # The shared finalise sub-tab under trx's configuration commits and reports as trx's did.
    notes <- omicQcNotes()
    held <- new.env(parent = emptyenv())
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(function(input, output, session) {
        mod_omic_qc_finalize_server("s4_finalize", fixture$workflow, "transcriptomics", "kit104", config)
    }, {
        session$setInputs(`s4_finalize-finalize_qc` = 1L)
        session$flushReact()
        held$text <- output$`s4_finalize-finalize_results`
    }))))
    history <- owner$getHistory()
    expect_identical(history, c("initial", "trx_raw_data_s4", "trx_qc_complete"))
    expect_identical(fixture$workflow$tab_status$quality_control, "complete")
    expect_identical(held$text, omicQcFinalizeText(sum(meta$assay_rows), NULL, history, config))
    expect_identical(notes$notes, "[message] QC complete! Proceed to Normalization.")
})

test_that("W6 metab's and trx's QC copies are deleted; their exports are wrappers; the trx profile has no QC hook", {
    ns <- asNamespace("MultiScholaR")
    deleted <- c("initializeMetabQcSubmodules", "buildMetabIntensityAssayTabsUi", "buildMetabIntensityFilterStats",
        "reportMetabIntensityFilterSuccess", "buildMetabDuplicateSummaryUi", "detectMetabDuplicateFeatures",
        "runMetabDuplicateResolutionObserver", "registerMetabDuplicateTableRenderers", "metabQcItsdColumns",
        "readMetabQcItsdState", "analyzeMetabQcItsdData", "runMetabQcItsdServerBody", "buildMetabQcS4DataSummaryUi",
        "validateMetabQcS4FinalizeState", "runMetabQcS4FinalizeWorkflow", "runMetabQcS4ServerBody",
        "trxQcTabsUi", "initializeTrxQcSubmodules", "registerTrxQcSummaryOutput", "trxQcFinalizeText",
        "trxQcFinalizeCollaborators", "runLipidQcS4FinalizeWorkflow", "reportLipidQcS4FinalizeSuccess",
        "registerLipidQcS4ProfileFinalize", "registerLipidIntensityRevertObserver")
    expect_identical(deleted[vapply(deleted, exists, logical(1), envir = ns, inherits = FALSE)], character())
    # trx's four QC figure functions are defined once, in the core.
    for (name in c("trxQcRawStats", "trxQcFilterText", "trxQcMetricsTable", "trxQcFilteringFigure")) {
        expect_match(utils::getSrcFilename(get(name, envir = ns)), "func_trx_qc_plots.R", fixed = TRUE, info = name)
    }
    profile <- trxOmicProfile()
    # KIT-104 W8: `class` and `read_state` too, once test-omic-gui-config.R pinned their values.
    expect_identical(intersect(c("qc_title", "qc_tabs", "qc_init", "register_qc_summary", "finalize", "read_metadata",
        "class", "read_state"), names(profile)), character())
    expect_identical(profile$omic_type, "transcriptomics")
    # Every finalise checks its alignment; trx's is its own design check (D10).
    for (omic in c("lipidomics", "metabolomics", "transcriptomics")) {
        expect_true(is.function(omicGuiConfig(omic)$qc$finalize_alignment), info = omic)
    }
    expect_identical(omicGuiConfig("transcriptomics")$qc$finalize_alignment, trxDesignAlignment)
})
