# DIA-002 and WFX-001 fixture builders. No objects are allocated when this helper loads.
# Existing S4 slot access is intentional in fixture construction/verification;
# these are the package's real classes, not a substitute data model.

#' Read an explicitly named fixture recipe.
#' @param recipe_id Unique recipe identifier.
#' @param repo_root Repository root; NULL resolves through here.
#' @return A recipe list.
#' @noRd
readWorkflowCheckpointRecipe <- function(recipe_id, repo_root = NULL) {
    if (is.null(repo_root)) repo_root <- here::here()
    path <- file.path(repo_root, "tests", "testdata", "workflow-checkpoints", "scenarios.json")
    manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
    matches <- Filter(\(x) identical(x$id, recipe_id), manifest$recipes)
    if (length(matches) != 1L) stop("fixture recipe must resolve uniquely", call. = FALSE)
    matches[[1L]]
}

#' Validate a large-fixture output root before allocation or writing.
#' @param output_root Existing explicit directory outside package source.
#' @param repo_root Package source root.
#' @return Canonical output directory, invisibly.
#' @noRd
validateWorkflowFixtureRoot <- function(output_root, repo_root) {
    if (!is.character(output_root) || length(output_root) != 1L ||
        is.na(output_root) || !dir.exists(output_root)) {
        stop("an existing explicit fixture output directory is required", call. = FALSE)
    }
    output_root <- normalizePath(output_root, mustWork = TRUE)
    repo_root <- normalizePath(repo_root, mustWork = TRUE)
    broad <- unique(c("/", path.expand("~"), dirname(repo_root)))
    if (output_root %in% broad || identical(output_root, repo_root) ||
        startsWith(output_root, paste0(repo_root, "/"))) {
        stop("large fixtures must be outside package source in a dedicated directory", call. = FALSE)
    }
    if (!identical(Sys.info()[["sysname"]], "Linux")) {
        stop("physical filesystem validation is currently Linux-only", call. = FALSE)
    }
    filesystem <- system2("stat", c("-f", "-c", "%T", "--", shQuote(output_root)),
        stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(filesystem, "status")) || length(filesystem) != 1L ||
        !filesystem %in% c("btrfs", "ext2/ext3", "xfs", "zfs")) {
        stop("fixture root is not a verified disk-backed filesystem", call. = FALSE)
    }
    invisible(output_root)
}

#' Check dimensions and enforce the small-test allocation boundary.
#' @param recipe Fixture recipe.
#' @param allow_large Explicit large-data opt-in.
#' @param output_root Existing dedicated disk-backed directory.
#' @param repo_root Repository root; NULL resolves through here.
#' @return The validated recipe.
#' @noRd
validateWorkflowCheckpointRecipe <- function(recipe, allow_large = FALSE,
                                             output_root = NULL, repo_root = NULL) {
    fields <- c("n_proteins", "n_samples", "peptides_per_protein", "forms_per_peptide",
        "history_depth", "annotation_chars", "seed")
    for (field in fields) {
        value <- recipe[[field]]
        if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
            !is.finite(value) || value < 1 || value != floor(value)) {
            stop(paste("invalid fixture dimension:", field), call. = FALSE)
        }
    }
    if (recipe$n_samples < 8 || recipe$n_samples %% 4 != 0 ||
        recipe$n_proteins < 8 || recipe$forms_per_peptide > 2) {
        stop("fixture requires >=8 proteins/samples, samples divisible by 4, <=2 forms",
            call. = FALSE)
    }
    rows <- with(recipe, n_proteins * n_samples * peptides_per_protein * forms_per_peptide)
    cells <- recipe$n_proteins * recipe$n_samples
    if (rows > .Machine$integer.max || cells > .Machine$integer.max) {
        stop("fixture exceeds supported index range", call. = FALSE)
    }
    if (rows > 200000 || cells > 100000 || recipe$annotation_chars > 256 ||
        recipe$history_depth > 16) {
        if (!isTRUE(allow_large)) stop("large fixture requires explicit opt-in", call. = FALSE)
        if (is.null(repo_root)) repo_root <- here::here()
        validateWorkflowFixtureRoot(output_root, repo_root)
    }
    recipe
}

#' Build balanced biological and technical replicate metadata.
#' @param n_samples Sample count divisible by four.
#' @return Ordered design data frame with rich metadata types.
#' @noRd
makeWorkflowCheckpointDesign <- function(n_samples) {
    stopifnot(n_samples >= 8, n_samples %% 4 == 0)
    samples <- sprintf("run_%04d", seq_len(n_samples))
    groups <- rep(c("control", "treated"), each = n_samples / 2)
    biological <- rep(rep(seq_len(n_samples / 4), each = 2), 2)
    data.frame(
        Run = samples,
        group = factor(groups, levels = c("control", "treated")),
        replicates = biological,
        tech_rep_group = paste(groups, biological, sep = "_"),
        batch = factor(rep(c("batch_1", "batch_2"), n_samples / 2),
            levels = c("batch_1", "batch_2", "unused_batch")),
        acquisition = ordered(rep(c("early", "late"), n_samples / 2),
            levels = c("early", "middle", "late")),
        excluded = FALSE,
        row.names = samples, stringsAsFactors = FALSE
    )
}

#' Encode synthetic peptide identities using only amino-acid letters.
#' @param indices Positive identity integers.
#' @return Unique synthetic peptide sequences, not biological annotations.
#' @noRd
workflowCheckpointSequences <- function(indices) {
    alphabet <- strsplit("ACDEFGHIKLMNPQRSTVWY", "", fixed = TRUE)[[1L]]
    value <- as.double(indices) - 1
    stopifnot(all(value >= 0), all(value < 20^6))
    parts <- vector("list", 6L)
    for (position in seq_len(6L)) {
        parts[[position]] <- alphabet[(value %% 20) + 1]
        value <- floor(value / 20)
    }
    paste0("M", do.call(paste0, parts), "K")
}

#' Build one deterministic precursor chunk, without retaining earlier chunks.
#' @param recipe Validated recipe.
#' @param protein_indices Explicit ordered protein subset.
#' @return DIA-NN-compatible long data with independent intensity noise.
#' @noRd
makeWorkflowPrecursorChunk <- function(recipe, protein_indices) {
    stopifnot(length(protein_indices) > 0L, !anyDuplicated(protein_indices),
        all(is.finite(protein_indices)), all(protein_indices == floor(protein_indices)),
        all(protein_indices >= 1L), all(protein_indices <= recipe$n_proteins))
    chunk_rows <- length(protein_indices) * recipe$n_samples *
        recipe$peptides_per_protein * recipe$forms_per_peptide
    if (!is.finite(chunk_rows) || chunk_rows > 200000) {
        stop("precursor chunk exceeds the 200000-row allocation boundary", call. = FALSE)
    }
    # Per-protein seed makes chunks reproducible independently of partitioning.
    chunks <- lapply(protein_indices, \(protein) {
        withr::with_seed(as.integer(recipe$seed + protein), {
            ids <- expand.grid(
                form = seq_len(recipe$forms_per_peptide),
                peptide = seq_len(recipe$peptides_per_protein),
                sample = seq_len(recipe$n_samples)
            )
            sequence <- workflowCheckpointSequences(
                (protein - 1) * recipe$peptides_per_protein + ids$peptide)
            modified <- ifelse(ids$form == 1, sequence,
                paste0("M[Oxidation]", substring(sequence, 2)))
            effect <- if (protein %% 10 == 0) 1 else 0
            log_quantity <- 22 + (protein %% 17) / 5 + ids$peptide / 8 +
                effect * (ids$sample > recipe$n_samples / 2) +
                stats::rnorm(nrow(ids), sd = 0.35)
            quantity <- 2^log_quantity
            # MAR plus intensity-dependent missingness; no imputation here.
            probability <- 0.02 + 0.10 * stats::plogis(22 - log_quantity)
            quantity[stats::runif(length(quantity)) < probability] <- NA_real_
            data.frame(
                Run = sprintf("run_%04d", ids$sample),
                Protein.Group = sprintf("PG_%06d", protein),
                Protein.Ids = sprintf("P%05d", protein),
                Stripped.Sequence = sequence, Modified.Sequence = modified,
                Precursor.Charge = as.integer(ids$form + 1),
                Precursor.Id = paste0(modified, "_", ids$form + 1),
                Precursor.Quantity = quantity,
                Precursor.Normalised = quantity * 1.03,
                Q.Value = stats::runif(nrow(ids), 0.0001, 0.009),
                Global.Q.Value = 0.001, Global.PG.Q.Value = 0.001,
                PG.Q.Value = 0.001, Proteotypic = 1L,
                stringsAsFactors = FALSE
            )
        })
    })
    result <- do.call(rbind, chunks)
    rownames(result) <- NULL
    result
}

#' Construct a small real peptide S4 fixture.
#' @param recipe Fixture recipe; large materialization is deliberately forbidden.
#' @return PeptideQuantitativeData in linear precursor units.
#' @noRd
makeWorkflowCheckpointPeptide <- function(recipe) {
    validateWorkflowCheckpointRecipe(recipe)
    MultiScholaR::PeptideQuantitativeDataDiann(
        peptide_data = makeWorkflowPrecursorChunk(recipe, seq_len(recipe$n_proteins)),
        design_matrix = makeWorkflowCheckpointDesign(recipe$n_samples),
        technical_replicate_id = "tech_rep_group",
        args = list(globalParameters = list(workflow_type = "DIA"),
            fixture = list(recipe = recipe$id, seed = recipe$seed,
                units = "linear precursor intensity", synthetic = TRUE))
    )
}

#' Construct one real protein S4 allocation fixture, not an IQ fit result.
#' @param recipe Fixture recipe.
#' @param generation Nonnegative synthetic quantity-generation index.
#' @param allow_large Explicit large-data opt-in.
#' @param output_root Dedicated external physical-disk directory.
#' @param repo_root Explicit source root for external service working directories.
#' @return Valid ProteinQuantitativeData with synthetic log2 quantities.
#' @noRd
makeWorkflowCheckpointProtein <- function(recipe, generation = 0L,
                                          allow_large = FALSE, output_root = NULL,
                                          repo_root = NULL) {
    validateWorkflowCheckpointRecipe(recipe, allow_large, output_root, repo_root)
    stopifnot(length(generation) == 1L, is.numeric(generation), is.finite(generation),
        generation >= 0, generation == floor(generation), generation <= 10000)
    design <- makeWorkflowCheckpointDesign(recipe$n_samples)
    proteins <- sprintf("PG_%06d", seq_len(recipe$n_proteins))
    values <- withr::with_seed(as.integer(recipe$seed + 100000 + generation), {
        baseline <- stats::rnorm(recipe$n_proteins, 23, 2)
        columns <- lapply(seq_len(recipe$n_samples), \(sample) {
            effect <- (seq_len(recipe$n_proteins) %% 10 == 0) *
                (sample > recipe$n_samples / 2)
            baseline + effect + stats::rnorm(recipe$n_proteins, sd = 0.35)
        })
        stats::setNames(columns, design$Run)
    })
    quant <- as.data.frame(c(list(Protein.Group = proteins), values),
        check.names = FALSE, stringsAsFactors = FALSE)
    annotation <- withr::with_seed(as.integer(recipe$seed + 200000), {
        alphabet <- c(letters, LETTERS, as.character(0:9))
        vapply(seq_len(recipe$n_proteins), \(i)
            paste(sample(alphabet, recipe$annotation_chars, replace = TRUE), collapse = ""),
            character(1))
    })
    identifiers <- data.frame(Protein.Group = proteins,
        Protein.Ids = sprintf("P%05d", seq_len(recipe$n_proteins)),
        gene_names = paste0("SYNTHETIC_", seq_len(recipe$n_proteins)),
        annotation = annotation, stringsAsFactors = FALSE)
    methods::new("ProteinQuantitativeData",
        protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = identifiers, design_matrix = design,
        sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
        args = list(globalParameters = list(workflow_type = "DIA"),
            fixture = list(recipe = recipe$id, generation = generation, synthetic = TRUE,
                units = "log2 protein intensity", estimator_executed = FALSE),
            nested = list(factor = factor(c("b", NA, "a"), levels = c("b", "a", "unused")),
                named = c(first = 1L, second = 2L)))
    )
}

#' Create an explicit codec-edge variant without claiming scientific validity.
#' @param object A small protein fixture.
#' @return S4 object preserving nonfinite values and rich metadata.
#' @noRd
makeWorkflowCheckpointEdges <- function(object) {
    quant <- methods::slot(object, "protein_quant_table")
    stopifnot(nrow(quant) >= 8L, ncol(quant) >= 3L, nrow(quant) * ncol(quant) <= 100000)
    quant[1:6, 2] <- c(NA_real_, NaN, Inf, -Inf, 0, -0)
    quant[7, -1] <- NA_real_
    quant[8, -1] <- 17
    attr(quant, "fixture_units") <- "log2; codec edge cases, not estimator input"
    methods::slot(object, "protein_quant_table") <- quant
    methods::slot(object, "args")$fixture$codec_edges <- TRUE
    object
}

#' Build deterministic uncertainty identity controls using the real EList class.
#' @param object Small protein fixture.
#' @param change Which aligned support/input change to introduce.
#' @return EList-shaped identity fixture; never claimed as real dpcQuant output.
#' @noRd
makeWorkflowUncertaintyControl <- function(object,
    change = c("none", "standard_error", "observations", "order")) {
    change <- match.arg(change)
    quant <- methods::slot(object, "protein_quant_table")
    stopifnot(nrow(quant) * ncol(quant) <= 100000)
    expression <- as.matrix(quant[-1])
    rownames(expression) <- quant[[1]]
    support <- matrix(0.2, nrow(expression), ncol(expression), dimnames = dimnames(expression))
    counts <- matrix(3L, nrow(expression), ncol(expression), dimnames = dimnames(expression))
    if (change == "standard_error") support[1, 1] <- 0.4
    if (change == "observations") counts[1, 1] <- 4L
    result <- methods::new("EList", list(
        E = expression,
        genes = data.frame(protein.id = rownames(expression), row.names = rownames(expression)),
        other = list(standard.error = support, n.observations = counts)
    ))
    if (change == "order") {
        index <- rev(seq_len(nrow(expression)))
        result$E <- result$E[index, , drop = FALSE]
        result$genes <- result$genes[index, , drop = FALSE]
        result$other <- lapply(result$other, \(x) x[index, , drop = FALSE])
    }
    attr(result, "fixture_only_not_dpcQuant") <- TRUE
    result
}

#' Read the actual current cgroup constraints for fixture preflight.
#' @return Limit and peak metadata; errors outside the required boundary.
#' @noRd
workflowFixtureCgroupReceipt <- function() {
    # Linux kernel interfaces have fixed paths; these are not project paths.
    membership <- readLines("/proc/self/cgroup", warn = FALSE)
    unified <- grep("^0::", membership, value = TRUE)
    if (length(unified) != 1L) stop("a cgroup v2 service is required")
    cgroup <- file.path("/sys/fs/cgroup", sub("^0::/", "", unified))
    read_counter <- function(name) {
        value <- readLines(file.path(cgroup, name), n = 1L, warn = FALSE)
        if (identical(value, "max")) return(Inf)
        as.numeric(value)
    }
    memory_max <- read_counter("memory.max")
    swap_max <- read_counter("memory.swap.max")
    if (!is.finite(memory_max) || memory_max > 8 * 1024^3 ||
        !identical(swap_max, 0)) {
        stop("preflight requires memory.max <=8 GiB and memory.swap.max=0")
    }
    list(cgroup = cgroup, memory_max = memory_max, swap_max = swap_max,
        memory_peak = read_counter("memory.peak"),
        swap_current = read_counter("memory.swap.current"))
}


#' Materialize and validate one recipe in a fresh resource-limited process.
#' @param recipe_id Recipe identifier.
#' @param output_root Explicit dedicated disk-backed parent.
#' @param repo_root Repository source root.
#' @return Small preflight receipt; no history or oracle is retained.
#' @noRd
preflightWorkflowCheckpointRecipe <- function(recipe_id, output_root, repo_root) {
    resource_limits <- workflowFixtureCgroupReceipt()
    validateWorkflowFixtureRoot(output_root, repo_root)
    recipe <- readWorkflowCheckpointRecipe(recipe_id, repo_root)
    validateWorkflowCheckpointRecipe(recipe, TRUE, output_root, repo_root)
    case_root <- file.path(output_root, recipe$id)
    if (file.exists(case_root)) stop("preflight output already exists; preserve its evidence")
    if (!dir.create(case_root)) stop("could not create recipe output directory")
    report_path <- file.path(case_root, "report.tsv")
    receipt_path <- file.path(case_root, "preflight.json")
    receipt <- list(recipe_id = recipe$id, seed = recipe$seed,
        kind = "fixture_resource_preflight_not_performance_comparison",
        resource_limits = resource_limits,
        helper_sha256 = digest::digest(file = file.path(repo_root, "tests", "testthat",
            "helper-workflow-checkpoint-fixtures.R"), algo = "sha256"),
        scenarios_sha256 = digest::digest(file = file.path(repo_root, "tests", "testdata",
            "workflow-checkpoints", "scenarios.json"), algo = "sha256"),
        status = "running", report_path = report_path, rows_written = 0,
        max_chunk_rows = 0, runtime_source = "pkgload", started = as.character(Sys.time()))
    on.exit({
        receipt$resources_at_exit <- workflowFixtureCgroupReceipt()
        receipt$finished <- as.character(Sys.time())
        receipt$report_bytes <- if (file.exists(report_path)) file.size(report_path) else 0
        jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE, pretty = TRUE)
    }, add = TRUE)
    # A killed process cannot run on.exit; leave a bound running receipt first.
    jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE, pretty = TRUE)
    tryCatch({
        proteins_per_chunk <- max(1L, floor(200000 /
            (recipe$n_samples * recipe$peptides_per_protein * recipe$forms_per_peptide)))
        for (start in seq.int(1L, recipe$n_proteins, by = proteins_per_chunk)) {
            indices <- seq.int(start, min(recipe$n_proteins, start + proteins_per_chunk - 1L))
            chunk <- makeWorkflowPrecursorChunk(recipe, indices)
            vroom::vroom_write(chunk, report_path, delim = "\t", append = start != 1L,
                col_names = start == 1L, progress = FALSE)
            receipt$rows_written <- receipt$rows_written + nrow(chunk)
            receipt$max_chunk_rows <- max(receipt$max_chunk_rows, nrow(chunk))
            rm(chunk)
        }
        gc()
        # Full input is opened alone, never alongside a full oracle/history.
        raw <- as.data.frame(vroom::vroom(report_path, delim = "\t",
            show_col_types = FALSE, progress = FALSE, altrep = FALSE))
        peptide <- MultiScholaR::PeptideQuantitativeDataDiann(
            peptide_data = raw, design_matrix = makeWorkflowCheckpointDesign(recipe$n_samples),
            technical_replicate_id = "tech_rep_group",
            args = list(fixture = list(recipe = recipe$id, synthetic = TRUE)))
        methods::validObject(peptide)
        stopifnot(nrow(raw) == recipe$precursor_rows)
        receipt$peptide_object_bytes <- as.numeric(utils::object.size(peptide))
        rm(raw, peptide)
        gc()
        protein <- makeWorkflowCheckpointProtein(recipe, 0L, TRUE, output_root, repo_root)
        methods::validObject(protein)
        receipt$protein_object_bytes <- as.numeric(utils::object.size(protein))
        rm(protein)
        gc()
        receipt$status <- "passed_resource_preflight"
        receipt
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        stop(error)
    })
}

#' Construct the finite WFX contrasting-operation qualification population.
#' @param repeats Number of complete four-sample blocks; bounded to 1000.
#' @return Input, design, parameters and hand-authored expected sample decisions.
#' @noRd
makeWorkflowExecutionFilterFixture <- function(repeats = 1L) {
    stopifnot(length(repeats) == 1L, is.numeric(repeats), is.finite(repeats),
        repeats >= 1L, repeats <= 1000L, repeats == floor(repeats))
    # A passes using NA/zero identities; B passes only by inclusion; C is removed.
    # D has four peptides on one protein: protein-only counting must fail.
    # Duplicate identities cross small reader batch boundaries.
    block <- data.frame(Run = c("A", "B", "C", "D", "A", "D", "A", "D", "D", "D"),
        Protein.Group = c("P1", "P1", "P1", "P1", "P2", "P1", "P1", "P1", "P1", "P1"),
        Stripped.Sequence = c("PEPTIDE", "PEPTIDE", "PEPTIDE", "D1", "PEPTIDE",
            "D2", "PEPTIDE", "D3", "D4", "D4"),
        Quantity = c(NA_real_, 1, 0, 1, 0, 2, NA, 3, 4, NA),
        Protein.Ids = c("P1; P2", rep("P1", 9)), stringsAsFactors = FALSE)
    blocks <- lapply(seq_len(repeats), \(i) {
        value <- block
        value$Run <- sprintf("%s_%04d", value$Run, i)
        value
    })
    input <- do.call(rbind, blocks)
    input$source_ordinal <- seq_len(nrow(input))
    samples <- unlist(lapply(seq_len(repeats), \(i) sprintf("%s_%04d", LETTERS[1:4], i)))
    input$Run <- factor(input$Run, levels = c(rev(samples), "unused"))
    design <- data.frame(Run = samples, group = rep(c("case", "control"), 2L * repeats))
    included <- sprintf("B_%04d", seq_len(repeats))
    list(input = input, design = design, threshold = 2L,
        inclusion = c(paste0(" ", included, " "), included, "absent", "", NA_character_),
        expected = list(counts_per_block = c(A = 2L, B = 1L, C = 1L, D = 4L),
            retained = samples[!startsWith(samples, "C_")],
            removed = samples[startsWith(samples, "C_")], absent_inclusion = "absent",
            retained_source_ordinal = input$source_ordinal[!startsWith(as.character(input$Run), "C_")],
            input_rows = 10L * repeats, retained_rows = 9L * repeats,
            output_run_type = "character", unused_level = "unused"))
}
# Physical pages expose automatic row names as logical-position strings, while
# full S4 assembly retains the native compact representation.
workflowPhysicalPageOracle <- function(table, offset = 0L) {
    if (!inherits(table, "tbl_df") && .row_names_info(table, type = 1L) <= 0L)
        rownames(table) <- as.character(offset + seq_len(nrow(table)))
    table
}
