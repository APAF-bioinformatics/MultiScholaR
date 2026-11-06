# KIT-099: KIT-092's reverify.R (protocol v1.7) with one change: the decision is given the member directory, the argument
# KIT-098 added to kitSingleJourneyDecision for the transcriptomics worker journals (protocol v1.8 section 10.5). Without it
# a trx member reads as worker_dispatch_missing (kit-099/reverify-attempt1.txt). Re-verifies one recorded single journey
# member offline under the current harness and decides it again on its recorded receipt. No journey runs.
# Usage, from the repository root, with the member's build library first on R_LIBS:
#   Rscript docs/qa/workflow-execution/kit/kit-099/reverify.R <stage> <copy root>
# <copy root>/runtime/<stage> is a reflinked copy of the member's stage directory (the runtime directory is only read,
# by the copy). The verifier writes its receipt into the copy. The copy's journey record points at the copy's store.
# Binding: the member's recorded binding must equal its case under the bindings it was launched with (the stage's
# source snapshot), and under a scratch copy of the current bindings (the copy root; the launch protocol hash and the
# launch KIT-005 canonical-form record, the two binding inputs v1.7 itself changed). So the member is judged under the
# current harness and oracle fingerprints with everything else as launched.
args <- commandArgs(TRUE)
stage <- args[[1L]]
copy_root <- normalizePath(args[[2L]])
repo <- normalizePath(".")
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
h <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
kit <- "docs/qa/workflow-execution/kit/bindings.json"
current <- h$kitBindings(file.path(repo, kit))
runtime <- current$root
stage_copy <- file.path(copy_root, "runtime", stage)
stopifnot(dir.exists(stage_copy), !identical(normalizePath(stage_copy), normalizePath(file.path(runtime, stage))))
launch <- h$kitBindings(file.path(stage_copy, "source", kit))
selection <- paste0("kit:journey-", sub("^kit-journey-", "", stage), ":candidate:1")
target <- h$kitSelectorStage(selection)
operation <- readRDS(file.path(stage_copy, target$id, "operation-receipt.rds"))
launch_case <- Filter(function(x) x$id == target$id, h$kitCases(launch, target$stage))[[1L]]
stopifnot(identical(operation$binding, launch_case$binding))
# The scratch bindings: the current bindings with the copy root and, from the launch bindings, the two binding inputs
# v1.7 itself changed: the protocol hash, and the KIT-005 baseline's canonical-form record (the generator recomputes
# those oracles from control-1's store at every bind, and the case digest carries `journey`).
raw <- jsonlite::read_json(file.path(repo, kit), simplifyVector = FALSE)
launch_raw <- jsonlite::read_json(file.path(stage_copy, "source", kit), simplifyVector = FALSE)
raw$harness$root <- file.path(copy_root, "runtime")
raw$harness$protocol_sha256 <- launch_raw$harness$protocol_sha256
raw$journey$baseline$oracle_canonical_form <- launch_raw$journey$baseline$oracle_canonical_form
scratch_repo <- file.path(copy_root, paste0("repo-", stage))
dir.create(dirname(file.path(scratch_repo, kit)), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(raw, file.path(scratch_repo, kit), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
scratch <- h$kitBindings(file.path(scratch_repo, kit))
case <- Filter(function(x) x$id == target$id, h$kitCases(scratch, target$stage))[[1L]]
stopifnot(identical(operation$binding, case$binding))
# The copy's journey record points at the copy's store and project.
journey_path <- file.path(stage_copy, target$id, "journey-result.rds")
journey <- readRDS(journey_path)
moved <- function(path) if (is.character(path)) sub(paste0(runtime, "/"), paste0(copy_root, "/runtime/"), path, fixed = TRUE) else path
journey$store <- moved(journey$store)
journey$project <- moved(journey$project)
for (step in names(journey$populations)) journey$populations[[step]]$store <- moved(journey$populations[[step]]$store)
stopifnot(startsWith(journey$store, copy_root), startsWith(journey$project, copy_root))
saveRDS(journey, journey_path)
unlink(file.path(stage_copy, target$id, c("verification.rds", "verify")), recursive = TRUE)
run_dir <- file.path(stage_copy, target$id, "verify"); dir.create(run_dir)
started <- Sys.time()
error <- tryCatch({ h$verifyKitMember(run_dir, scratch_repo, selection, api); NULL }, error = function(e) conditionMessage(e))
verification <- readRDS(file.path(stage_copy, target$id, "verification.rds"))
receipt <- readRDS(file.path(stage_copy, "stage.rds"))
recorded <- readRDS(file.path(runtime, stage, "decision.rds"))
receipt$cases[[target$id]]$verifier_completion <- list(accepted = is.null(error) && identical(verification$status, "complete"),
    reproduced = verification$reproduced, reason = if (!is.null(error)) paste("Error:", error))
decision <- h$kitSingleJourneyDecision(receipt, list(case), scratch, member_dir = file.path(stage_copy, target$id))
saveRDS(list(stage = stage, verification = verification, decision = decision, recorded_decision = recorded,
    launch_protocol_sha256 = launch$protocol_sha256, current_protocol_sha256 = current$protocol_sha256),
    file.path(copy_root, paste0(stage, "-reverify.rds")))
products <- function(populations) sum(vapply(populations, function(p) length(p$fingerprint$products), integer(1)))
cat(sprintf("REVERIFY %s | recorded %s%s exact %s | now %s%s exact %s in_memory %s | oracle mismatches %d | additional products %d | operation identity products %d | canonical form differences %s | judged steps %s | %.0f s\n",
    stage, recorded$status, if (length(recorded$failing_rule)) paste0(" (", paste(recorded$failing_rule, collapse = ","), ")") else "",
    format(recorded$member$exact), decision$status,
    if (length(decision$failing_rule)) paste0(" (", paste(decision$failing_rule, collapse = ","), ")") else "",
    format(decision$member$exact), format(decision$member$in_memory), length(verification$oracle_mismatches),
    length(verification$additional_products), length(verification$operation_identity_products),
    paste(verification$canonical_form_differences, collapse = ",") %||% "", paste(names(verification$populations), collapse = ","),
    as.numeric(difftime(Sys.time(), started, units = "secs"))))
if (!is.null(error)) cat("  verifier error:", error, "\n")
if (length(verification$oracle_mismatches)) cat("  mismatches:", paste(verification$oracle_mismatches, collapse = " "), "\n")
