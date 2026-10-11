# KIT-011 (audit PX-02): the screen lifecycle member on the kit_011 build (production default, on-disk
# state from import) compared state by state with the KIT-001/KIT-005 screen member produced on the
# kit_001 build (`kit-lifecycle-complete.kit-005-screen`), with the comparator KIT-005's control-2 used.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-011/screen-cross-build-check.R
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
reference <- readRDS(file.path(b$root, "kit-lifecycle-complete.kit-005-screen", "member", "journey-result.rds"))
candidate <- readRDS(file.path(b$root, "kit-lifecycle-complete", "member", "journey-result.rds"))
stage <- readRDS(file.path(b$root, "kit-lifecycle-complete", "stage.rds"))
stopifnot(identical(reference$status, "complete"), identical(candidate$status, "complete"))
strict <- h$kitReproduceMismatches(reference$populations, candidate$populations)
mismatches <- h$kitReproduceMismatches(reference$populations, candidate$populations, ignore_products = h$kitOperationIdentityProducts)
# The reference screen member ran on the kit_001 source before the RUV sweep rules of 2026-10-10 (install
# browser-smoke.sMOGTK, commit 2ff38edf), so its `normalization_science` sweep record differs from every later
# build's while the RUV-corrected state itself is identical; the KIT-005 full-size oracles carry the rules.
explained <- "normalisation_ruv/products/normalization_science"
op <- stage$cases$member$operation
record <- list(schema = "kit.screen_cross_build_check.v1", ticket = "KIT-011", date = format(Sys.Date()),
    reference = list(stage = "kit-lifecycle-complete.kit-005-screen", build = "kit_001 (browser-smoke.sMOGTK era)",
        populations = sum(vapply(reference$populations, function(p) !is.null(p$fingerprint), logical(1)))),
    candidate = list(stage = "kit-lifecycle-complete", build = b$lifecycle$build, install = b$builds[[b$lifecycle$build]]$install_member,
        decision = stage$decision$status, elapsed_seconds = round(op$elapsed_seconds), peak_gib = round(op$peak_bytes / 1024^3, 3),
        anon_peak_gib = round(op$sampled_anon_bytes / 1024^3, 3), resident_end_gib = round(op$resident_end_bytes / 1024^3, 3),
        populations = sum(vapply(candidate$populations, function(p) !is.null(p$fingerprint), logical(1))),
        phase_peak_anon_gib = lapply(op$phase_peak_anon_bytes, function(x) round(x / 1024^3, 3))),
    strict_mismatches = strict, mismatches_excluding_operation_identity_products = mismatches,
    explained = explained, unexplained = setdiff(mismatches, explained), identical_states = !length(setdiff(mismatches, explained)))
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-011/screen-cross-build-check.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
cat("decision", stage$decision$status, "| populations", record$candidate$populations, "| strict mismatches", length(strict),
    "| after identity products", length(mismatches), "| unexplained", length(record$unexplained), "\n")
if (length(strict)) print(strict)
