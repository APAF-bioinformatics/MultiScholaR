# KIT-167: fixtures shared by the reference mapping test files. The engine is replaced by a stub that writes a
# recorded-form DIAMOND table (one hit per protein, strain entity i against reference entry i), so the whole mapping
# path runs without DIAMOND; one lane runs the real engine (section 8.6).

#' A strain FASTA of `n` RefSeq-style entries holding the first `n` sequences of a shipped reference
k167StrainFromReference <- function(upid = "UP000000750", n = 60L) {
    proteome <- readReferenceProteome(referenceOrganismFiles(upid)$fasta)[seq_len(n), ]
    path <- tempfile(fileext = ".fasta")
    writeLines(paste0(">ref|WP_", sprintf("%09d", seq_len(n)), ".1| MULTISPECIES: protein [Streptococcus]\n",
        proteome$sequence), path)
    list(path = path, reference_accessions = proteome$accession)
}

#' A stub engine: `version` prints 2.2.8, `makedb` does nothing, `blastp` writes `n` hits pairing entity i with
#' reference entry i in both directions (identity 99.5, both coverages 100)
k167StubDiamond <- function(n = 60L) {
    testthat::skip_on_os("windows")
    path <- file.path(tempfile("k167-stub-"), "diamond")
    dir.create(dirname(path))
    writeLines(c("#!/bin/sh",
        "case \"$1\" in",
        "  version) echo 'diamond version 2.2.8' ;;",
        "  makedb) exit 0 ;;",
        "  blastp)",
        "    out=''; query=''",
        "    while [ $# -gt 0 ]; do case \"$1\" in --out) out=\"$2\"; shift ;; --query) query=\"$2\"; shift ;; esac; shift; done",
        "    if [ \"$query\" = strain.fasta ]; then a=s; b=r; else a=r; b=s; fi",
        "    : > \"$out\"; i=1",
        sprintf("    while [ $i -le %d ]; do", as.integer(n)),
        "      printf '%s%07d\\t%s%07d\\t99.5\\t100\\t0\\t0\\t1\\t100\\t1\\t100\\t1e-50\\t200\\t500\\t100\\t100\\t100\\t100\\n' $a $i $b $i >> \"$out\"",
        "      i=$((i+1))",
        "    done ;;",
        "esac"), path)
    Sys.chmod(path, "0755")
    path
}

#' Replace the engine check for the caller's lifetime with one returning the stub
k167UseStubEngine <- function(stub, env = parent.frame()) {
    testthat::local_mocked_bindings(checkDiamond = function(...) list(path = stub, version = "2.2.8",
        sha256 = digest::digest(file = stub, algo = "sha256"), platform = "linux"), .env = env)
}

#' A session owner on a temporary store, with reference states on
k167Owner <- function(env = parent.frame()) {
    root <- tempfile("k167-owner-")
    dir.create(root)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "k167", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    withr::defer(closeDIAStateAdapter(manager), envir = env)
    manager
}

#' Network calls made to fail (a curated mapping and a curated-reference annotation make none)
k167NoNetwork <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(req_perform = function(...) stop("network call made"), .package = "httr2",
        .env = env)
    testthat::local_mocked_bindings(curl_fetch_memory = function(...) stop("network call made"),
        curl_fetch_disk = function(...) stop("network call made"), .package = "curl", .env = env)
}
