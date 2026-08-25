args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop("usage: <contract.json> <memory|artifact> <output.json>", call. = FALSE)
}

contract_path <- normalizePath(args[[1L]], mustWork = TRUE)
backend <- match.arg(args[[2L]], c("memory", "artifact"))
output_path <- args[[3L]]
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1L]])
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

devtools::load_all(repo_root, quiet = TRUE)
environment <- new.env(parent = .GlobalEnv)
sys.source(
    file.path(repo_root, "tools", "profiling", "omics_workload_contract.R"),
    envir = environment
)
adapter_path <- file.path(
    repo_root,
    "tools",
    "profiling",
    "omics_workload_lipidomics.R"
)
contract <- environment$omicsWorkloadReadContract(contract_path)
adapter <- environment$omicsWorkloadLoadAdapter(adapter_path, contract)
do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
set.seed(as.integer(contract$rng$seed))
run_dir <- tempfile("lipid-eviction-worker-")
dir.create(run_dir, recursive = TRUE)
on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)
prepared <- adapter$prepare(list(contract = contract, run_dir = run_dir))
payload_digest <- digest::digest(
    file = prepared$payload_path,
    algo = "sha256",
    serialize = FALSE
)
truth_digest <- digest::digest(
    file = prepared$truth_path,
    algo = "sha256",
    serialize = FALSE
)
if (!identical(payload_digest, contract$expected_digests$payload_sha256) ||
    !identical(truth_digest, contract$expected_digests$truth_sha256)) {
    stop("frozen lipidomics workload digest mismatch", call. = FALSE)
}
data <- utils::read.delim(
    prepared$payload_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = "NA"
)
assay_names <- names(contract$assay_mix)
assays <- lapply(assay_names, \(assay_name) {
    assay <- data[data$assay == assay_name, , drop = FALSE]
    assay$assay <- NULL
    row.names(assay) <- NULL
    assay
})
names(assays) <- assay_names
data_tbl <- assays
data_cln <- unserialize(serialize(assays, NULL, version = 3L))
source_bytes <- as.numeric(object.size(data_tbl)) +
    as.numeric(object.size(data_cln))
artifact_bytes <- 0
if (backend == "artifact") {
    artifact_dir <- file.path(run_dir, "artifacts")
    dir.create(artifact_dir)
    for (index in seq_along(assays)) {
        encoded <- encodeArtifactTable(
            assays[[index]],
            owner = paste("lipidomics eviction", assay_names[[index]])
        )
        path <- file.path(artifact_dir, sprintf("assay-%04d.parquet", index))
        do.call(
            arrow::write_parquet,
            c(
                list(x = encoded$payload, sink = path),
                artifactParquetWriteArgs(encoded)
            )
        )
        artifact_bytes <- artifact_bytes + as.numeric(file.info(path)$size)
    }
    rss_before <- as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]])
    data_tbl <- NULL
    data_cln <- NULL
    assays <- NULL
    data <- NULL
    invisible(gc(full = TRUE))
    invisible(gc(full = TRUE))
} else {
    rss_before <- as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]])
}
Sys.sleep(0.25)
rss_after <- as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]])
result <- list(
    schema = "multischolar.lipidomics_eviction_measurement",
    schema_version = "1.0.0",
    workload_id = contract$workload_id,
    contract_digest = environment$omicsWorkloadDigest(contract),
    privacy_classification = contract$privacy$classification,
    backend = backend,
    payload_sha256 = payload_digest,
    truth_sha256 = truth_digest,
    feature_count = as.integer(contract$dimensions$feature_count),
    sample_count = as.integer(contract$dimensions$sample_count),
    assay_count = as.integer(contract$dimensions$assay_count),
    source_bytes = source_bytes,
    artifact_bytes = artifact_bytes,
    rss_before_eviction = rss_before,
    rss_after_eviction = rss_after,
    retained_source_bytes = if (backend == "artifact") 0 else source_bytes,
    generated_claim_boundary = paste(
        "resource and eviction evidence only;",
        "not parser or biological validation"
    )
)
jsonlite::write_json(
    result,
    output_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    digits = 17
)
