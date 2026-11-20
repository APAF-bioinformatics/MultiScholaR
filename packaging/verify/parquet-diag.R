# run-win-k147 diagnostic (instrumentation, not package code): do the corpus Parquet files composed on this platform
# differ from the committed (Linux-made) ones in content or only in encoding? Run with the bundle's R.
# Usage: Rscript --vanilla parquet-diag.R <checkout> <out dir>
# Prints DIAG lines to the log, no length cut.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) == 3L) arguments <- arguments[-1L]
checkout <- normalizePath(arguments[[1L]], mustWork = TRUE)
out <- arguments[[2L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
say <- function(...) cat("DIAG", sprintf(...), "\n")
loadNamespace("MultiScholaR")

say("R %s; arrow %s; MultiScholaR %s", getRversion(), utils::packageVersion("arrow"), utils::packageVersion("MultiScholaR"))
say("Sys.getlocale(): %s", Sys.getlocale())
l10n <- l10n_info()
say("l10n_info(): %s", paste(names(l10n), unlist(l10n), sep = "=", collapse = "; "))
say("LANG=%s LC_ALL=%s LC_COLLATE=%s", Sys.getenv("LANG"), Sys.getenv("LC_ALL"), Sys.getenv("LC_COLLATE"))
info <- arrow::arrow_info()
say("arrow capabilities: %s", paste(names(info$capabilities), info$capabilities, sep = "=", collapse = " "))

# 1(a): the checked-out manifest-listed json
json <- file.path(checkout, "tests", "testdata", "integration", "corpus", "proteomics", "integration_input.json")
bytes <- readBin(json, "raw", file.size(json))
say("1a %s: %d bytes, %d CR bytes, %d LF bytes, sha256 %s", basename(json), length(bytes), sum(bytes == as.raw(13L)),
    sum(bytes == as.raw(10L)), digest::digest(file = json, algo = "sha256"))
lf <- bytes[bytes != as.raw(13L)]
say("1a the same file with CR bytes removed: sha256 %s", digest::digest(lf, algo = "sha256", serialize = FALSE))
committed_manifest <- jsonlite::read_json(file.path(checkout, "tests", "testdata", "integration", "corpus", "manifest.json"))
say("1a committed sha256 of proteomics/integration_input.json: %s", committed_manifest$files[["proteomics/integration_input.json"]]$sha256)

# compose
h <- new.env(parent = globalenv())
sys.source(file.path(checkout, "tests", "testdata", "integration", "compose_integration_corpus.R"), h)
fresh_dir <- file.path(out, "fresh")
h$composeIntegrationCorpus(fresh_dir, repo = checkout)
committed_dir <- file.path(checkout, "tests", "testdata", "integration", "corpus")
say("composed into %s", fresh_dir)

write_parquet <- get("writeIntegrationInputParquet", envir = asNamespace("MultiScholaR"))
for (name in c("samples.parquet", "features_proteomics.parquet", "matrix_proteomics.parquet")) {
    a <- file.path(committed_dir, "proteomics", name)
    b <- file.path(fresh_dir, "proteomics", name)
    say("== %s: committed %d bytes sha256 %s; fresh %d bytes sha256 %s", name, file.size(a),
        digest::digest(file = a, algo = "sha256"), file.size(b), digest::digest(file = b, algo = "sha256"))
    ta <- as.data.frame(arrow::read_parquet(a))
    tb <- as.data.frame(arrow::read_parquet(b))
    say("%s dims committed %s fresh %s; names identical %s", name, paste(dim(ta), collapse = "x"),
        paste(dim(tb), collapse = "x"), identical(names(ta), names(tb)))
    if (!identical(names(ta), names(tb))) {
        say("%s names committed: %s", name, paste(names(ta), collapse = ","))
        say("%s names fresh: %s", name, paste(names(tb), collapse = ","))
    }
    say("%s whole table identical(): %s", name, identical(ta, tb))
    for (column in intersect(names(ta), names(tb))) {
        x <- ta[[column]]
        y <- tb[[column]]
        same <- identical(x, y)
        if (is.character(x) && is.character(y)) {
            ex <- paste(names(table(Encoding(x))), table(Encoding(x)), collapse = " ")
            ey <- paste(names(table(Encoding(y))), table(Encoding(y)), collapse = " ")
            differ <- which(x != y | is.na(x) != is.na(y))
            same_set <- setequal(x, y)
            say("%s column %s (text): identical %s; Encoding committed [%s] fresh [%s]; same values as a set %s; rows differing %d%s",
                name, column, same, ex, ey, same_set, length(differ),
                if (length(differ)) sprintf("; first: %s", paste(sprintf("row %d committed '%s' fresh '%s'", head(differ, 5),
                    x[head(differ, 5)], y[head(differ, 5)]), collapse = " | ")) else "")
            if (same_set && length(differ)) say("%s column %s: same values, different row order", name, column)
        } else if (is.numeric(x) && is.numeric(y)) {
            d <- abs(x - y)
            say("%s column %s (%s/%s): identical %s; NA patterns identical %s; differing %d; largest abs difference %s",
                name, column, typeof(x), typeof(y), same, identical(is.na(x), is.na(y)), sum(d > 0, na.rm = TRUE),
                format(suppressWarnings(max(d, na.rm = TRUE)), digits = 17))
        } else {
            say("%s column %s: types %s / %s; identical %s", name, column, class(x)[[1L]], class(y)[[1L]], same)
        }
    }
    # Decisive: the committed content re-written here with the package's writer.
    again <- file.path(out, paste0("rewritten-", name))
    write_parquet(ta, again)
    say("%s committed content re-written on this platform: sha256 %s (committed %s; equal %s)", name,
        digest::digest(file = again, algo = "sha256"), digest::digest(file = a, algo = "sha256"),
        identical(digest::digest(file = again, algo = "sha256"), digest::digest(file = a, algo = "sha256")))
    tb_again <- file.path(out, paste0("rewritten-fresh-", name))
    write_parquet(tb, tb_again)
    say("%s fresh content re-written on this platform: sha256 %s (fresh %s; equal %s)", name,
        digest::digest(file = tb_again, algo = "sha256"), digest::digest(file = b, algo = "sha256"),
        identical(digest::digest(file = tb_again, algo = "sha256"), digest::digest(file = b, algo = "sha256")))
}

# sort evidence: the composer's sort (compose_integration_corpus.R:361), the generator's order (kit_omics_corpus.R:55)
recipe <- h$integrationCorpusRecipe(h$INTEGRATION_CORPUS_SEED)
batches <- unique(unlist(lapply(recipe$omics, `[[`, "batches")))
probe <- c("B1", "b1", "B10", "B2", "batch_B1", "batch_b2", "KO_10", "KO_2", "KO_1", "WT_1", "ko_1", "_x", "x_1", "x1", "X1",
    "trx_1", "trx_10", "TRX_2", "prot-1", "prot_1", "prot.1", "a b", "ab", "A")
show <- function(x) paste(x, collapse = ",")
say("sort batches default: %s; radix: %s", show(sort(batches)), show(sort(batches, method = "radix")))
say("sort probe default: %s", show(sort(probe)))
say("sort probe radix:   %s", show(sort(probe, method = "radix")))
old <- Sys.getlocale("LC_COLLATE")
invisible(Sys.setlocale("LC_COLLATE", "C"))
say("sort probe C locale: %s", show(sort(probe)))
invisible(Sys.setlocale("LC_COLLATE", old))
say("default equals radix for probe: %s; for batches: %s", identical(sort(probe), sort(probe, method = "radix")),
    identical(sort(batches), sort(batches, method = "radix")))
say("done")
