# KIT-112 D8: the large pair, store against store, and the heads' DA manifests (KIT-111's compare-stores.R plus the manifest
# comparison, D6). KIT-111's header follows.
# KIT-111 D8: the large pair, store against store. Each member's store is fingerprinted with the harness's canonical
# form (kitFingerprintStore, tools/profiling/workflow_kit.R), and every state's class, rows, components, support,
# summary_sha256 and products, and every head product's fingerprint, must be identical (kitStateFingerprintMatches'
# fields). Product ids are store-local, so head products are compared as the sorted multiset of their fingerprints.
# Rendered pages (kind `rendered_asset`) carry per-run content (times) and are compared by count only, as the harness's
# in-process recorders leave them out (kitStepProducts(rendered = FALSE)); checked on 107l against 108l, two exact
# in-process members, where the one differing head product is the rendered report page.
# Usage (from the repository root): Rscript compare-stores.R <library> <baseline member dir> <candidate member dir> <out tsv>
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
.libPaths(c(args[[1L]], .libPaths()))
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
env <- new.env(parent = globalenv())
sys.source("tools/profiling/workflow_kit.R", env)
fingerprint <- function(member) env$kitFingerprintStore(readRDS(file.path(member, "journey-result.rds"))$store, api)
a <- fingerprint(args[[2L]])
b <- fingerprint(args[[3L]])
digest_of <- function(x) digest::digest(x, algo = "sha256")
rows <- list()
add <- function(item, field, same) rows[[length(rows) + 1L]] <<- data.frame(item = item, field = field, identical = same)
for (state in union(names(a$states), names(b$states))) {
    x <- a$states[[state]]
    y <- b$states[[state]]
    if (is.null(x) || is.null(y)) { add(state, "present", FALSE); next }
    for (field in c("class", "rows", "components", "support", "summary_sha256"))
        add(state, field, identical(x[[field]], y[[field]]))
    add(state, "products", identical(sort(vapply(x$products, digest_of, "")), sort(vapply(y$products, digest_of, ""))))
}
rendered <- function(products) vapply(products, function(p) identical(p$kind, "rendered_asset"), logical(1))
add("head", "product fingerprints (not rendered)", identical(sort(unname(vapply(a$products[!rendered(a$products)], digest_of, ""))),
    sort(unname(vapply(b$products[!rendered(b$products)], digest_of, "")))))
add("head", "rendered page count", identical(sum(rendered(a$products)), sum(rendered(b$products))))
# KIT-112: the heads' `da_results` manifests, by contrast table, parameters, model and each contrast's value, in the
# harness's canonical form (`kitCanonical`: saved references replaced by their content, tables by their pages).
manifest <- function(member) {
    dir <- readRDS(file.path(member, "journey-result.rds"))$store
    store <- env$kitOpenStore(dir, api)
    on.exit(api("closeWorkflowEntityStore")(store), add = TRUE)
    snapshot <- api("hydrateSavedWorkflowCheckpoint")(store, readRDS(file.path(dir, "dia-state-head.rds")))
    if (is.null(snapshot$da_results)) return(NULL)
    value <- api("hydrateSavedWorkflowEntity")(store, snapshot$da_results)
    memo <- new.env(parent = emptyenv())
    list(contrasts = digest_of(value$contrasts), parameters = digest_of(value$parameters), model = digest_of(value$model),
        source = digest_of(env$kitCanonical(value$source, store, api, memo = memo)),
        results = vapply(value$results, function(r) digest_of(env$kitCanonical(r, store, api, memo = memo)), ""))
}
ma <- manifest(args[[2L]])
mb <- manifest(args[[3L]])
add("da_results", "present", identical(is.null(ma), is.null(mb)))
if (!is.null(ma) && !is.null(mb)) {
    for (field in c("contrasts", "parameters", "model", "source")) add("da_results", field, identical(ma[[field]], mb[[field]]))
    add("da_results", "contrast names", identical(names(ma$results), names(mb$results)))
    for (name in union(names(ma$results), names(mb$results))) add("da_results", paste("result", name), identical(ma$results[[name]], mb$results[[name]]))
}
table <- do.call(rbind, rows)
utils::write.table(table, args[[4L]], sep = "\t", quote = FALSE, row.names = FALSE)
cat("STATES", length(union(names(a$states), names(b$states))), "COMPARED", nrow(table), "DIFFERENT", sum(!table$identical), "\n")
if (any(!table$identical)) print(table[!table$identical, ], row.names = FALSE)
