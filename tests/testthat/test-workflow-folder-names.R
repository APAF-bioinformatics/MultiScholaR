# KIT-052 option A (protocol v1.11 item 6): the folders the package writes have short names (store `.we`, staging and
# request root `.wx`, worker request `w-<16>`, session folder `.n-<8>`); a folder written before keeps its name and opens.

test_that("a new store is .we; a root holding .workflow-entities opens that folder and makes no other", {
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "folder-names", "DIA")
    expect_identical(basename(store$directory), ".we")
    expect_identical(store$directory, normalizePath(workflowEntityDirectory(root), winslash = "/"))
    expect_false(dir.exists(file.path(root, ".workflow-entities")))
    identity <- store$identity
    closeWorkflowEntityStore(store)
    # Compatibility: the same store under its earlier name, as a project from before protocol v1.11 holds it.
    expect_true(file.rename(file.path(root, ".we"), file.path(root, ".workflow-entities")))
    reopened <- newWorkflowEntityStore(root, "folder-names", "DIA")
    expect_identical(basename(reopened$directory), ".workflow-entities")
    expect_identical(reopened$identity, identity)
    expect_false(dir.exists(file.path(root, ".we")))
    closeWorkflowEntityStore(reopened)
})

test_that("workflowExportPaths makes .wx in a fresh folder and uses .workflow-export where it exists", {
    id <- "w-0123456789abcdef"
    fresh <- withr::local_tempdir()
    paths <- workflowExportPaths(fresh, id)
    expect_identical(basename(paths$staging_root), ".wx")
    expect_identical(paths$request_dir, file.path(paths$staging_root, id))
    expect_true(dir.exists(paths$request_dir))
    expect_false(dir.exists(file.path(fresh, ".workflow-export")))
    # Compatibility: a folder holding the earlier staging root (a request or publication left pending there).
    earlier <- withr::local_tempdir()
    dir.create(file.path(earlier, ".workflow-export"))
    paths <- workflowExportPaths(earlier, id)
    expect_identical(basename(paths$staging_root), ".workflow-export")
    expect_identical(paths$request_dir, file.path(paths$staging_root, id))
    expect_false(dir.exists(file.path(earlier, ".wx")))
})

test_that("a worker request id is w-<16> or the earlier wfx-<32>, and both give the 12-character unit key", {
    current <- "w-0123456789abcdef"
    earlier <- paste0("wfx-", strrep("0123456789abcdef", 2L))
    expect_true(isWorkflowRequestId(current))
    expect_true(isWorkflowRequestId(earlier))
    refused <- list(paste0("w-", strrep("a", 15L)), paste0("w-", strrep("a", 17L)), paste0("w-", strrep("A", 16L)),
        paste0("wfx-", strrep("a", 16L)), paste0("wfx-", strrep("a", 31L)), paste0("gui-", strrep("a", 32L)), "wfx-test",
        NA_character_, c(current, current), character(), 1, NULL)
    for (id in refused) expect_false(isWorkflowRequestId(id), info = paste(format(id), collapse = " "))
    # The earlier units were `<stem>-w<characters 5 to 16 of the id>`; an earlier request keeps its unit names.
    expect_identical(workflowRequestUnitKey(earlier), substr(earlier, 5L, 16L))
    expect_identical(workflowRequestUnitKey(current), "0123456789ab")
})

test_that("a new session folder is .n-<8> under its root, and the replaced-folder rule removes one and no user folder", {
    root <- withr::local_tempdir()
    folders <- vapply(1:5, function(i) newWorkflowSessionFolder(root), "")
    expect_true(all(grepl("^[.]n-[a-f0-9]{8}$", basename(folders))))
    expect_identical(unique(dirname(folders)), root)
    expect_false(any(file.exists(folders)))
    paths <- list(source_dir = root)
    owned <- folders[[1L]]
    kept <- folders[[2L]]
    for (folder in c(owned, kept)) dir.create(folder)
    expect_true(removeDIAReplacedSessionFolder(owned, kept, paths))
    expect_false(dir.exists(owned))
    expect_false(removeDIAReplacedSessionFolder(kept, kept, paths))
    expect_true(dir.exists(kept))
    # A user's folders, including ones that only resemble the name, are never touched.
    for (name in c("my-results", ".n-0123abc", ".n-0123abcde", ".n-0123ABCD", "n-0123abcd")) {
        user <- file.path(root, name)
        dir.create(user)
        expect_false(removeDIAReplacedSessionFolder(user, kept, paths), info = name)
        expect_true(dir.exists(user), info = name)
    }
    # Nor an application-named folder outside the project's source directory.
    elsewhere <- file.path(withr::local_tempdir(), basename(owned))
    dir.create(elsewhere)
    expect_false(removeDIAReplacedSessionFolder(elsewhere, kept, paths))
    expect_true(dir.exists(elsewhere))
})

test_that("a native archive's store members sit under .we or .workflow-entities, never both or neither", {
    members <- function(prefix) c("manifest.json", "session-references.rds",
        paste0(prefix, c("/store.json", "/dia-state-head.rds", "/entities/dia-state-abc/commit.json",
            "/entities/dia-state-abc/descriptor.rds", "/entities/dia-state-abc/payload.qs2")),
        paste0("files/", strrep("a", 64L)))
    for (prefix in c(".we", ".workflow-entities")) {
        expect_true(checkDIANativeMemberNames(members(prefix)), info = prefix)
        expect_identical(workflowNativeArchiveStorePrefix(members(prefix)), prefix)
    }
    for (prefix in c(".wx", ".wex", "we", ".workflow-export")) {
        expect_error(checkDIANativeMemberNames(members(prefix)), class = "saved_entity_path")
    }
    both <- unique(c(members(".we"), members(".workflow-entities")))
    expect_true(checkDIANativeMemberNames(both))
    expect_error(workflowNativeArchiveStorePrefix(both), class = "saved_entity_corrupt")
    expect_error(workflowNativeArchiveStorePrefix(c("manifest.json", "session-references.rds", paste0("files/", strrep("a", 64L)))),
        class = "saved_entity_corrupt")
})
