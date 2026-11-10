# Audit A111-06 (2026-11-11). The checkpoint fixture helper reads tests/testdata/workflow-checkpoints/scenarios.json
# through here::here(). Under the R CMD check layout no project root exists for here to find, so every block that
# builds a checkpoint fixture errored there and ran only under load_all. The shipped tests/testdata sits beside
# tests/testthat in both layouts, so here's root is declared from this file: the repository root in a source tree
# (unchanged), the check directory under R CMD check. The helper itself is bound by the journey harness and is not
# edited.
# testthat runs setup files with tests/testthat as the working directory.
if (file.exists(file.path("..", "testdata", "workflow-checkpoints", "scenarios.json")) &&
    file.exists("setup-here-root.R")) {
    suppressMessages(here::i_am("tests/testthat/setup-here-root.R"))
}
