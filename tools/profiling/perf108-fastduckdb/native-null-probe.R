# One isolated public-C provider probe; no installed namespace is patched.
root <- commandArgs(TRUE)[[1L]]
library <- commandArgs(TRUE)[[2L]]
stopifnot(!dir.exists(root), dir.create(root, recursive = TRUE))
.libPaths(c(library, .libPaths()))
loadNamespace("MultiScholaR")
for (name in c("workflow_dia_exact_sum.c", "duckdb.h"))
    stopifnot(file.copy(file.path("src", name), file.path(root, name)))
writeLines(c('#include <R.h>', '#include <Rinternals.h>', '#include <R_ext/Rdynload.h>',
    'extern const char *workflowRegisterDIAOrderedSum(void *, const char *);',
    'void R_init_nullprobe(DllInfo *dll) {',
    ' (void)dll; R_RegisterCCallable("MultiScholaR", "dia_ordered_sum_repair_probe_v1", (DL_FUNC)&workflowRegisterDIAOrderedSum);',
    '}'), file.path(root, "nullprobe.c"))
writeLines("CFLAGS = -std=c11 -O2 -fPIC -fno-fast-math -Wall -Wextra -Werror -Wno-cast-function-type",
    file.path(root, "Makevars"))
status <- withr::with_dir(root, withr::with_envvar(c(
    R_MAKEVARS_USER = file.path(root, "Makevars"), PKG_LIBS = "-ldl"),
    system2(file.path(R.home("bin"), "R"), c("CMD", "SHLIB", "-o", "nullprobe.so",
        "nullprobe.c", "workflow_dia_exact_sum.c"))))
stopifnot(status == 0L)
dyn.load(file.path(root, "nullprobe.so"), local = TRUE)
checks <- withr::with_envvar(c(MULTISCHOLAR_PERF108_NATIVE_NULL = "true",
    MULTISCHOLAR_PERF108_SUM_SYMBOL = "dia_ordered_sum_repair_probe_v1"),
    testthat::test_file("tests/testthat/test-perf108-native-null.R", reporter = "summary"))
checks <- as.data.frame(checks)
saveRDS(checks, file.path(root, "checks.rds"))
stopifnot(!any(checks$error | checks$failed > 0 | checks$skipped))
