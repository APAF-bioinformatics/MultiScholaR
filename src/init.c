#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP workflow_reclaim_pages(SEXP);
extern SEXP workflow_release_file_pages(SEXP, SEXP, SEXP, SEXP);
extern SEXP workflow_allocator_snapshot(void);
static const R_CallMethodDef call_methods[] = {
    {"workflow_reclaim_pages", (DL_FUNC) &workflow_reclaim_pages, 1},
    {"workflow_release_file_pages", (DL_FUNC) &workflow_release_file_pages, 4},
    {"workflow_allocator_snapshot", (DL_FUNC) &workflow_allocator_snapshot, 0},
    {NULL, NULL, 0}
};

void attribute_visible R_init_MultiScholaR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
