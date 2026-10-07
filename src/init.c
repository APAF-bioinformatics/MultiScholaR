#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern const char *workflowRegisterDIAOrderedSum(void *, const char *);
extern SEXP workflow_reclaim_pages(SEXP);
extern SEXP workflow_release_file_pages(SEXP, SEXP, SEXP, SEXP);
extern SEXP workflow_allocator_snapshot(void);
extern SEXP workflow_dia_ordered_new(void);
extern SEXP workflow_dia_ordered_push(SEXP, SEXP, SEXP);
extern SEXP workflow_dia_ordered_close(SEXP);
static const R_CallMethodDef call_methods[] = {
    {"workflow_reclaim_pages", (DL_FUNC) &workflow_reclaim_pages, 1},
    {"workflow_release_file_pages", (DL_FUNC) &workflow_release_file_pages, 4},
    {"workflow_allocator_snapshot", (DL_FUNC) &workflow_allocator_snapshot, 0},
    {"workflow_dia_ordered_new", (DL_FUNC) &workflow_dia_ordered_new, 0},
    {"workflow_dia_ordered_push", (DL_FUNC) &workflow_dia_ordered_push, 3},
    {"workflow_dia_ordered_close", (DL_FUNC) &workflow_dia_ordered_close, 1},
    {NULL, NULL, 0}
};

void attribute_visible R_init_MultiScholaR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_RegisterCCallable("MultiScholaR", "dia_ordered_sum_v1", (DL_FUNC) &workflowRegisterDIAOrderedSum);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
