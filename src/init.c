#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

extern SEXP C_multischolar_malloc_trim(void);

static const R_CallMethodDef call_methods[] = {
    {"C_multischolar_malloc_trim", (DL_FUNC) &C_multischolar_malloc_trim, 0},
    {NULL, NULL, 0}
};

void R_init_MultiScholaR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
