#include <R.h>
#include <Rinternals.h>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

SEXP C_multischolar_malloc_trim(void)
{
#if defined(__GLIBC__)
    return Rf_ScalarLogical(malloc_trim(0) != 0);
#else
    return Rf_ScalarLogical(0);
#endif
}
