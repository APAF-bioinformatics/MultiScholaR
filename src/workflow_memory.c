#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#if defined(__linux__)
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#endif
#if defined(__linux__) && defined(__GLIBC__)
#include <malloc.h>
#endif

/* Capability probe or best-effort release of free pages; never live R objects. */
SEXP workflow_reclaim_pages(SEXP release) {
    if (TYPEOF(release) != LGLSXP || XLENGTH(release) != 1 ||
        LOGICAL(release)[0] == NA_LOGICAL)
        Rf_error("release must be one non-missing logical value");
#if defined(__linux__) && defined(__GLIBC__)
    return Rf_ScalarLogical(LOGICAL(release)[0] ? malloc_trim(0) != 0 : TRUE);
#else
    return Rf_ScalarLogical(NA_LOGICAL);
#endif
}

/* Advice is file-specific and never alters bytes; callers own the file/range. */
SEXP workflow_release_file_pages(SEXP path, SEXP offset, SEXP length, SEXP sync) {
    if (TYPEOF(path) != STRSXP || XLENGTH(path) != 1 || STRING_ELT(path, 0) == NA_STRING ||
        TYPEOF(sync) != LGLSXP || XLENGTH(sync) != 1 || LOGICAL(sync)[0] == NA_LOGICAL ||
        !Rf_isNumeric(offset) || XLENGTH(offset) != 1 || !Rf_isNumeric(length) || XLENGTH(length) != 1)
        Rf_error("invalid file-page release arguments");
    double start = Rf_asReal(offset), bytes = Rf_asReal(length);
    if (!R_FINITE(start) || !R_FINITE(bytes) || start < 0 || bytes < 0 ||
        start != floor(start) || bytes != floor(bytes) || start + bytes > 9007199254740991.0)
        Rf_error("invalid file-page release range");
#if defined(__linux__)
    int fd = open(Rf_translateCharUTF8(STRING_ELT(path, 0)), O_RDONLY | O_CLOEXEC | O_NOFOLLOW);
    if (fd < 0) Rf_error("cannot open owned cache-release file: %s", strerror(errno));
    struct stat info;
    if (fstat(fd, &info) || !S_ISREG(info.st_mode) || start > (double)info.st_size ||
        (bytes > 0 && start + bytes > (double)info.st_size)) {
        close(fd);
        Rf_error("cache-release file or range differs");
    }
    if (LOGICAL(sync)[0] && fdatasync(fd)) {
        int code = errno;
        close(fd);
        Rf_error("cannot synchronize completed output prefix: %s", strerror(code));
    }
    int code = posix_fadvise(fd, (off_t)start, (off_t)bytes, POSIX_FADV_DONTNEED);
    int closed = close(fd);
    if (code) Rf_error("file cache advice failed: %s", strerror(code));
    if (closed) Rf_error("cannot close cache-release descriptor");
    return Rf_ScalarLogical(TRUE);
#else
    return Rf_ScalarLogical(NA_LOGICAL);
#endif
}

/* Bounded allocator counters. These overlap R allocations and are not a heap census. */
SEXP workflow_allocator_snapshot(void) {
#if defined(__linux__) && defined(__GLIBC__)
    struct mallinfo2 info = mallinfo2();
    SEXP values = PROTECT(Rf_allocVector(REALSXP, 5));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
    const char *labels[] = {"arena_bytes", "in_use_arena_bytes", "free_arena_bytes",
        "mmap_bytes", "main_top_free_bytes"};
    REAL(values)[0] = (double)info.arena;
    REAL(values)[1] = (double)info.uordblks;
    REAL(values)[2] = (double)info.fordblks;
    REAL(values)[3] = (double)info.hblkhd;
    REAL(values)[4] = (double)info.keepcost;
    for (int i = 0; i < 5; ++i) SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
    Rf_setAttrib(values, R_NamesSymbol, names);
    UNPROTECT(2);
    return values;
#else
    return R_NilValue;
#endif
}
