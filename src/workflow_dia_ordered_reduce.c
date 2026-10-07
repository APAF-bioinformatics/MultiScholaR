/* LGPL-3.0-or-later. Sequential DIA reduction of explicitly ordered groups.
 * R calls stay on its main thread; long-double totals never cross through R
 * until a complete group's final double is emitted. */
#include <R.h>
#include <Rinternals.h>
#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#if defined(__linux__) && defined(__x86_64__) && LDBL_MANT_DIG == 64 && LDBL_MAX_EXP == 16384
#define GROUP_BYTES ((size_t)67108864)
typedef struct { uint64_t hash; size_t offset; } FormSlot;
typedef struct {
    char *key[3], *arena;
    size_t key_capacity[3], arena_capacity, arena_used, capacity, forms;
    FormSlot *slots;
    long double raw, normalised;
    double ordinal, input_rows, groups, peak_bytes;
    int raw_observed, normalised_observed, source_count, active, finished, failed;
} OrderedReducer;

static size_t state_bytes(const OrderedReducer *s) {
    return sizeof(*s) + s->arena_capacity + s->capacity * sizeof(FormSlot) +
        s->key_capacity[0] + s->key_capacity[1] + s->key_capacity[2];
}
static void fail(OrderedReducer *s, const char *message) {
    s->failed = 1;
    Rf_error("%s", message);
}
static void admit(OrderedReducer *s, size_t extra) {
    size_t current = state_bytes(s);
    if (extra > GROUP_BYTES || current > GROUP_BYTES - extra)
        fail(s, "ordered group state exceeds 64 MiB");
    if ((double)(current + extra) > s->peak_bytes) s->peak_bytes = (double)(current + extra);
}
static void release_state(SEXP pointer) {
    OrderedReducer *s = R_ExternalPtrAddr(pointer);
    if (!s) return;
    for (int key = 0; key < 3; key++) free(s->key[key]);
    free(s->arena);
    free(s->slots);
    free(s);
    R_ClearExternalPtr(pointer);
}
static OrderedReducer *get_state(SEXP pointer) {
    if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrTag(pointer) != Rf_install("dia_ordered_reducer_v1"))
        Rf_error("invalid ordered reducer handle");
    OrderedReducer *s = R_ExternalPtrAddr(pointer);
    if (!s) Rf_error("ordered reducer is closed");
    return s;
}
SEXP workflow_dia_ordered_new(void) {
    OrderedReducer *s = calloc(1, sizeof(*s));
    if (!s) Rf_error("cannot allocate ordered reducer");
    SEXP pointer = PROTECT(R_MakeExternalPtr(s, Rf_install("dia_ordered_reducer_v1"), R_NilValue));
    R_RegisterCFinalizerEx(pointer, release_state, TRUE);
    UNPROTECT(1);
    return pointer;
}
SEXP workflow_dia_ordered_close(SEXP pointer) {
    if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrTag(pointer) != Rf_install("dia_ordered_reducer_v1"))
        Rf_error("invalid ordered reducer handle");
    release_state(pointer);
    return R_NilValue;
}
static void reset_group(OrderedReducer *s) {
    if (s->arena_capacity > 65536) {
        free(s->arena); s->arena = NULL; s->arena_capacity = 0;
    }
    if (s->capacity > 1024) {
        free(s->slots); s->slots = NULL; s->capacity = 0;
    }
    if (s->slots) memset(s->slots, 0, s->capacity * sizeof(FormSlot));
    for (int key = 0; key < 3; key++) if (s->key_capacity[key] > 4096) {
        free(s->key[key]); s->key[key] = NULL; s->key_capacity[key] = 0;
    }
    s->arena_used = s->forms = 0;
    s->raw = s->normalised = 0.0L;
    s->raw_observed = s->normalised_observed = s->source_count = 0;
    s->ordinal = 0;
    s->active = 0;
}
static void copy_keys(OrderedReducer *s, const char **keys) {
    for (int i = 0; i < 3; i++) {
        size_t size = strlen(keys[i]) + 1;
        if (size > s->key_capacity[i]) {
            admit(s, size); /* Include possible realloc overlap conservatively. */
            char *next = realloc(s->key[i], size);
            if (!next) fail(s, "cannot allocate group key");
            s->key[i] = next; s->key_capacity[i] = size;
        }
        memcpy(s->key[i], keys[i], size);
    }
    s->active = 1;
}
static void grow_slots(OrderedReducer *s) {
    size_t capacity = s->capacity ? s->capacity * 2 : 16;
    admit(s, capacity * sizeof(FormSlot));
    FormSlot *slots = calloc(capacity, sizeof(FormSlot));
    if (!slots) fail(s, "cannot allocate group forms");
    for (size_t i = 0; i < s->capacity; i++) if (s->slots[i].hash) {
        size_t at = s->slots[i].hash & (capacity - 1);
        while (slots[at].hash) at = (at + 1) & (capacity - 1);
        slots[at] = s->slots[i];
    }
    free(s->slots); s->slots = slots; s->capacity = capacity;
}
static void add_form(OrderedReducer *s, const char *value) {
    uint64_t hash = UINT64_C(14695981039346656037);
    for (const unsigned char *p = (const unsigned char *)value; *p; p++) {
        hash ^= *p; hash *= UINT64_C(1099511628211);
    }
    if (!hash) hash = 1;
    if (!s->capacity || s->forms * 10 >= s->capacity * 7) grow_slots(s);
    size_t at = hash & (s->capacity - 1);
    while (s->slots[at].hash) {
        if (s->slots[at].hash == hash && !strcmp(s->arena + s->slots[at].offset, value)) return;
        at = (at + 1) & (s->capacity - 1);
    }
    size_t length = strlen(value) + 1;
    if (length > GROUP_BYTES - s->arena_used) fail(s, "ordered group forms exceed 64 MiB");
    size_t needed = s->arena_used + length;
    if (needed > s->arena_capacity) {
        size_t capacity = s->arena_capacity ? s->arena_capacity : 64;
        while (capacity < needed) capacity *= 2;
        admit(s, capacity);
        char *arena = realloc(s->arena, capacity);
        if (!arena) fail(s, "cannot allocate group form values");
        s->arena = arena; s->arena_capacity = capacity;
    }
    memcpy(s->arena + s->arena_used, value, length);
    s->slots[at].hash = hash; s->slots[at].offset = s->arena_used;
    s->arena_used = needed; s->forms++;
}
static void emit_group(OrderedReducer *s, SEXP out, int row) {
    for (int key = 0; key < 3; key++)
        SET_STRING_ELT(VECTOR_ELT(out, key), row, Rf_mkCharCE(s->key[key], CE_UTF8));
    REAL(VECTOR_ELT(out, 3))[row] = !s->raw_observed ? NA_REAL : s->raw > DBL_MAX ? INFINITY : (double)s->raw;
    REAL(VECTOR_ELT(out, 4))[row] = !s->normalised_observed ? NA_REAL : s->normalised > DBL_MAX ? INFINITY : (double)s->normalised;
    INTEGER(VECTOR_ELT(out, 5))[row] = (int)s->forms;
    INTEGER(VECTOR_ELT(out, 6))[row] = s->raw_observed;
    INTEGER(VECTOR_ELT(out, 7))[row] = s->normalised_observed;
    INTEGER(VECTOR_ELT(out, 8))[row] = s->source_count;
    s->groups++;
}
SEXP workflow_dia_ordered_push(SEXP pointer, SEXP data, SEXP final) {
    OrderedReducer *s = get_state(pointer);
    if (s->failed || s->finished) Rf_error("ordered reducer is failed or finished");
    if (TYPEOF(final) != LGLSXP || XLENGTH(final) != 1 || LOGICAL(final)[0] == NA_LOGICAL ||
        TYPEOF(data) != VECSXP || XLENGTH(data) != 7) fail(s, "invalid ordered reducer batch");
    const char *names[] = {"sample", "protein", "peptide", "raw", "normalised", "form", "source_order"};
    SEXP fields = Rf_getAttrib(data, R_NamesSymbol);
    if (TYPEOF(fields) != STRSXP || XLENGTH(fields) != 7) fail(s, "invalid ordered reducer fields");
    R_xlen_t n = XLENGTH(VECTOR_ELT(data, 0));
    if (n > 10000 || (LOGICAL(final)[0] && n)) fail(s, "finish requires an empty batch; transfers are limited to 10000 rows");
    for (int col = 0; col < 7; col++) {
        int type = col < 3 || col == 5 ? STRSXP : REALSXP;
        if (STRING_ELT(fields, col) == NA_STRING || strcmp(CHAR(STRING_ELT(fields, col)), names[col]) ||
            TYPEOF(VECTOR_ELT(data, col)) != type || XLENGTH(VECTOR_ELT(data, col)) != n)
            fail(s, "ordered reducer schema differs");
    }
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 9));
    for (int col = 0; col < 9; col++) SET_VECTOR_ELT(out, col,
        Rf_allocVector(col < 3 ? STRSXP : col < 5 ? REALSXP : INTSXP, n ? n : 1));
    int emitted = 0;
    for (R_xlen_t row = 0; row < n; row++) {
        const char *keys[3];
        for (int col = 0; col < 3; col++) {
            SEXP value = STRING_ELT(VECTOR_ELT(data, col), row);
            if (value == NA_STRING) fail(s, "missing ordered group key");
            keys[col] = Rf_translateCharUTF8(value);
        }
        int comparison = 0;
        if (s->active) for (int col = 0; col < 3 && !comparison; col++) comparison = strcmp(keys[col], s->key[col]);
        if (comparison < 0) fail(s, "group keys are not in ascending binary order");
        if (s->active && comparison > 0) { emit_group(s, out, emitted++); reset_group(s); }
        if (!s->active) copy_keys(s, keys);
        double ordinal = REAL(VECTOR_ELT(data, 6))[row];
        if (!isfinite(ordinal) || ordinal < 1 || ordinal > 9007199254740991.0 ||
            floor(ordinal) != ordinal || ordinal <= s->ordinal) fail(s, "source order is not strictly increasing within group");
        if (s->source_count == INT_MAX) fail(s, "group count exceeds integer range");
        s->ordinal = ordinal; s->source_count++; s->input_rows++;
        double raw = REAL(VECTOR_ELT(data, 3))[row], norm = REAL(VECTOR_ELT(data, 4))[row];
        if (!ISNA(raw)) {
            if (!isfinite(raw) || raw < 0) fail(s, "DIA sum requires finite nonnegative linear quantities");
            s->raw += (long double)raw; s->raw_observed++;
        }
        if (!ISNA(norm)) {
            if (!isfinite(norm) || norm < 0) fail(s, "DIA sum requires finite nonnegative linear quantities");
            s->normalised += (long double)norm; s->normalised_observed++;
        }
        SEXP form = STRING_ELT(VECTOR_ELT(data, 5), row);
        if (form != NA_STRING) add_form(s, Rf_translateCharUTF8(form));
    }
    if (LOGICAL(final)[0]) {
        if (s->active) { emit_group(s, out, emitted++); reset_group(s); }
        s->finished = 1;
    }
    for (int col = 0; col < 9; col++) SET_VECTOR_ELT(out, col, Rf_lengthgets(VECTOR_ELT(out, col), emitted));
    const char *output_names[] = {"sample", "protein", "peptide", "Peptide.RawQuantity", "Peptide.Normalised",
        "peptidoform_count", "raw_observed", "normalised_observed", "source_count"};
    SEXP output_fields = PROTECT(Rf_allocVector(STRSXP, 9));
    for (int i = 0; i < 9; i++) SET_STRING_ELT(output_fields, i, Rf_mkChar(output_names[i]));
    Rf_setAttrib(out, R_NamesSymbol, output_fields);
    SEXP row_names = PROTECT(Rf_allocVector(INTSXP, emitted ? 2 : 0));
    if (emitted) { INTEGER(row_names)[0] = NA_INTEGER; INTEGER(row_names)[1] = -emitted; }
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
    SEXP cls = PROTECT(Rf_mkString("data.frame"));
    Rf_setAttrib(out, R_ClassSymbol, cls);
    if (LOGICAL(final)[0]) {
        SEXP stats = PROTECT(Rf_allocVector(REALSXP, 4));
        REAL(stats)[0] = s->input_rows; REAL(stats)[1] = s->groups;
        REAL(stats)[2] = s->peak_bytes; REAL(stats)[3] = s->input_rows > 0 ? 1 : 0;
        SEXP labels = PROTECT(Rf_allocVector(STRSXP, 4));
        const char *stats_names[] = {"input_rows", "groups", "max_group_state_bytes", "max_active_groups"};
        for (int i = 0; i < 4; i++) SET_STRING_ELT(labels, i, Rf_mkChar(stats_names[i]));
        Rf_setAttrib(stats, R_NamesSymbol, labels);
        Rf_setAttrib(out, Rf_install("ordered_reducer_statistics"), stats);
        UNPROTECT(2);
    }
    UNPROTECT(4);
    return out;
}
#else
SEXP workflow_dia_ordered_new(void) { Rf_error("ordered reducer requires qualified Linux x86_64 long double"); return R_NilValue; }
SEXP workflow_dia_ordered_push(SEXP a, SEXP b, SEXP c) { (void)a; (void)b; (void)c; return workflow_dia_ordered_new(); }
SEXP workflow_dia_ordered_close(SEXP a) { (void)a; return workflow_dia_ordered_new(); }
#endif
