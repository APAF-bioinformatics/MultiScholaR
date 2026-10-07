/* MultiScholaR LGPL-3.0-or-later scientific provider for fastduckdb ABI v1.
 * The callback body is copied byte-for-byte from dia_rsum.c and checked at build.
 * Only registration plumbing differs from the retained experimental extension.
 */
#include <float.h>
#if defined(__linux__) && defined(__x86_64__) && LDBL_MANT_DIG == 64 && LDBL_MAX_EXP == 16384
#include <dlfcn.h>
#include <math.h>
#include "duckdb.h"
static __typeof__(&duckdb_create_data_chunk) native_duckdb_create_data_chunk;
#define duckdb_create_data_chunk native_duckdb_create_data_chunk
static __typeof__(&duckdb_destroy_data_chunk) native_duckdb_destroy_data_chunk;
#define duckdb_destroy_data_chunk native_duckdb_destroy_data_chunk
static __typeof__(&duckdb_data_chunk_set_size) native_duckdb_data_chunk_set_size;
#define duckdb_data_chunk_set_size native_duckdb_data_chunk_set_size
static __typeof__(&duckdb_vector_get_column_type) native_duckdb_vector_get_column_type;
#define duckdb_vector_get_column_type native_duckdb_vector_get_column_type
static __typeof__(&duckdb_create_selection_vector) native_duckdb_create_selection_vector;
#define duckdb_create_selection_vector native_duckdb_create_selection_vector
static __typeof__(&duckdb_destroy_selection_vector) native_duckdb_destroy_selection_vector;
#define duckdb_destroy_selection_vector native_duckdb_destroy_selection_vector
static __typeof__(&duckdb_selection_vector_get_data_ptr) native_duckdb_selection_vector_get_data_ptr;
#define duckdb_selection_vector_get_data_ptr native_duckdb_selection_vector_get_data_ptr
static __typeof__(&duckdb_vector_copy_sel) native_duckdb_vector_copy_sel;
#define duckdb_vector_copy_sel native_duckdb_vector_copy_sel
static __typeof__(&duckdb_vector_size) native_duckdb_vector_size;
#define duckdb_vector_size native_duckdb_vector_size
static __typeof__(&duckdb_create_struct_type) native_duckdb_create_struct_type;
#define duckdb_create_struct_type native_duckdb_create_struct_type
static __typeof__(&duckdb_struct_vector_get_child) native_duckdb_struct_vector_get_child;
#define duckdb_struct_vector_get_child native_duckdb_struct_vector_get_child
static __typeof__(&duckdb_create_list_type) native_duckdb_create_list_type;
#define duckdb_create_list_type native_duckdb_create_list_type
static __typeof__(&duckdb_create_logical_type) native_duckdb_create_logical_type;
#define duckdb_create_logical_type native_duckdb_create_logical_type
static __typeof__(&duckdb_create_scalar_function) native_duckdb_create_scalar_function;
#define duckdb_create_scalar_function native_duckdb_create_scalar_function
static __typeof__(&duckdb_data_chunk_get_size) native_duckdb_data_chunk_get_size;
#define duckdb_data_chunk_get_size native_duckdb_data_chunk_get_size
static __typeof__(&duckdb_data_chunk_get_vector) native_duckdb_data_chunk_get_vector;
#define duckdb_data_chunk_get_vector native_duckdb_data_chunk_get_vector
static __typeof__(&duckdb_destroy_logical_type) native_duckdb_destroy_logical_type;
#define duckdb_destroy_logical_type native_duckdb_destroy_logical_type
static __typeof__(&duckdb_destroy_scalar_function) native_duckdb_destroy_scalar_function;
#define duckdb_destroy_scalar_function native_duckdb_destroy_scalar_function
static __typeof__(&duckdb_list_vector_get_child) native_duckdb_list_vector_get_child;
#define duckdb_list_vector_get_child native_duckdb_list_vector_get_child
static __typeof__(&duckdb_register_scalar_function) native_duckdb_register_scalar_function;
#define duckdb_register_scalar_function native_duckdb_register_scalar_function
static __typeof__(&duckdb_scalar_function_add_parameter) native_duckdb_scalar_function_add_parameter;
#define duckdb_scalar_function_add_parameter native_duckdb_scalar_function_add_parameter
static __typeof__(&duckdb_scalar_function_set_error) native_duckdb_scalar_function_set_error;
#define duckdb_scalar_function_set_error native_duckdb_scalar_function_set_error
static __typeof__(&duckdb_scalar_function_set_function) native_duckdb_scalar_function_set_function;
#define duckdb_scalar_function_set_function native_duckdb_scalar_function_set_function
static __typeof__(&duckdb_scalar_function_set_name) native_duckdb_scalar_function_set_name;
#define duckdb_scalar_function_set_name native_duckdb_scalar_function_set_name
static __typeof__(&duckdb_scalar_function_set_return_type) native_duckdb_scalar_function_set_return_type;
#define duckdb_scalar_function_set_return_type native_duckdb_scalar_function_set_return_type
static __typeof__(&duckdb_scalar_function_set_special_handling) native_duckdb_scalar_function_set_special_handling;
#define duckdb_scalar_function_set_special_handling native_duckdb_scalar_function_set_special_handling
static __typeof__(&duckdb_validity_row_is_valid) native_duckdb_validity_row_is_valid;
#define duckdb_validity_row_is_valid native_duckdb_validity_row_is_valid
static __typeof__(&duckdb_validity_set_row_invalid) native_duckdb_validity_set_row_invalid;
#define duckdb_validity_set_row_invalid native_duckdb_validity_set_row_invalid
static __typeof__(&duckdb_vector_ensure_validity_writable) native_duckdb_vector_ensure_validity_writable;
#define duckdb_vector_ensure_validity_writable native_duckdb_vector_ensure_validity_writable
static __typeof__(&duckdb_vector_get_data) native_duckdb_vector_get_data;
#define duckdb_vector_get_data native_duckdb_vector_get_data
static __typeof__(&duckdb_vector_get_validity) native_duckdb_vector_get_validity;
#define duckdb_vector_get_validity native_duckdb_vector_get_validity

/* One engine reference is retained after the first successful registration;
 * registered callbacks use the pointers below for the life of the process. */
static void *retained_engine = NULL;

static void dia_reset_engine_pointers(void) {
    native_duckdb_create_data_chunk = NULL;
    native_duckdb_destroy_data_chunk = NULL;
    native_duckdb_data_chunk_set_size = NULL;
    native_duckdb_vector_get_column_type = NULL;
    native_duckdb_create_selection_vector = NULL;
    native_duckdb_destroy_selection_vector = NULL;
    native_duckdb_selection_vector_get_data_ptr = NULL;
    native_duckdb_vector_copy_sel = NULL;
    native_duckdb_vector_size = NULL;
    native_duckdb_create_struct_type = NULL;
    native_duckdb_struct_vector_get_child = NULL;
    native_duckdb_create_list_type = NULL;
    native_duckdb_create_logical_type = NULL;
    native_duckdb_create_scalar_function = NULL;
    native_duckdb_data_chunk_get_size = NULL;
    native_duckdb_data_chunk_get_vector = NULL;
    native_duckdb_destroy_logical_type = NULL;
    native_duckdb_destroy_scalar_function = NULL;
    native_duckdb_list_vector_get_child = NULL;
    native_duckdb_register_scalar_function = NULL;
    native_duckdb_scalar_function_add_parameter = NULL;
    native_duckdb_scalar_function_set_error = NULL;
    native_duckdb_scalar_function_set_function = NULL;
    native_duckdb_scalar_function_set_name = NULL;
    native_duckdb_scalar_function_set_return_type = NULL;
    native_duckdb_scalar_function_set_special_handling = NULL;
    native_duckdb_validity_row_is_valid = NULL;
    native_duckdb_validity_set_row_invalid = NULL;
    native_duckdb_vector_ensure_validity_writable = NULL;
    native_duckdb_vector_get_data = NULL;
    native_duckdb_vector_get_validity = NULL;
}

/* Drop this call's engine reference. Pointers are reset only when no earlier
 * registration depends on them; a repeat call resolves identical addresses. */
static const char *dia_abandon_engine(void *library, const char *message) {
    if (!retained_engine) dia_reset_engine_pointers();
    dlclose(library);
    return message;
}

static void dia_rsum(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector lists = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector child = duckdb_list_vector_get_child(lists);
    duckdb_list_entry *entries = duckdb_vector_get_data(lists);
    double *values = duckdb_vector_get_data(child);
    double *result = duckdb_vector_get_data(output);
    uint64_t *list_validity = duckdb_vector_get_validity(lists);
    uint64_t *value_validity = duckdb_vector_get_validity(child);
    duckdb_vector_ensure_validity_writable(output);
    uint64_t *output_validity = duckdb_vector_get_validity(output);
    idx_t rows = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < rows; row++) {
        long double total = 0.0L;
        bool observed = false;
        if (duckdb_validity_row_is_valid(list_validity, row)) {
            idx_t end = entries[row].offset + entries[row].length;
            for (idx_t at = entries[row].offset; at < end; at++) {
                if (!duckdb_validity_row_is_valid(value_validity, at)) continue;
                double value = values[at];
                if (!isfinite(value) || value < 0) {
                    duckdb_scalar_function_set_error(info, "DIA sum requires finite nonnegative linear quantities");
                    return;
                }
                total += (long double)value;
                observed = true;
            }
        }
        if (!observed) {
            duckdb_validity_set_row_invalid(output_validity, row);
        } else {
            result[row] = total > DBL_MAX ? INFINITY : (double)total;
        }
    }
}


/* One ordered list owns both streams; each stream retains the original sum's
 * sequential arithmetic and independent missing-value rule. No partial sums. */
static void dia_paired_rsum(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector lists = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector pairs = duckdb_list_vector_get_child(lists);
    duckdb_list_entry *entries = duckdb_vector_get_data(lists);
    uint64_t *list_validity = duckdb_vector_get_validity(lists);
    uint64_t *pair_validity = duckdb_vector_get_validity(pairs);
    idx_t rows = duckdb_data_chunk_get_size(input);
    for (idx_t column = 0; column < 2; column++) {
        duckdb_vector child = duckdb_struct_vector_get_child(pairs, column);
        duckdb_vector destination = duckdb_struct_vector_get_child(output, column);
        double *values = duckdb_vector_get_data(child);
        double *result = duckdb_vector_get_data(destination);
        uint64_t *value_validity = duckdb_vector_get_validity(child);
        duckdb_vector_ensure_validity_writable(destination);
        uint64_t *output_validity = duckdb_vector_get_validity(destination);
        for (idx_t row = 0; row < rows; row++) {
            long double total = 0.0L;
            bool observed = false;
            if (duckdb_validity_row_is_valid(list_validity, row)) {
                idx_t end = entries[row].offset + entries[row].length;
                for (idx_t at = entries[row].offset; at < end; at++) {
                    if (!duckdb_validity_row_is_valid(pair_validity, at) ||
                        !duckdb_validity_row_is_valid(value_validity, at)) continue;
                    double value = values[at];
                    if (!isfinite(value) || value < 0) {
                        duckdb_scalar_function_set_error(info, "DIA sum requires finite nonnegative linear quantities");
                        return;
                    }
                    total += (long double)value;
                    observed = true;
                }
            }
            if (!observed) duckdb_validity_set_row_invalid(output_validity, row);
            else result[row] = total > DBL_MAX ? INFINITY : (double)total;
        }
    }
}

/* DuckDB 1.5.5 flattens a constant LIST's entries but not its child. Only
 * shared-entry inputs can arrive through that branch. Ordinary grouped lists
 * have distinct entries and recursively flattened children; keep their path.
 * Copy one shared list using the public API, which expands child encodings,
 * then reuse that normalized entry for every row. Arithmetic is unchanged. */
static void dia_apply_flat_list(duckdb_function_info info, duckdb_data_chunk input,
                                duckdb_vector output, duckdb_scalar_function_t callback) {
    idx_t rows = duckdb_data_chunk_get_size(input);
    duckdb_vector lists = duckdb_data_chunk_get_vector(input, 0);
    duckdb_list_entry *entries = duckdb_vector_get_data(lists);
    uint64_t *validity = duckdb_vector_get_validity(lists);
    if (!rows || !duckdb_validity_row_is_valid(validity, 0) || !entries[0].length) {
        callback(info, input, output);
        return;
    }
    for (idx_t row = 1; row < rows; row++) {
        if (!duckdb_validity_row_is_valid(validity, row) ||
            entries[row].offset != entries[0].offset || entries[row].length != entries[0].length) {
            callback(info, input, output);
            return;
        }
    }
    if (rows > duckdb_vector_size()) {
        duckdb_scalar_function_set_error(info, "DIA scalar chunk exceeds the engine vector size");
        return;
    }
    duckdb_logical_type type = duckdb_vector_get_column_type(lists);
    duckdb_data_chunk normalized = duckdb_create_data_chunk(&type, 1);
    duckdb_destroy_logical_type(&type);
    if (!normalized) {
        duckdb_scalar_function_set_error(info, "cannot allocate normalized DIA list input");
        return;
    }
    duckdb_selection_vector selection = duckdb_create_selection_vector(1);
    duckdb_selection_vector_get_data_ptr(selection)[0] = 0;
    duckdb_vector flat = duckdb_data_chunk_get_vector(normalized, 0);
    duckdb_vector_copy_sel(lists, flat, selection, 1, 0, 0);
    duckdb_destroy_selection_vector(selection);
    duckdb_list_entry *flat_entries = duckdb_vector_get_data(flat);
    for (idx_t row = 1; row < rows; row++) flat_entries[row] = flat_entries[0];
    duckdb_data_chunk_set_size(normalized, rows);
    callback(info, normalized, output);
    duckdb_destroy_data_chunk(&normalized);
}

static void dia_rsum_flat_input(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    dia_apply_flat_list(info, input, output, dia_rsum);
}

static void dia_paired_rsum_flat_input(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    dia_apply_flat_list(info, input, output, dia_paired_rsum);
}

const char *workflowRegisterDIAOrderedSum(void *connection, const char *engine_library) {
    void *library = dlopen(engine_library, RTLD_NOW | RTLD_LOCAL);
    if (!library) return "cannot retain the scoped DuckDB engine";
    if (retained_engine && library != retained_engine) {
        dlclose(library);
        return "the scoped DuckDB engine changed within this process";
    }
    native_duckdb_create_struct_type = (__typeof__(native_duckdb_create_struct_type))dlsym(library, "duckdb_create_struct_type");
    if (!native_duckdb_create_struct_type) return dia_abandon_engine(library, "missing public C API: duckdb_create_struct_type");
    native_duckdb_struct_vector_get_child = (__typeof__(native_duckdb_struct_vector_get_child))dlsym(library, "duckdb_struct_vector_get_child");
    if (!native_duckdb_struct_vector_get_child) return dia_abandon_engine(library, "missing public C API: duckdb_struct_vector_get_child");
    native_duckdb_create_list_type = (__typeof__(native_duckdb_create_list_type))dlsym(library, "duckdb_create_list_type");
    if (!native_duckdb_create_list_type) return dia_abandon_engine(library, "missing public C API: duckdb_create_list_type");
    native_duckdb_create_logical_type = (__typeof__(native_duckdb_create_logical_type))dlsym(library, "duckdb_create_logical_type");
    if (!native_duckdb_create_logical_type) return dia_abandon_engine(library, "missing public C API: duckdb_create_logical_type");
    native_duckdb_create_scalar_function = (__typeof__(native_duckdb_create_scalar_function))dlsym(library, "duckdb_create_scalar_function");
    if (!native_duckdb_create_scalar_function) return dia_abandon_engine(library, "missing public C API: duckdb_create_scalar_function");
    native_duckdb_data_chunk_get_size = (__typeof__(native_duckdb_data_chunk_get_size))dlsym(library, "duckdb_data_chunk_get_size");
    if (!native_duckdb_data_chunk_get_size) return dia_abandon_engine(library, "missing public C API: duckdb_data_chunk_get_size");
    native_duckdb_data_chunk_get_vector = (__typeof__(native_duckdb_data_chunk_get_vector))dlsym(library, "duckdb_data_chunk_get_vector");
    if (!native_duckdb_data_chunk_get_vector) return dia_abandon_engine(library, "missing public C API: duckdb_data_chunk_get_vector");
    native_duckdb_destroy_logical_type = (__typeof__(native_duckdb_destroy_logical_type))dlsym(library, "duckdb_destroy_logical_type");
    if (!native_duckdb_destroy_logical_type) return dia_abandon_engine(library, "missing public C API: duckdb_destroy_logical_type");
    native_duckdb_destroy_scalar_function = (__typeof__(native_duckdb_destroy_scalar_function))dlsym(library, "duckdb_destroy_scalar_function");
    if (!native_duckdb_destroy_scalar_function) return dia_abandon_engine(library, "missing public C API: duckdb_destroy_scalar_function");
    native_duckdb_list_vector_get_child = (__typeof__(native_duckdb_list_vector_get_child))dlsym(library, "duckdb_list_vector_get_child");
    if (!native_duckdb_list_vector_get_child) return dia_abandon_engine(library, "missing public C API: duckdb_list_vector_get_child");
    native_duckdb_register_scalar_function = (__typeof__(native_duckdb_register_scalar_function))dlsym(library, "duckdb_register_scalar_function");
    if (!native_duckdb_register_scalar_function) return dia_abandon_engine(library, "missing public C API: duckdb_register_scalar_function");
    native_duckdb_scalar_function_add_parameter = (__typeof__(native_duckdb_scalar_function_add_parameter))dlsym(library, "duckdb_scalar_function_add_parameter");
    if (!native_duckdb_scalar_function_add_parameter) return dia_abandon_engine(library, "missing public C API: duckdb_scalar_function_add_parameter");
    native_duckdb_scalar_function_set_error = (__typeof__(native_duckdb_scalar_function_set_error))dlsym(library, "duckdb_scalar_function_set_error");
    if (!native_duckdb_scalar_function_set_error) return dia_abandon_engine(library, "missing public C API: duckdb_scalar_function_set_error");
    native_duckdb_scalar_function_set_function = (__typeof__(native_duckdb_scalar_function_set_function))dlsym(library, "duckdb_scalar_function_set_function");
    if (!native_duckdb_scalar_function_set_function) return dia_abandon_engine(library, "missing public C API: duckdb_scalar_function_set_function");
    native_duckdb_scalar_function_set_name = (__typeof__(native_duckdb_scalar_function_set_name))dlsym(library, "duckdb_scalar_function_set_name");
    if (!native_duckdb_scalar_function_set_name) return dia_abandon_engine(library, "missing public C API: duckdb_scalar_function_set_name");
    native_duckdb_scalar_function_set_return_type = (__typeof__(native_duckdb_scalar_function_set_return_type))dlsym(library, "duckdb_scalar_function_set_return_type");
    if (!native_duckdb_scalar_function_set_return_type) return dia_abandon_engine(library, "missing public C API: duckdb_scalar_function_set_return_type");
    native_duckdb_scalar_function_set_special_handling = (__typeof__(native_duckdb_scalar_function_set_special_handling))dlsym(library, "duckdb_scalar_function_set_special_handling");
    if (!native_duckdb_scalar_function_set_special_handling) return dia_abandon_engine(library, "missing public C API: duckdb_scalar_function_set_special_handling");
    native_duckdb_validity_row_is_valid = (__typeof__(native_duckdb_validity_row_is_valid))dlsym(library, "duckdb_validity_row_is_valid");
    if (!native_duckdb_validity_row_is_valid) return dia_abandon_engine(library, "missing public C API: duckdb_validity_row_is_valid");
    native_duckdb_validity_set_row_invalid = (__typeof__(native_duckdb_validity_set_row_invalid))dlsym(library, "duckdb_validity_set_row_invalid");
    if (!native_duckdb_validity_set_row_invalid) return dia_abandon_engine(library, "missing public C API: duckdb_validity_set_row_invalid");
    native_duckdb_vector_ensure_validity_writable = (__typeof__(native_duckdb_vector_ensure_validity_writable))dlsym(library, "duckdb_vector_ensure_validity_writable");
    if (!native_duckdb_vector_ensure_validity_writable) return dia_abandon_engine(library, "missing public C API: duckdb_vector_ensure_validity_writable");
    native_duckdb_vector_get_data = (__typeof__(native_duckdb_vector_get_data))dlsym(library, "duckdb_vector_get_data");
    if (!native_duckdb_vector_get_data) return dia_abandon_engine(library, "missing public C API: duckdb_vector_get_data");
    native_duckdb_vector_get_validity = (__typeof__(native_duckdb_vector_get_validity))dlsym(library, "duckdb_vector_get_validity");
    if (!native_duckdb_vector_get_validity) return dia_abandon_engine(library, "missing public C API: duckdb_vector_get_validity");
    native_duckdb_create_data_chunk = (__typeof__(native_duckdb_create_data_chunk))dlsym(library, "duckdb_create_data_chunk");
    if (!native_duckdb_create_data_chunk) return dia_abandon_engine(library, "missing public C API: duckdb_create_data_chunk");
    native_duckdb_destroy_data_chunk = (__typeof__(native_duckdb_destroy_data_chunk))dlsym(library, "duckdb_destroy_data_chunk");
    if (!native_duckdb_destroy_data_chunk) return dia_abandon_engine(library, "missing public C API: duckdb_destroy_data_chunk");
    native_duckdb_data_chunk_set_size = (__typeof__(native_duckdb_data_chunk_set_size))dlsym(library, "duckdb_data_chunk_set_size");
    if (!native_duckdb_data_chunk_set_size) return dia_abandon_engine(library, "missing public C API: duckdb_data_chunk_set_size");
    native_duckdb_vector_get_column_type = (__typeof__(native_duckdb_vector_get_column_type))dlsym(library, "duckdb_vector_get_column_type");
    if (!native_duckdb_vector_get_column_type) return dia_abandon_engine(library, "missing public C API: duckdb_vector_get_column_type");
    native_duckdb_create_selection_vector = (__typeof__(native_duckdb_create_selection_vector))dlsym(library, "duckdb_create_selection_vector");
    if (!native_duckdb_create_selection_vector) return dia_abandon_engine(library, "missing public C API: duckdb_create_selection_vector");
    native_duckdb_destroy_selection_vector = (__typeof__(native_duckdb_destroy_selection_vector))dlsym(library, "duckdb_destroy_selection_vector");
    if (!native_duckdb_destroy_selection_vector) return dia_abandon_engine(library, "missing public C API: duckdb_destroy_selection_vector");
    native_duckdb_selection_vector_get_data_ptr = (__typeof__(native_duckdb_selection_vector_get_data_ptr))dlsym(library, "duckdb_selection_vector_get_data_ptr");
    if (!native_duckdb_selection_vector_get_data_ptr) return dia_abandon_engine(library, "missing public C API: duckdb_selection_vector_get_data_ptr");
    native_duckdb_vector_copy_sel = (__typeof__(native_duckdb_vector_copy_sel))dlsym(library, "duckdb_vector_copy_sel");
    if (!native_duckdb_vector_copy_sel) return dia_abandon_engine(library, "missing public C API: duckdb_vector_copy_sel");
    native_duckdb_vector_size = (__typeof__(native_duckdb_vector_size))dlsym(library, "duckdb_vector_size");
    if (!native_duckdb_vector_size) return dia_abandon_engine(library, "missing public C API: duckdb_vector_size");
    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_logical_type quantity = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type list = duckdb_create_list_type(quantity);
    duckdb_scalar_function_set_name(function, "dia_r_sum");
    duckdb_scalar_function_add_parameter(function, list);
    duckdb_scalar_function_set_return_type(function, quantity);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, dia_rsum_flat_input);
    duckdb_state status = duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);
    duckdb_destroy_logical_type(&list);
    if (status == DuckDBSuccess) {
        duckdb_logical_type children[] = {quantity, quantity};
        const char *names[] = {"raw", "normalised"};
        duckdb_logical_type pair = duckdb_create_struct_type(children, names, 2);
        list = duckdb_create_list_type(pair);
        function = duckdb_create_scalar_function();
        duckdb_scalar_function_set_name(function, "dia_paired_r_sum");
        duckdb_scalar_function_add_parameter(function, list);
        duckdb_scalar_function_set_return_type(function, pair);
        duckdb_scalar_function_set_special_handling(function);
        duckdb_scalar_function_set_function(function, dia_paired_rsum_flat_input);
        status = duckdb_register_scalar_function(connection, function);
        duckdb_destroy_scalar_function(&function);
        duckdb_destroy_logical_type(&list);
        duckdb_destroy_logical_type(&pair);
    }
    duckdb_destroy_logical_type(&quantity);
    if (status != DuckDBSuccess) return dia_abandon_engine(library, "DIA ordered sum registration failed");
    if (retained_engine) dlclose(library); else retained_engine = library;
    return NULL;
}
#else
const char *workflowRegisterDIAOrderedSum(void *connection, const char *engine_library) {
    (void)connection; (void)engine_library;
    return "DIA native sum requires the qualified Linux x86_64 long-double runtime";
}
#endif
