/* DIA-025 experimental DuckDB scalar: source-ordered, R-compatible double sum.
 * Build as an isolated C API extension, never as a production package default.
 * The caller supplies LIST(DOUBLE ORDER BY immutable source ordinal).
 * No R objects or callbacks are touched on DuckDB execution threads.
 */
#define DUCKDB_EXTENSION_NAME dia_rsum
#include "duckdb_extension.h"
#include <float.h>
#include <math.h>

#if FLT_RADIX != 2 || LDBL_MANT_DIG != 64 || LDBL_MAX_EXP != 16384
#error "DIA experiment requires the frozen R build's 64-bit long-double accumulator"
#endif

DUCKDB_EXTENSION_EXTERN

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

DUCKDB_EXTENSION_ENTRYPOINT(duckdb_connection connection, duckdb_extension_info info,
                           struct duckdb_extension_access *access) {
    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_logical_type quantity = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type list = duckdb_create_list_type(quantity);
    duckdb_scalar_function_set_name(function, "dia_r_sum");
    duckdb_scalar_function_add_parameter(function, list);
    duckdb_scalar_function_set_return_type(function, quantity);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, dia_rsum);
    duckdb_state status = duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);
    duckdb_destroy_logical_type(&list);
    duckdb_destroy_logical_type(&quantity);
    if (status != DuckDBSuccess) access->set_error(info, "DIA sum registration failed");
    return status == DuckDBSuccess;
}
