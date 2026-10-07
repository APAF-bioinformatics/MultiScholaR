/* PERF-108 feasibility only: public C API, no R driver internals. */
#include <R.h>
#include <Rinternals.h>
#include <dlfcn.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "duckdb.h"
#if FLT_RADIX != 2 || LDBL_MANT_DIG != 64 || LDBL_MAX_EXP != 16384
#error "Requires the accepted 64-bit long-double accumulator"
#endif
static __typeof__(&duckdb_close) p_duckdb_close;
#define duckdb_close p_duckdb_close
static __typeof__(&duckdb_connect) p_duckdb_connect;
#define duckdb_connect p_duckdb_connect
static __typeof__(&duckdb_create_config) p_duckdb_create_config;
#define duckdb_create_config p_duckdb_create_config
static __typeof__(&duckdb_create_list_type) p_duckdb_create_list_type;
#define duckdb_create_list_type p_duckdb_create_list_type
static __typeof__(&duckdb_create_logical_type) p_duckdb_create_logical_type;
#define duckdb_create_logical_type p_duckdb_create_logical_type
static __typeof__(&duckdb_create_scalar_function) p_duckdb_create_scalar_function;
#define duckdb_create_scalar_function p_duckdb_create_scalar_function
static __typeof__(&duckdb_data_chunk_get_size) p_duckdb_data_chunk_get_size;
#define duckdb_data_chunk_get_size p_duckdb_data_chunk_get_size
static __typeof__(&duckdb_data_chunk_get_vector) p_duckdb_data_chunk_get_vector;
#define duckdb_data_chunk_get_vector p_duckdb_data_chunk_get_vector
static __typeof__(&duckdb_destroy_config) p_duckdb_destroy_config;
#define duckdb_destroy_config p_duckdb_destroy_config
static __typeof__(&duckdb_destroy_data_chunk) p_duckdb_destroy_data_chunk;
#define duckdb_destroy_data_chunk p_duckdb_destroy_data_chunk
static __typeof__(&duckdb_destroy_logical_type) p_duckdb_destroy_logical_type;
#define duckdb_destroy_logical_type p_duckdb_destroy_logical_type
static __typeof__(&duckdb_destroy_prepare) p_duckdb_destroy_prepare;
#define duckdb_destroy_prepare p_duckdb_destroy_prepare
static __typeof__(&duckdb_destroy_result) p_duckdb_destroy_result;
#define duckdb_destroy_result p_duckdb_destroy_result
static __typeof__(&duckdb_destroy_scalar_function) p_duckdb_destroy_scalar_function;
#define duckdb_destroy_scalar_function p_duckdb_destroy_scalar_function
static __typeof__(&duckdb_disconnect) p_duckdb_disconnect;
#define duckdb_disconnect p_duckdb_disconnect
static __typeof__(&duckdb_execute_prepared_streaming) p_duckdb_execute_prepared_streaming;
#define duckdb_execute_prepared_streaming p_duckdb_execute_prepared_streaming
static __typeof__(&duckdb_result_is_streaming) p_duckdb_result_is_streaming;
#define duckdb_result_is_streaming p_duckdb_result_is_streaming
static __typeof__(&duckdb_fetch_chunk) p_duckdb_fetch_chunk;
#define duckdb_fetch_chunk p_duckdb_fetch_chunk
static __typeof__(&duckdb_free) p_duckdb_free;
#define duckdb_free p_duckdb_free
static __typeof__(&duckdb_library_version) p_duckdb_library_version;
#define duckdb_library_version p_duckdb_library_version
static __typeof__(&duckdb_list_vector_get_child) p_duckdb_list_vector_get_child;
#define duckdb_list_vector_get_child p_duckdb_list_vector_get_child
static __typeof__(&duckdb_open_ext) p_duckdb_open_ext;
#define duckdb_open_ext p_duckdb_open_ext
static __typeof__(&duckdb_prepare) p_duckdb_prepare;
#define duckdb_prepare p_duckdb_prepare
static __typeof__(&duckdb_query) p_duckdb_query;
#define duckdb_query p_duckdb_query
static __typeof__(&duckdb_register_scalar_function) p_duckdb_register_scalar_function;
#define duckdb_register_scalar_function p_duckdb_register_scalar_function
static __typeof__(&duckdb_result_error) p_duckdb_result_error;
#define duckdb_result_error p_duckdb_result_error
static __typeof__(&duckdb_scalar_function_add_parameter) p_duckdb_scalar_function_add_parameter;
#define duckdb_scalar_function_add_parameter p_duckdb_scalar_function_add_parameter
static __typeof__(&duckdb_scalar_function_set_error) p_duckdb_scalar_function_set_error;
#define duckdb_scalar_function_set_error p_duckdb_scalar_function_set_error
static __typeof__(&duckdb_scalar_function_set_function) p_duckdb_scalar_function_set_function;
#define duckdb_scalar_function_set_function p_duckdb_scalar_function_set_function
static __typeof__(&duckdb_scalar_function_set_name) p_duckdb_scalar_function_set_name;
#define duckdb_scalar_function_set_name p_duckdb_scalar_function_set_name
static __typeof__(&duckdb_scalar_function_set_return_type) p_duckdb_scalar_function_set_return_type;
#define duckdb_scalar_function_set_return_type p_duckdb_scalar_function_set_return_type
static __typeof__(&duckdb_scalar_function_set_special_handling) p_duckdb_scalar_function_set_special_handling;
#define duckdb_scalar_function_set_special_handling p_duckdb_scalar_function_set_special_handling
static __typeof__(&duckdb_set_config) p_duckdb_set_config;
#define duckdb_set_config p_duckdb_set_config
static __typeof__(&duckdb_validity_row_is_valid) p_duckdb_validity_row_is_valid;
#define duckdb_validity_row_is_valid p_duckdb_validity_row_is_valid
static __typeof__(&duckdb_validity_set_row_invalid) p_duckdb_validity_set_row_invalid;
#define duckdb_validity_set_row_invalid p_duckdb_validity_set_row_invalid
static __typeof__(&duckdb_value_double) p_duckdb_value_double;
#define duckdb_value_double p_duckdb_value_double
static __typeof__(&duckdb_value_is_null) p_duckdb_value_is_null;
#define duckdb_value_is_null p_duckdb_value_is_null
static __typeof__(&duckdb_vector_ensure_validity_writable) p_duckdb_vector_ensure_validity_writable;
#define duckdb_vector_ensure_validity_writable p_duckdb_vector_ensure_validity_writable
static __typeof__(&duckdb_vector_get_data) p_duckdb_vector_get_data;
#define duckdb_vector_get_data p_duckdb_vector_get_data
static __typeof__(&duckdb_vector_get_validity) p_duckdb_vector_get_validity;
#define duckdb_vector_get_validity p_duckdb_vector_get_validity

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

static duckdb_state register_sum(duckdb_connection connection) {
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
    return status;
}

SEXP fastduckdb_public_api_probe(SEXP library) {
    void *handle = dlopen(CHAR(STRING_ELT(library, 0)), RTLD_NOW | RTLD_LOCAL);
    if (!handle) Rf_error("could not retain DuckDB DLL: %s", dlerror());
    char failure_message[512] = "";
    duckdb_database database = NULL;
    duckdb_connection connection = NULL;
    duckdb_config config = NULL;
    duckdb_result result = {0};
    duckdb_prepared_statement prepared = NULL;
    duckdb_data_chunk chunk = NULL;
    char *open_error = NULL;
    int checks = 0;
#define REQUIRE(condition, message) do { if (!(condition)) { snprintf(failure_message, sizeof(failure_message), "%s", message); goto cleanup; } } while (0)
    p_duckdb_close = (__typeof__(p_duckdb_close))dlsym(handle, "duckdb_close");
    REQUIRE(p_duckdb_close != NULL, "missing public C symbol duckdb_close");
    p_duckdb_connect = (__typeof__(p_duckdb_connect))dlsym(handle, "duckdb_connect");
    REQUIRE(p_duckdb_connect != NULL, "missing public C symbol duckdb_connect");
    p_duckdb_create_config = (__typeof__(p_duckdb_create_config))dlsym(handle, "duckdb_create_config");
    REQUIRE(p_duckdb_create_config != NULL, "missing public C symbol duckdb_create_config");
    p_duckdb_create_list_type = (__typeof__(p_duckdb_create_list_type))dlsym(handle, "duckdb_create_list_type");
    REQUIRE(p_duckdb_create_list_type != NULL, "missing public C symbol duckdb_create_list_type");
    p_duckdb_create_logical_type = (__typeof__(p_duckdb_create_logical_type))dlsym(handle, "duckdb_create_logical_type");
    REQUIRE(p_duckdb_create_logical_type != NULL, "missing public C symbol duckdb_create_logical_type");
    p_duckdb_create_scalar_function = (__typeof__(p_duckdb_create_scalar_function))dlsym(handle, "duckdb_create_scalar_function");
    REQUIRE(p_duckdb_create_scalar_function != NULL, "missing public C symbol duckdb_create_scalar_function");
    p_duckdb_data_chunk_get_size = (__typeof__(p_duckdb_data_chunk_get_size))dlsym(handle, "duckdb_data_chunk_get_size");
    REQUIRE(p_duckdb_data_chunk_get_size != NULL, "missing public C symbol duckdb_data_chunk_get_size");
    p_duckdb_data_chunk_get_vector = (__typeof__(p_duckdb_data_chunk_get_vector))dlsym(handle, "duckdb_data_chunk_get_vector");
    REQUIRE(p_duckdb_data_chunk_get_vector != NULL, "missing public C symbol duckdb_data_chunk_get_vector");
    p_duckdb_destroy_config = (__typeof__(p_duckdb_destroy_config))dlsym(handle, "duckdb_destroy_config");
    REQUIRE(p_duckdb_destroy_config != NULL, "missing public C symbol duckdb_destroy_config");
    p_duckdb_destroy_data_chunk = (__typeof__(p_duckdb_destroy_data_chunk))dlsym(handle, "duckdb_destroy_data_chunk");
    REQUIRE(p_duckdb_destroy_data_chunk != NULL, "missing public C symbol duckdb_destroy_data_chunk");
    p_duckdb_destroy_logical_type = (__typeof__(p_duckdb_destroy_logical_type))dlsym(handle, "duckdb_destroy_logical_type");
    REQUIRE(p_duckdb_destroy_logical_type != NULL, "missing public C symbol duckdb_destroy_logical_type");
    p_duckdb_destroy_prepare = (__typeof__(p_duckdb_destroy_prepare))dlsym(handle, "duckdb_destroy_prepare");
    REQUIRE(p_duckdb_destroy_prepare != NULL, "missing public C symbol duckdb_destroy_prepare");
    p_duckdb_destroy_result = (__typeof__(p_duckdb_destroy_result))dlsym(handle, "duckdb_destroy_result");
    REQUIRE(p_duckdb_destroy_result != NULL, "missing public C symbol duckdb_destroy_result");
    p_duckdb_destroy_scalar_function = (__typeof__(p_duckdb_destroy_scalar_function))dlsym(handle, "duckdb_destroy_scalar_function");
    REQUIRE(p_duckdb_destroy_scalar_function != NULL, "missing public C symbol duckdb_destroy_scalar_function");
    p_duckdb_disconnect = (__typeof__(p_duckdb_disconnect))dlsym(handle, "duckdb_disconnect");
    REQUIRE(p_duckdb_disconnect != NULL, "missing public C symbol duckdb_disconnect");
    p_duckdb_execute_prepared_streaming = (__typeof__(p_duckdb_execute_prepared_streaming))dlsym(handle, "duckdb_execute_prepared_streaming");
    REQUIRE(p_duckdb_execute_prepared_streaming != NULL, "missing public C symbol duckdb_execute_prepared_streaming");
    p_duckdb_fetch_chunk = (__typeof__(p_duckdb_fetch_chunk))dlsym(handle, "duckdb_fetch_chunk");
    REQUIRE(p_duckdb_fetch_chunk != NULL, "missing public C symbol duckdb_fetch_chunk");
    p_duckdb_free = (__typeof__(p_duckdb_free))dlsym(handle, "duckdb_free");
    REQUIRE(p_duckdb_free != NULL, "missing public C symbol duckdb_free");
    p_duckdb_library_version = (__typeof__(p_duckdb_library_version))dlsym(handle, "duckdb_library_version");
    REQUIRE(p_duckdb_library_version != NULL, "missing public C symbol duckdb_library_version");
    p_duckdb_list_vector_get_child = (__typeof__(p_duckdb_list_vector_get_child))dlsym(handle, "duckdb_list_vector_get_child");
    REQUIRE(p_duckdb_list_vector_get_child != NULL, "missing public C symbol duckdb_list_vector_get_child");
    p_duckdb_open_ext = (__typeof__(p_duckdb_open_ext))dlsym(handle, "duckdb_open_ext");
    REQUIRE(p_duckdb_open_ext != NULL, "missing public C symbol duckdb_open_ext");
    p_duckdb_prepare = (__typeof__(p_duckdb_prepare))dlsym(handle, "duckdb_prepare");
    REQUIRE(p_duckdb_prepare != NULL, "missing public C symbol duckdb_prepare");
    p_duckdb_query = (__typeof__(p_duckdb_query))dlsym(handle, "duckdb_query");
    REQUIRE(p_duckdb_query != NULL, "missing public C symbol duckdb_query");
    p_duckdb_register_scalar_function = (__typeof__(p_duckdb_register_scalar_function))dlsym(handle, "duckdb_register_scalar_function");
    REQUIRE(p_duckdb_register_scalar_function != NULL, "missing public C symbol duckdb_register_scalar_function");
    p_duckdb_result_error = (__typeof__(p_duckdb_result_error))dlsym(handle, "duckdb_result_error");
    REQUIRE(p_duckdb_result_error != NULL, "missing public C symbol duckdb_result_error");
    p_duckdb_scalar_function_add_parameter = (__typeof__(p_duckdb_scalar_function_add_parameter))dlsym(handle, "duckdb_scalar_function_add_parameter");
    REQUIRE(p_duckdb_scalar_function_add_parameter != NULL, "missing public C symbol duckdb_scalar_function_add_parameter");
    p_duckdb_scalar_function_set_error = (__typeof__(p_duckdb_scalar_function_set_error))dlsym(handle, "duckdb_scalar_function_set_error");
    REQUIRE(p_duckdb_scalar_function_set_error != NULL, "missing public C symbol duckdb_scalar_function_set_error");
    p_duckdb_scalar_function_set_function = (__typeof__(p_duckdb_scalar_function_set_function))dlsym(handle, "duckdb_scalar_function_set_function");
    REQUIRE(p_duckdb_scalar_function_set_function != NULL, "missing public C symbol duckdb_scalar_function_set_function");
    p_duckdb_scalar_function_set_name = (__typeof__(p_duckdb_scalar_function_set_name))dlsym(handle, "duckdb_scalar_function_set_name");
    REQUIRE(p_duckdb_scalar_function_set_name != NULL, "missing public C symbol duckdb_scalar_function_set_name");
    p_duckdb_scalar_function_set_return_type = (__typeof__(p_duckdb_scalar_function_set_return_type))dlsym(handle, "duckdb_scalar_function_set_return_type");
    REQUIRE(p_duckdb_scalar_function_set_return_type != NULL, "missing public C symbol duckdb_scalar_function_set_return_type");
    p_duckdb_scalar_function_set_special_handling = (__typeof__(p_duckdb_scalar_function_set_special_handling))dlsym(handle, "duckdb_scalar_function_set_special_handling");
    REQUIRE(p_duckdb_scalar_function_set_special_handling != NULL, "missing public C symbol duckdb_scalar_function_set_special_handling");
    p_duckdb_set_config = (__typeof__(p_duckdb_set_config))dlsym(handle, "duckdb_set_config");
    REQUIRE(p_duckdb_set_config != NULL, "missing public C symbol duckdb_set_config");
    p_duckdb_validity_row_is_valid = (__typeof__(p_duckdb_validity_row_is_valid))dlsym(handle, "duckdb_validity_row_is_valid");
    REQUIRE(p_duckdb_validity_row_is_valid != NULL, "missing public C symbol duckdb_validity_row_is_valid");
    p_duckdb_validity_set_row_invalid = (__typeof__(p_duckdb_validity_set_row_invalid))dlsym(handle, "duckdb_validity_set_row_invalid");
    REQUIRE(p_duckdb_validity_set_row_invalid != NULL, "missing public C symbol duckdb_validity_set_row_invalid");
    p_duckdb_value_double = (__typeof__(p_duckdb_value_double))dlsym(handle, "duckdb_value_double");
    REQUIRE(p_duckdb_value_double != NULL, "missing public C symbol duckdb_value_double");
    p_duckdb_value_is_null = (__typeof__(p_duckdb_value_is_null))dlsym(handle, "duckdb_value_is_null");
    REQUIRE(p_duckdb_value_is_null != NULL, "missing public C symbol duckdb_value_is_null");
    p_duckdb_vector_ensure_validity_writable = (__typeof__(p_duckdb_vector_ensure_validity_writable))dlsym(handle, "duckdb_vector_ensure_validity_writable");
    REQUIRE(p_duckdb_vector_ensure_validity_writable != NULL, "missing public C symbol duckdb_vector_ensure_validity_writable");
    p_duckdb_vector_get_data = (__typeof__(p_duckdb_vector_get_data))dlsym(handle, "duckdb_vector_get_data");
    REQUIRE(p_duckdb_vector_get_data != NULL, "missing public C symbol duckdb_vector_get_data");
    p_duckdb_vector_get_validity = (__typeof__(p_duckdb_vector_get_validity))dlsym(handle, "duckdb_vector_get_validity");
    REQUIRE(p_duckdb_vector_get_validity != NULL, "missing public C symbol duckdb_vector_get_validity");
    p_duckdb_result_is_streaming = (__typeof__(p_duckdb_result_is_streaming))dlsym(handle, "duckdb_result_is_streaming");
    REQUIRE(p_duckdb_result_is_streaming != NULL, "missing streaming inspection symbol");
    REQUIRE(duckdb_create_config(&config) == DuckDBSuccess, "config failed");
    REQUIRE(duckdb_set_config(config, "threads", "8") == DuckDBSuccess, "threads failed");
    REQUIRE(duckdb_set_config(config, "memory_limit", "12GiB") == DuckDBSuccess, "memory limit failed");
    REQUIRE(duckdb_set_config(config, "allow_unsigned_extensions", "false") == DuckDBSuccess, "unsigned policy failed");
    REQUIRE(duckdb_set_config(config, "temp_directory", "") == DuckDBSuccess, "temp directory failed");
    REQUIRE(duckdb_set_config(config, "max_temp_directory_size", "0B") == DuckDBSuccess, "no spill failed");
    REQUIRE(duckdb_open_ext(NULL, &database, config, &open_error) == DuckDBSuccess, "database open failed");
    REQUIRE(duckdb_connect(database, &connection) == DuckDBSuccess, "connect failed");
    REQUIRE(register_sum(connection) == DuckDBSuccess, "direct registration failed");
    REQUIRE(duckdb_query(connection, "SELECT dia_r_sum([1e16,1,1]::DOUBLE[])", &result) == DuckDBSuccess, "threshold query failed");
    REQUIRE(duckdb_value_double(&result, 0, 0) == 10000000000000002.0, "threshold differs");
    checks++;
    duckdb_destroy_result(&result);
    REQUIRE(duckdb_query(connection, "SELECT dia_r_sum([NULL]::DOUBLE[])", &result) == DuckDBSuccess, "missing query failed");
    REQUIRE(duckdb_value_is_null(&result, 0, 0), "all missing differs");
    checks++;
    duckdb_destroy_result(&result);
    REQUIRE(duckdb_query(connection, "SELECT dia_r_sum([-1]::DOUBLE[])", &result) == DuckDBError, "negative accepted");
    checks++;
    duckdb_destroy_result(&result);
    REQUIRE(duckdb_query(connection, "SELECT dia_r_sum(['Infinity'::DOUBLE])", &result) == DuckDBError, "infinity accepted");
    checks++;
    duckdb_destroy_result(&result);
    REQUIRE(duckdb_query(connection, "SELECT CASE WHEN current_setting('allow_unsigned_extensions') THEN 1 ELSE 0 END", &result) == DuckDBSuccess, "policy query failed");
    REQUIRE(duckdb_value_double(&result, 0, 0) == 0.0, "unsigned loading enabled");
    checks++;
    duckdb_destroy_result(&result);
    REQUIRE(duckdb_prepare(connection, "SELECT range::DOUBLE AS value FROM range(100001)", &prepared) == DuckDBSuccess, "prepare failed");
    REQUIRE(duckdb_execute_prepared_streaming(prepared, &result) == DuckDBSuccess, "stream failed");
    REQUIRE(duckdb_result_is_streaming(result), "query silently materialized its result");
    checks++;
    idx_t rows = 0, chunks = 0;
    double total = 0;
    while ((chunk = duckdb_fetch_chunk(result)) != NULL) {
        idx_t n = duckdb_data_chunk_get_size(chunk);
        REQUIRE(n > 0 && n <= 10000, "result chunk exceeded bound");
        double *values = duckdb_vector_get_data(duckdb_data_chunk_get_vector(chunk, 0));
        for (idx_t i = 0; i < n; i++) total += values[i];
        rows += n;
        chunks++;
        duckdb_destroy_data_chunk(&chunk);
    }
    REQUIRE(!duckdb_result_error(&result), "fetch failed");
    REQUIRE(rows == 100001 && chunks > 1 && total == 5000050000.0, "chunked values differ");
    checks++;
cleanup:
    if (chunk) duckdb_destroy_data_chunk(&chunk);
    if (p_duckdb_destroy_result) duckdb_destroy_result(&result);
    if (prepared) duckdb_destroy_prepare(&prepared);
    if (connection) duckdb_disconnect(&connection);
    if (database) duckdb_close(&database);
    if (config) duckdb_destroy_config(&config);
    if (open_error) duckdb_free(open_error);
    dlclose(handle);
    if (failure_message[0]) Rf_error("public C API feasibility failed: %s", failure_message);
    return Rf_ScalarInteger(checks);
}
