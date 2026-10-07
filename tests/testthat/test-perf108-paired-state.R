test_that("paired ordered sums preserve independent sequential arithmetic", {
    skip_if(Sys.getenv("MULTISCHOLAR_PERF108_STRUCTURAL") != "A")
    con <- fastduckdb::open_query_fast()
    on.exit(fastduckdb::close_query_fast(con))
    fastduckdb::register_query_native_fast(con, "MultiScholaR", "dia_ordered_sum_v1")
    query <- function(sql) fastduckdb::query_data_fast(con, sql)
    values <- list(c(1e16, 1, 1), c(NA_real_, 0, 1), rep(NA_real_, 3),
        c(.Machine$double.xmin / 2, .Machine$double.xmin / 2, 0),
        c(.Machine$double.xmax, .Machine$double.xmax, NA_real_),
        rep(c(1e16, 1, 1), 10000L))
    for (x in values) {
        y <- rev(x)
        frame <- data.frame(ordinal = seq_along(x), measurement_a = x, measurement_b = y)
        fastduckdb::execute_query_fast(con, "DROP TABLE IF EXISTS fixture")
        fastduckdb::execute_query_fast(con, paste0("CREATE TEMP TABLE fixture ",
            "(ordinal INTEGER, measurement_a DOUBLE, measurement_b DOUBLE)"))
        for (start in seq.int(1L, nrow(frame), by = 10000L)) {
            end <- min(nrow(frame), start + 9999L)
            fastduckdb::append_query_table_fast(con, "fixture", frame[start:end, ])
        }
        actual <- query(paste0("SELECT p.raw AS a, p.normalised AS b, a_ref, b_ref FROM (SELECT ",
            "dia_paired_r_sum(list(struct_pack(raw := measurement_a, normalised := measurement_b) ",
            "ORDER BY ordinal)) p, dia_r_sum(list(measurement_a ORDER BY ordinal)) a_ref, ",
            "dia_r_sum(list(measurement_b ORDER BY ordinal)) b_ref FROM fixture) q"))
        expect_identical(actual$a, actual$a_ref)
        expect_identical(actual$b, actual$b_ref)
        expected <- function(z) if (all(is.na(z))) NA_real_ else sum(z, na.rm = TRUE)
        expect_identical(actual$a, expected(x))
        expect_identical(actual$b, expected(y))
    }
    empty <- query("SELECT p.raw AS a, p.normalised AS b FROM
        (SELECT dia_paired_r_sum(NULL::STRUCT(raw DOUBLE, normalised DOUBLE)[]) p) q")
    expect_true(is.na(empty$a) && is.na(empty$b))
    for (bad in c("-1", "'Infinity'::DOUBLE", "'NaN'::DOUBLE")) {
        expect_error(query(paste0("SELECT p.raw FROM (SELECT dia_paired_r_sum(",
            "[struct_pack(raw := ", bad, ", normalised := 1::DOUBLE)]) p) q")),
            "finite nonnegative")
        expect_identical(query("SELECT 1::DOUBLE AS recovered")$recovered, 1)
    }
    fastduckdb::close_query_fast(con)
    expect_false(fastduckdb::query_is_valid_fast(con))
    expect_error(query("SELECT 1"), "closed")
})
