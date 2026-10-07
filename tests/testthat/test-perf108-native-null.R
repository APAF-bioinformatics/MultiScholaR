test_that("native ordered sums preserve missingness for compact and grouped lists", {
    skip_if(Sys.getenv("MULTISCHOLAR_PERF108_NATIVE_NULL") != "true", "bound native repair check")
    scope <- fastduckdb::open_query_fast()
    on.exit(fastduckdb::close_query_fast(scope))
    symbol <- Sys.getenv("MULTISCHOLAR_PERF108_SUM_SYMBOL", "dia_ordered_sum_v1")
    fastduckdb::register_query_native_fast(scope, "MultiScholaR", symbol)
    collect <- function(sql) fastduckdb::query_data_fast(scope, sql)
    cases <- list("NULL::DOUBLE[]" = NA_real_, "[]::DOUBLE[]" = NA_real_,
        "[NULL]::DOUBLE[]" = NA_real_, "[NULL,NULL]::DOUBLE[]" = NA_real_,
        "list_resize([]::DOUBLE[],65)" = NA_real_, "[0,0]::DOUBLE[]" = 0,
        "[0,NULL,2]::DOUBLE[]" = 2, "[1e16,1,1]::DOUBLE[]" = sum(c(1e16, 1, 1)),
        "list_slice([NULL,1e16,1,1,NULL]::DOUBLE[],2,4)" = sum(c(1e16, 1, 1)),
        "list_filter([NULL,0,1,2]::DOUBLE[], x -> x > 0)" = 3,
        "list_resize([1]::DOUBLE[],20005,1::DOUBLE)" = 20005)
    for (input in names(cases)) {
        value <- collect(paste0("SELECT dia_r_sum(", input, ") AS x"))$x
        expect_identical(value, cases[[input]], info = input)
        expect_identical(collect(paste0("SELECT dia_r_sum(", input, ") IS NULL AS x"))$x,
            is.na(cases[[input]]), info = input)
    }
    # Repeated constant children and dictionary/filter paths across engine chunks.
    for (sql in c(
        "SELECT range::INTEGER AS i, dia_r_sum(list_resize([]::DOUBLE[],65)) AS x FROM range(5001)",
        "SELECT range::INTEGER AS i, dia_r_sum(CASE WHEN range%2=0 THEN [NULL,NULL]::DOUBLE[] ELSE [1,NULL]::DOUBLE[] END) AS x FROM range(5001) WHERE range%2=0 ORDER BY i")) {
        result <- fastduckdb::send_query_fast(scope, sql)
        part <- fastduckdb::fetch_query_fast(result)
        expect_identical(part$x, rep(NA_real_, nrow(part)))
        expect_identical(typeof(part$x), "double")
        fastduckdb::clear_query_fast(result)
    }
    grouped <- collect(paste(
        "SELECT g, dia_r_sum(list(x ORDER BY o)) AS x FROM (VALUES",
        "(1,1,1e16::DOUBLE),(1,2,1::DOUBLE),(1,3,1::DOUBLE),",
        "(2,1,NULL::DOUBLE),(2,2,NULL::DOUBLE),",
        "(3,1,0::DOUBLE),(3,2,NULL::DOUBLE),(3,3,2::DOUBLE)) t(g,o,x)",
        "GROUP BY g ORDER BY g"))
    expect_identical(grouped, data.frame(g = 1:3, x = c(sum(c(1e16, 1, 1)), NA_real_, 2)))
    for (value in c("-1", "'Infinity'::DOUBLE", "'NaN'::DOUBLE")) {
        expect_error(collect(paste0("SELECT dia_r_sum([", value, "]::DOUBLE[])")),
            "finite nonnegative")
        expect_identical(collect("SELECT dia_r_sum([1]::DOUBLE[]) AS x")$x, 1)
    }
    paired <- collect(paste(
        "SELECT p.raw AS raw, p.normalised AS normalised FROM (SELECT",
        "dia_paired_r_sum([struct_pack(raw := NULL::DOUBLE, normalised := 0::DOUBLE),",
        "struct_pack(raw := NULL::DOUBLE, normalised := 2::DOUBLE)]) AS p)"))
    expect_identical(paired, data.frame(raw = NA_real_, normalised = 2))
})
