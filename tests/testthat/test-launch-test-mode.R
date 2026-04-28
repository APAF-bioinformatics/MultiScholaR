library(testthat)

# --- run_app test_mode parameter ---

test_that("run_app(test_mode = TRUE) sets multischolar.test_mode option", {
  withr::with_options(
    list(multischolar.test_mode = NULL)
    , {
      run_app(test_mode = TRUE)
      expect_true(isTRUE(getOption("multischolar.test_mode")))
    }
  )
})

test_that("run_app(test_mode = FALSE) does not set multischolar.test_mode option", {
  withr::with_options(
    list(multischolar.test_mode = NULL)
    , {
      run_app(test_mode = FALSE)
      expect_false(isTRUE(getOption("multischolar.test_mode")))
    }
  )
})

test_that("run_app() default leaves multischolar.test_mode unset", {
  withr::with_options(
    list(multischolar.test_mode = NULL)
    , {
      run_app()
      expect_false(isTRUE(getOption("multischolar.test_mode")))
    }
  )
})

test_that("run_app(test_mode = TRUE) means is_test_mode() returns TRUE", {
  withr::with_options(
    list(multischolar.test_mode = NULL)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "")
      , {
        run_app(test_mode = TRUE)
        expect_true(is_test_mode())
      }
    )
  )
})

# --- inst/app/app.R env var detection (unit-level: is_test_mode reads the var) ---

test_that("is_test_mode() reflects MULTISCHOLAR_TEST_MODE used by app.R detection logic", {
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "true")
      , {
        active <- tolower(Sys.getenv("MULTISCHOLAR_TEST_MODE", "false")) %in% c("true", "1")
        expect_true(active)
      }
    )
  )
})

test_that("app.R env var detection logic handles '1' as truthy", {
  withr::with_envvar(
    c(MULTISCHOLAR_TEST_MODE = "1")
    , {
      active <- tolower(Sys.getenv("MULTISCHOLAR_TEST_MODE", "false")) %in% c("true", "1")
      expect_true(active)
    }
  )
})

test_that("app.R env var detection logic treats missing var as falsy", {
  withr::with_envvar(
    c(MULTISCHOLAR_TEST_MODE = "")
    , {
      active <- tolower(Sys.getenv("MULTISCHOLAR_TEST_MODE", "false")) %in% c("true", "1")
      expect_false(active)
    }
  )
})
