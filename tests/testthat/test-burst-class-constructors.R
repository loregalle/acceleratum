# new_burst ----
test_that("new_burst attaches class and attributes without validating contents", {
  df <- make_burst_df()
  out <- new_burst(df, data_col = "burst", axes = "xyz", ts_col = NULL)
  expect_s3_class(out, "aclrtm_burst")
  expect_equal(attr(out, "data_col"), "burst")
  expect_equal(attr(out, "axes"), "xyz")
  expect_null(attr(out, "ts_col"))
})

test_that("new_burst does not validate axes against actual burst contents", {
  df <- make_burst_df(axes = "xyz")
  # deliberately mismatched axes string vs actual 3-column data — should NOT error here
  expect_no_error(new_burst(df, data_col = "burst", axes = "xy"))
})

test_that("new_burst enforces basic argument types via stopifnot", {
  df <- make_burst_df()
  expect_error(new_burst(list(), data_col = "burst", axes = "xyz"))
  expect_error(new_burst(df, data_col = 1, axes = "xyz"))
  expect_error(new_burst(df, data_col = "burst", axes = 1))
  expect_error(new_burst(df, data_col = "burst", axes = "xyz", ts_col = 1))
})


# validate_burst ----
test_that("validate_burst passes through a valid object", {
  x <- new_burst(make_burst_df(), data_col = "burst", axes = "xyz")
  expect_identical(validate_burst(x), x)
})

test_that("validate_burst rejects wrong class", {
  expect_error(validate_burst(data.frame(x = 1)), "aclrtm_burst")
})

test_that("validate_burst rejects malformed data_col/axes/ts_col attributes", {
  df <- make_burst_df()
  x1 <- new_burst(df, data_col = "burst", axes = "xyz")
  attr(x1, "data_col") <- 5
  expect_error(validate_burst(x1), "data_col.*must be a single character")

  x2 <- new_burst(df, data_col = "burst", axes = "xyz")
  attr(x2, "axes") <- NULL
  expect_error(validate_burst(x2), "axes.*must be a single character")

  x3 <- new_burst(df, data_col = "burst", axes = "xyz")
  attr(x3, "ts_col") <- c("a", "b")
  expect_error(validate_burst(x3), "ts_col.*must be NULL or a single")
})

test_that("validate_burst errors when data_col or ts_col are missing from df", {
  df <- make_burst_df()
  x <- new_burst(df, data_col = "nope", axes = "xyz")
  expect_error(validate_burst(x), "not found in `x`")

  x2 <- new_burst(df, data_col = "burst", axes = "xyz", ts_col = "missing_ts")
  expect_error(validate_burst(x2), "not found in `x`")
})

test_that("validate_burst rejects a non-list burst column", {
  df <- data.frame(id = 1:3, burst = 1:3)
  x <- new_burst(df, data_col = "burst", axes = "x")
  expect_error(validate_burst(x), "must be a list")
})

test_that("validate_burst rejects non-matrix elements", {
  df <- make_burst_df()
  df$burst[[2]] <- 1:12  # plain vector, not a matrix
  x <- new_burst(df, data_col = "burst", axes = "xyz")
  expect_error(validate_burst(x), "must be matrices.*row\\(s\\): 2")
})

test_that("validate_burst rejects non-numeric matrix elements", {
  df <- make_burst_df()
  df$burst[[2]] <- matrix(letters[1:12], ncol = 3)
  x <- new_burst(df, data_col = "burst", axes = "xyz")
  expect_error(validate_burst(x), "must be numeric.*row\\(s\\): 2")
})

test_that("validate_burst rejects matrices with wrong ncol for axes", {
  df <- make_burst_df(axes = "xyz")
  df$burst[[3]] <- matrix(1:8, ncol = 2)  # 2 cols, but axes = "xyz" (n_axes = 3)
  x <- new_burst(df, data_col = "burst", axes = "xyz")
  expect_error(validate_burst(x), "must have 3 column\\(s\\).*row\\(s\\): 3")
})

test_that("validate_burst warns on empty bursts", {
  df <- make_burst_df()
  df$burst[[1]] <- matrix(numeric(0), ncol = 3)
  x <- new_burst(df, data_col = "burst", axes = "xyz")
  expect_warning(validate_burst(x), "Empty burst.*row\\(s\\): 1")
})

test_that("validate_burst warns on ragged (unequal nrow) bursts", {
  df <- make_burst_df(n = 3, axes = "xyz", nrow_each = 4)
  df$burst[[2]] <- matrix(1:9, ncol = 3)  # 3 rows instead of 4
  x <- new_burst(df, data_col = "burst", axes = "xyz")
  expect_warning(validate_burst(x), "Ragged bursts.*row\\(s\\) 2")
})

test_that("validate_burst warns on duplicate timestamps", {
  df <- make_burst_df()
  df$ts <- c(1, 1, 2)
  x <- new_burst(df, data_col = "burst", axes = "xyz", ts_col = "ts")
  expect_warning(validate_burst(x), "duplicate values")
})

# burst ----
test_that("burst validates df, data_col, and ts_col arguments", {
  df <- make_burst_df()
  expect_error(burst(list(), "burst"), "must be a data.frame")
  expect_error(burst(df, 5), "must be a single character string")
  expect_error(burst(df, "nope"), "not found in `df`")
  expect_error(burst(df, "burst", ts_col = "nope"), "not found in `df`")
})

## character column parsing (bench) ----

test_that("burst parses real space-separated burst data from `bench`", {
  data(bench, package = "acceleratum")
  out <- burst(bench, data_col = "accelerations_raw", ts_col = "timestamp")
  expect_s3_class(out, "aclrtm_burst")
  expect_true(is.list(out[["accelerations_raw"]]))
  expect_true(all(vapply(out[["accelerations_raw"]], is.matrix, logical(1))))
  # every parsed burst should have exactly 3 columns (xyz), varying nrow
  ncols <- vapply(out[["accelerations_raw"]], ncol, integer(1))
  expect_true(all(ncols == 3))
})

test_that("burst infers axes as xyz from bench's raw burst lengths", {
  data(bench, package = "acceleratum")
  expect_message(
    out <- burst(bench, data_col = "accelerations_raw"),
    "inferred axes = \"xyz\""
  )
  expect_equal(attr(out, "axes"), "xyz")
})

## list-of-vectors parsing ----

test_that("burst parses a list column of plain numeric vectors", {
  df <- data.frame(id = 1:2)
  df$burst <- list(1:9, 10:18)  # length 9, divisible by 3
  out <- burst(df, "burst", axes = "xyz")
  expect_equal(dim(out$burst[[1]]), c(3, 3))
})

test_that("burst leaves list-of-matrices columns as-is, only fixing colnames", {
  df <- make_burst_df(axes = "xyz")
  out <- burst(df, "burst", axes = "xyz")
  expect_equal(colnames(out$burst[[1]]), c("x", "y", "z"))
})

## errors ----

test_that("burst errors when a vector length isn't divisible by the number of axes", {
  df <- data.frame(id = 1)
  df$burst <- list(1:10)
  expect_error(burst(df, "burst", axes = "xyz"), "not divisible")
})

test_that("burst propagates validate_burst errors (e.g. ragged input)", {
  df <- make_burst_df(n = 2, axes = "xyz", nrow_each = 4)
  df$burst[[2]] <- matrix(1:3, ncol = 3)  # ragged (1 row instead of 4)
  # burst() calls validate_burst() internally, which only warns for ragged —
  # confirm it surfaces, not that we re-derive every validate_burst case here
  expect_warning(burst(df, "burst", axes = "xyz"), "Ragged bursts")
})
