# accelerometery_to_burst ----

test_that("errors on non-accelerometery input", {
  expect_error(accelerometery_to_burst(matrix(1:9, 3, 3)),
               "must be an aclrtm_accelerometery")
})

test_that("errors when neither burst_length nor burst_size supplied", {
  x <- make_test_accel(n = 10)
  expect_error(accelerometery_to_burst(x), "must be supplied")
})

test_that("errors on invalid ts_col / data_col arguments", {
  x <- make_test_accel(n = 10)
  expect_error(accelerometery_to_burst(x, burst_size = 4, ts_col = 1),
               "ts_col.*must be a single character")
  expect_error(accelerometery_to_burst(x, burst_size = 4, data_col = c("a", "b")),
               "data_col.*must be a single character")
})

test_that("burst_size correctly determines number and size of bursts", {
  x <- make_test_accel(n = 10, sampling_rate = NULL, start_time = NULL)
  out <- accelerometery_to_burst(x, burst_size = 4)
  expect_s3_class(out, "aclrtm_burst")
  # 10 rows / burst_size 4 -> bursts of 4, 4, 2
  expect_equal(nrow(out), 3)
  expect_equal(vapply(out$burst, nrow, integer(1)), c(4, 4, 2))
})

test_that("burst_length converts to burst_size via sampling_rate", {
  x <- make_test_accel(n = 20, sampling_rate = 5)  # 5 Hz
  out <- accelerometery_to_burst(x, burst_length = 1)  # 1s * 5Hz = 5 samples/burst
  expect_equal(nrow(out), 4)                     # 20 / 5 = 4 bursts
  expect_equal(nrow(out$burst[[1]]), 5)
})

test_that("burst_length takes precedence over burst_size when both supplied", {
  x <- make_test_accel(n = 20, sampling_rate = 5)
  out <- accelerometery_to_burst(x, burst_length = 1, burst_size = 999)
  expect_equal(nrow(out$burst[[1]]), 5)  # derived from burst_length, not burst_size
})

test_that("burst_length requires sampling_rate on x", {
  x <- make_test_accel(n = 10, sampling_rate = NULL)
  expect_error(accelerometery_to_burst(x, burst_length = 1),
               "requires `sampling_rate`")
})

test_that("burst_length shorter than one sample errors", {
  x <- make_test_accel(n = 10, sampling_rate = 1)  # 1 Hz
  expect_error(accelerometery_to_burst(x, burst_length = 0.1),
               "shorter than one sample")
})

test_that("burst_size must be a single positive integer", {
  x <- make_test_accel(n = 10)
  expect_error(accelerometery_to_burst(x, burst_size = 0), "positive integer")
  expect_error(accelerometery_to_burst(x, burst_size = c(1, 2)), "positive integer")
  expect_error(accelerometery_to_burst(x, burst_size = 3.4),
               "not an integer value")
})

test_that("burst matrices carry the correct axes and no leftover attributes", {
  x <- make_test_accel(n = 10, axes = "xyz")
  out <- accelerometery_to_burst(x, burst_size = 4)
  expect_equal(attr(out, "axes"), "xyz")
  first_burst <- out$burst[[1]]
  expect_equal(colnames(first_burst), c("x", "y", "z"))
  expect_null(attr(first_burst, "sampling_rate"))
  expect_null(attr(first_burst, "start_time"))
})

test_that("ts_col and data_col names are respected", {
  x <- make_test_accel(n = 10, sampling_rate = 5)
  out <- accelerometery_to_burst(x, burst_size = 4, ts_col = "t", data_col = "acc")
  expect_true(all(c("t", "acc") %in% names(out)))
  expect_equal(attr(out, "data_col"), "acc")
  expect_equal(attr(out, "ts_col"), "t")
})

test_that("timestamp column: full sequence when both start_time and sampling_rate present", {
  st <- as.POSIXct("2024-01-01", tz = "UTC")
  x <- make_test_accel(n = 10, sampling_rate = 5, start_time = st)
  out <- accelerometery_to_burst(x, burst_size = 5)
  # burst 1 starts at sample 1 (t=0s), burst 2 starts at sample 6 (t=1s)
  expect_equal(out$timestamp, st + c(0, 1))
})

test_that("timestamp column: only first entry set when start_time present but not sampling_rate", {
  st <- as.POSIXct("2024-01-01", tz = "UTC")
  x <- make_test_accel(n = 10, sampling_rate = NULL, start_time = st)
  expect_message(out <- accelerometery_to_burst(x, burst_size = 5),
                 "Remaining timestamps will be NA")
  expect_equal(out$timestamp[1], st)
  expect_true(is.na(out$timestamp[2]))
})

test_that("timestamp column: starts at 0 when sampling_rate present but not start_time", {
  x <- make_test_accel(n = 10, sampling_rate = 5, start_time = NULL)
  expect_message(out <- accelerometery_to_burst(x, burst_size = 5),
                 "will start at 0")
  expect_equal(out$timestamp, c(0, 1))
})

test_that("no timestamp column when neither start_time nor sampling_rate present", {
  x <- make_test_accel(n = 10, sampling_rate = NULL, start_time = NULL)
  out <- accelerometery_to_burst(x, burst_size = 5)
  expect_false("timestamp" %in% names(out))
})

## a2b alias ----

test_that("a2b is equivalent to accelerometery_to_burst", {
  x <- make_test_accel(n = 10, sampling_rate = 5)
  expect_identical(
    a2b(x, burst_size = 4),
    accelerometery_to_burst(x, burst_size = 4)
  )
})

# burst_to_accelerometery ----

## input validation ----

test_that("errors on non-burst input", {
  expect_error(burst_to_accelerometery(data.frame(x = 1)),
               "must be an aclrtm_burst")
})

test_that("errors on invalid sampling_rate", {
  x <- new_burst(make_burst_df(n = 2), data_col = "burst", axes = "xyz")
  expect_error(burst_to_accelerometery(x, sampling_rate = -1),
               "positive finite number")
  expect_error(burst_to_accelerometery(x, sampling_rate = c(1, 2)),
               "positive finite number")
})

test_that("errors on invalid start_time", {
  x <- new_burst(make_burst_df(n = 2), data_col = "burst", axes = "xyz")
  expect_error(burst_to_accelerometery(x, start_time = c(1, 2)),
               "POSIXct or non-negative numeric")
  expect_error(burst_to_accelerometery(x, start_time = -5),
               "POSIXct or non-negative numeric")
})

## stacking behavior ----

test_that("bursts are stacked into a single matrix with correct axes", {
  df <- make_burst_df(n = 3, axes = "xyz", nrow_each = 4)
  x <- new_burst(df, data_col = "burst", axes = "xyz")
  out <- burst_to_accelerometery(x)
  expect_s3_class(out, "aclrtm_accelerometery")
  expect_equal(nrow(out), 12)  # 3 bursts * 4 rows
  expect_equal(colnames(out), c("x", "y", "z"))
})

## no ts_col: user-supplied values pass straight through ----

test_that("without ts_col, sampling_rate and start_time come only from the user", {
  df <- make_burst_df(n = 2, nrow_each = 4)
  x <- new_burst(df, data_col = "burst", axes = "xyz")

  out <- burst_to_accelerometery(x, sampling_rate = 10,
                                 start_time = as.POSIXct("2024-01-01", tz = "UTC"))
  expect_equal(attr(out, "sampling_rate"), 10)
  expect_equal(attr(out, "start_time"), as.POSIXct("2024-01-01", tz = "UTC"))
})

test_that("without ts_col and no user-supplied values, attributes stay NULL", {
  df <- make_burst_df(n = 2, nrow_each = 4)
  x <- new_burst(df, data_col = "burst", axes = "xyz")
  out <- burst_to_accelerometery(x)
  expect_null(attr(out, "sampling_rate"))
  expect_null(attr(out, "start_time"))
})

## ts_col present: start_time derivation ----

test_that("start_time is taken from the first timestamp, preserving POSIXct", {
  ts <- as.POSIXct(c(0, 1, 2), origin = "1970-01-01", tz = "UTC")
  df <- add_ts_col(make_burst_df(n = 3, nrow_each = 5), ts)
  x <- new_burst(df, data_col = "burst", axes = "xyz", ts_col = "ts")
  out <- burst_to_accelerometery(x)
  expect_equal(attr(out, "start_time"), ts[[1]])
  expect_s3_class(attr(out, "start_time"), "POSIXct")
})

test_that("start_time stays numeric when timestamps are numeric", {
  df <- add_ts_col(make_burst_df(n = 3, nrow_each = 5), c(0, 1, 2))
  x <- new_burst(df, data_col = "burst", axes = "xyz", ts_col = "ts")
  out <- burst_to_accelerometery(x)
  expect_equal(attr(out, "start_time"), 0)
  expect_false(inherits(attr(out, "start_time"), "POSIXct"))
})

## ts_col present: sampling_rate estimation ----

test_that("sampling_rate is estimated as the median implied rate across bursts", {
  # 3 bursts of 5 rows, 1 second apart -> implied rate 5 Hz each -> median 5
  df <- add_ts_col(make_burst_df(n = 3, nrow_each = 5), c(0, 1, 2))
  x <- new_burst(df, data_col = "burst", axes = "xyz", ts_col = "ts")
  out <- burst_to_accelerometery(x)
  expect_equal(attr(out, "sampling_rate"), 5)
})

test_that("sampling_rate estimation handles unequal burst gaps via median", {
  # nrows 4,4,4 over gaps 1,2 -> implied rates 4/1=4, 4/2=2 -> median = 3
  df <- make_burst_df(n = 3, nrow_each = 4)
  df <- add_ts_col(df, c(0, 1, 3))
  x <- new_burst(df, data_col = "burst", axes = "xyz", ts_col = "ts")
  out <- burst_to_accelerometery(x)
  expect_equal(attr(out, "sampling_rate"), 3)
})

test_that("user-supplied sampling_rate close to the estimate is used without warning", {
  df <- add_ts_col(make_burst_df(n = 3, nrow_each = 5), c(0, 1, 2))  # implied rate 5
  x <- new_burst(df, data_col = "burst", axes = "xyz", ts_col = "ts")
  expect_no_warning(out <- burst_to_accelerometery(x, sampling_rate = 5.1))
  expect_equal(attr(out, "sampling_rate"), 5.1)
})

test_that("user-supplied sampling_rate far from the estimate warns but is still used", {
  df <- add_ts_col(make_burst_df(n = 3, nrow_each = 5), c(0, 1, 2))  # implied rate 5
  x <- new_burst(df, data_col = "burst", axes = "xyz", ts_col = "ts")
  expect_warning(out <- burst_to_accelerometery(x, sampling_rate = 10),
                 "differs from the rate implied")
  expect_equal(attr(out, "sampling_rate"), 10)  # user value wins despite the warning
})

test_that("only the first check_n (max 5) bursts are used for the consistency check", {
  # 8 bursts, all gaps = 1s except a deliberately different last gap;
  # check_n = min(5, 7) = 5, so the odd final gap shouldn't affect the warning
  df <- make_burst_df(n = 8, nrow_each = 5)
  ts <- c(0:6, 100)  # last gap is huge (94s) -> would tank the estimate if included
  df <- add_ts_col(df, ts)
  x <- new_burst(df, data_col = "burst", axes = "xyz", ts_col = "ts")
  # implied rate for first 5 gaps is a clean 5 Hz; user supplies 5 -> should NOT warn
  expect_no_warning(burst_to_accelerometery(x, sampling_rate = 5))
})

## ts_col present: single-burst edge case ----

test_that("a single burst with ts_col warns and cannot estimate sampling_rate", {
  df <- add_ts_col(make_burst_df(n = 1, nrow_each = 5), as.POSIXct("2024-01-01", tz = "UTC"))
  x <- new_burst(df, data_col = "burst", axes = "xyz", ts_col = "ts")
  expect_warning(out <- burst_to_accelerometery(x), "Only one burst present")
  expect_null(attr(out, "sampling_rate"))
})

test_that("a single burst with ts_col keeps a user-supplied sampling_rate", {
  df <- add_ts_col(make_burst_df(n = 1, nrow_each = 5), as.POSIXct("2024-01-01", tz = "UTC"))
  x <- new_burst(df, data_col = "burst", axes = "xyz", ts_col = "ts")
  expect_warning(out <- burst_to_accelerometery(x, sampling_rate = 20),
                 "Only one burst present")
  expect_equal(attr(out, "sampling_rate"), 20)
})

## b2a alias ----

test_that("b2a is equivalent to burst_to_accelerometery", {
  df <- add_ts_col(make_burst_df(n = 3, nrow_each = 5), c(0, 1, 2))
  x <- new_burst(df, data_col = "burst", axes = "xyz", ts_col = "ts")
  expect_identical(b2a(x), burst_to_accelerometery(x))
})
