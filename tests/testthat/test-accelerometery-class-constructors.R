test_that("numeric method builds a matrix with inferred axes", {
  x <- 1:12
  expect_message(out <- accelerometery(x), "inferred axes = \"xyz\"")
  expect_s3_class(out, "aclrtm_accelerometery")
  expect_equal(colnames(out), c("x", "y", "z"))
  expect_equal(dim(out), c(4, 3))
})

test_that("numeric method respects explicit axes", {
  x <- 1:8
  out <- accelerometery(x, axes = "yz")
  expect_equal(colnames(out), c("y", "z"))
  expect_equal(dim(out), c(4, 2))
})

test_that("numeric method errors when length isn't divisible by number of axes", {
  expect_error(accelerometery(1:10, axes = "xyz"), "not divisible")
})

test_that("numeric method attaches sampling_rate and start_time", {
  out <- accelerometery(1:6, axes = "xy", sampling_rate = 4, start_time = as.POSIXct(0))
  expect_equal(attr(out, "sampling_rate"), 4)
  expect_equal(attr(out, "start_time"), as.POSIXct(0))
})

# accelerometery.matrix ----

test_that("matrix method infers axes when no column names are present", {
  x <- matrix(1:12, 4, 3)
  expect_message(out <- accelerometery(x), "inferred axes = \"xyz\"")
  expect_equal(colnames(out), c("x", "y", "z"))
})

test_that("matrix method keeps valid existing column names", {
  x <- matrix(1:8, 4, 2, dimnames = list(NULL, c("y", "z")))
  out <- accelerometery(x)
  expect_equal(colnames(out), c("y", "z"))
})

test_that("matrix method warns and overrides invalid existing column names", {
  x <- matrix(1:12, 4, 3, dimnames = list(NULL, letters[1:3]))
  expect_warning(out <- suppressMessages(accelerometery(x)),
                 "not valid axis labels")
  expect_equal(colnames(out), c("x", "y", "z"))
})

test_that("matrix method warns when axes overrides existing column names", {
  x <- matrix(1:12, 4, 3, dimnames = list(NULL, c("x", "y", "z")))
  expect_warning(accelerometery(x, axes = "zyx"), "overridden by")
})

test_that("matrix method errors on column count mismatch with axes", {
  x <- matrix(1:8, 4, 2)
  expect_error(accelerometery(x, axes = "xyz"), "but `x` has")
})

test_that("matrix method rejects matrices outside 1-3 columns", {
  x <- matrix(1:16, 4, 4)
  expect_error(accelerometery(x), "1, 2, or 3 columns")
})

# accelerometery.data.frame ----

test_that("data.frame method derives start_time and sampling_rate from ts_col", {
  df <- data.frame(
    x = 1:4, y = 5:8, z = 9:12,
    timestamp = as.POSIXct(seq(.5, 2, .5), tz = "UTC")
  )
  expect_message(out <- accelerometery(df, ts_col = "timestamp"),
                 "estimated as 2 Hz")
  expect_equal(attr(out, "start_time"), as.POSIXct(0.5, tz = "UTC"))
  expect_equal(attr(out, "sampling_rate"), 2) # 1 / 0.5s spacing
})

test_that("data.frame method lets explicit sampling_rate override the estimate", {
  df <- data.frame(
    x = 1:4, timestamp = as.POSIXct(seq(.5, 2, .5), tz = "UTC")
  )
  out <- accelerometery(df, sampling_rate = 99, ts_col = "timestamp")
  expect_equal(attr(out, "sampling_rate"), 99)
})

test_that("data.frame method reorders unsorted rows by ts_col", {
  df <- data.frame(x = c(3, 1, 2), timestamp = c(3, 1, 2))
  expect_message(out <- accelerometery(df, ts_col = "timestamp",
                                       sampling_rate = 1), "rearranged")
  expect_equal(as.numeric(out[, "x"]), c(1, 2, 3))
})

test_that("data.frame method errors with >3 numeric columns and no axes", {
  df <- data.frame(a = 1:4, b = 1:4, c = 1:4, d = 1:4)
  expect_error(accelerometery(df), "at most 3")
})

test_that("data.frame method errors when ts_col is invalid", {
  df <- data.frame(x = 1:4, timestamp = 1:4)
  expect_error(accelerometery(df, ts_col = "nope"), "not found")
  expect_error(accelerometery(df, ts_col = 5), "out of range")
})

# Validation ----
test_that("validate_accelerometery passes through a valid object", {
  x <- new_accelerometery(matrix(1:6, 3, 2, dimnames = list(NULL, c("x","y"))))
  expect_identical(validate_accelerometery(x), x)
})

test_that("validate_accelerometery rejects wrong class", {
  x <- matrix(1:6, 3, 2)
  expect_error(validate_accelerometery(x), "aclrtm_accelerometery.*object")
})

test_that("validate_accelerometery rejects bad column counts", {
  x <- new_accelerometery(matrix(1:12, 3, 4))
  expect_error(validate_accelerometery(x), "1, 2, or 3 columns")
})

test_that("validate_accelerometery rejects invalid column names", {
  x <- new_accelerometery(matrix(1:6, 3, 2, dimnames = list(NULL, c("x","w"))))
  expect_error(validate_accelerometery(x), "subset of c\\('x','y','z'\\)")
})

test_that("validate_accelerometery rejects an invalid sampling_rate", {
  x <- new_accelerometery(matrix(1:6, 3, 2), sampling_rate = -1)
  expect_error(validate_accelerometery(x), "positive finite number")
})

test_that("validate_accelerometery rejects an invalid start_time", {
  x <- new_accelerometery(matrix(1:6, 3, 2), start_time = c(1, 2))
  expect_error(validate_accelerometery(x), "single POSIXct or positive numeric")
})

