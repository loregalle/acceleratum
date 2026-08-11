# print ----
test_that("print.aclrtm_accelerometery displays header and matrix", {
  x <- make_test_accel(n = 3)
  expect_snapshot(print(x))
})

test_that("print.aclrtm_accelerometery handles missing attributes gracefully", {
  x <- new_accelerometery(matrix(1:6, 3, 2, dimnames = list(NULL, c("x","y"))))
  expect_snapshot(print(x))
})

# head and tail ----
test_that("head preserves sampling_rate and start_time unchanged", {
  x <- make_test_accel(n = 10)
  h <- head(x, 3)
  expect_equal(nrow(h), 3)
  expect_equal(attr(h, "start_time"), attr(x, "start_time"))
  expect_equal(attr(h, "sampling_rate"), attr(x, "sampling_rate"))
})

test_that("tail shifts start_time forward when sampling_rate is known", {
  x <- make_test_accel(n = 10, sampling_rate = 5,
                       start_time = as.POSIXct("2024-01-01", tz = "UTC"))
  t <- tail(x, 3)
  n_dropped <- 10 - 3
  expect_equal(attr(t, "start_time"),
               attr(x, "start_time") + n_dropped / 5)
})

test_that("tail drops start_time (with message) when sampling_rate is unknown", {
  x <- new_accelerometery(matrix(1:20, 10, 2, dimnames = list(NULL, c("x","y"))),
                          start_time = as.POSIXct("2024-01-01", tz = "UTC"))
  expect_message(t <- tail(x, 3), "cannot be")
  expect_null(attr(t, "start_time"))
})

test_that("tail returns unchanged object when nothing is dropped", {
  x <- make_test_accel(n = 5)
  t <- tail(x, 10)  # n >= nrow(x), nothing dropped
  expect_identical(t, x)
})

# square bracket indexing ----
test_that("column-only subsetting (missing i) keeps attributes untouched", {
  x <- make_test_accel(n = 5)
  sub <- x[, "x"]
  # single column with drop = FALSE default -> stays a matrix
  expect_equal(attr(sub, "sampling_rate"), attr(x, "sampling_rate"))
  expect_equal(attr(sub, "start_time"), attr(x, "start_time"))
})

test_that("sequential row subsetting shifts start_time", {
  x <- make_test_accel(n = 10, sampling_rate = 5,
                       start_time = as.POSIXct("2024-01-01", tz = "UTC"))
  sub <- x[4:10, ]
  expect_equal(attr(sub, "start_time"), attr(x, "start_time") + 3 / 5)
})

test_that("non-sequential row subsetting drops start_time and sampling_rate, with message", {
  x <- make_test_accel(n = 10)
  expect_message(sub <- x[c(1, 3, 5), ], "not sequential")
  expect_null(attr(sub, "start_time"))
  expect_null(attr(sub, "sampling_rate"))
})

test_that("empty row selection returns object with no attributes", {
  x <- make_test_accel(n = 10)
  sub <- x[x[, "x"] > 1e6, ]  # matches nothing
  expect_null(attr(sub, "start_time"))
  expect_null(attr(sub, "sampling_rate"))
})

test_that("subsetting to a single row falls back to a plain vector when drop=TRUE", {
  x <- make_test_accel(n = 5)
  sub <- x[1, , drop = TRUE]
  expect_false(inherits(sub, "aclrtm_accelerometery"))
})

# plotting ----
test_that("plot runs without error for default axes", {
  x <- make_test_accel(n = 10)
  expect_no_error(plot(x))
})

test_that("plot errors when none of the requested axes are present", {
  x <- make_test_accel(n = 10, axes = "xy")
  expect_error(plot(x, axes = "z"), "None of the requested axes")
})

test_that("plot subsets rows via y", {
  x <- make_test_accel(n = 10)
  expect_no_error(plot(x, y = 1:5))
})

test_that("plot returns x invisibly", {
  x <- make_test_accel(n = 5)
  # capture without printing to the graphics device inside the test output
  grDevices::pdf(NULL)
  on.exit(dev.off())
  expect_identical(withVisible(plot(x))$visible, FALSE)
})
