# regular data ----

test_that("regular 8 Hz data: 5 s bursts -> 3 s bursts", {
  b <- make_burst()                       # 12 bursts x 40 rows = 60 s
  expect_no_warning(out <- reburst(b, 3))

  expect_s3_class(out, "aclrtm_burst")
  expect_identical(attr(out, "data_col"), "accelerations_raw")
  expect_identical(attr(out, "axes"), "xyz")
  expect_identical(attr(out, "ts_col"), "timestamp")

  expect_equal(nrow(out), 20L)
  expect_equal(out$timestamp, seq(0, 57, by = 3))
  expect_true(all(n_rows(out) == 24L))
  expect_true(all(vapply(out$accelerations_raw, ncol, integer(1L)) == 3L))
  expect_equal(colnames(out$accelerations_raw[[1]]), c("x", "y", "z"))

  # first row of bin j is global row 24 * j
  first_x <- vapply(out$accelerations_raw, function(m) m[1, "x"], numeric(1L))
  expect_equal(first_x, 24 * (0:19))

  # no row lost, duplicated or reordered
  expect_equal(do.call(rbind, out$accelerations_raw),
               do.call(rbind, b$accelerations_raw))
})

test_that("5 s bursts -> 10 s bursts", {
  b <- make_burst()
  expect_no_warning(out <- reburst(b, 10))
  expect_equal(out$timestamp, seq(0, 50, by = 10))
  expect_true(all(n_rows(out) == 80L))
})

test_that("reburst to the same length is an identity", {
  b <- make_burst()
  expect_no_warning(out <- reburst(b, 5))
  expect_equal(out$timestamp, b$timestamp)
  expect_equal(out$accelerations_raw, b$accelerations_raw)
})

test_that("supplying the true sampling_rate gives the same result as estimating it", {
  b <- make_burst()
  expect_equal(reburst(b, 3, sampling_rate = 8), reburst(b, 3))
})

test_that("POSIXct timestamps are preserved in class and time zone", {
  start <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  b <- make_burst(ts = start + (0:11) * 5)
  out <- reburst(b, 3)

  expect_s3_class(out$timestamp, "POSIXct")
  expect_identical(attr(out$timestamp, "tzone"), "UTC")
  expect_equal(out$timestamp[1], start)
  expect_equal(as.numeric(diff(out$timestamp), units = "secs"), rep(3, 19))
})

test_that("`offset` shifts the bins and drops the incomplete leading bin", {
  b <- make_burst()
  expect_warning(out <- reburst(b, 3, offset = 1), "1 of 21 bins dropped")

  expect_equal(out$timestamp[1], 1)
  expect_equal(unname(out$accelerations_raw[[1]][1, "x"]), 8)  # row at t = 1 s

  nr <- n_rows(out)
  expect_true(all(head(nr, -1L) == 24L))
  expect_equal(tail(nr, 1L), 16L)   # last bin kept although incomplete
})

test_that("numeric `offset` is an offset in seconds from the first burst", {
  start <- 1740149727
  b <- make_burst(ts = start + (0:11) * 5)

  out <- suppressWarnings(reburst(b, 5, offset = 1))
  expect_equal(out$timestamp[1], start + 1 - 5)   # one partial bin before it
  out <- suppressWarnings(reburst(b, 5, offset = 1, min_fill = 0))
  expect_true(all(out$timestamp %% 5 == (start + 1) %% 5))

  # offset of 0 equals the default
  expect_equal(reburst(b, 3, offset = 0), reburst(b, 3))
})

test_that("POSIXct `offset` is absolute", {
  start <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  b <- make_burst(ts = start + (0:11) * 5)

  out <- suppressWarnings(reburst(b, 3, offset = start + 1))
  expect_true(all(as.numeric(out$timestamp) %% 3 == (as.numeric(start) + 1) %% 3))

  # POSIXct offset with non-POSIXct timestamps is an error
  expect_error(reburst(make_burst(), 3, offset = start), "not")
  expect_error(reburst(make_burst(), 3, offset = "a"), "must be a single number")
  expect_error(reburst(make_burst(), 3, offset = NULL), "must be a single number")
})

# imperfect data ----

test_that("bursts of varying length are stitched on a continuous time axis", {
  rows <- c(40L, 39L, 41L, 40L, 40L, 38L, 42L, rep(40L, 5))   # still 480 rows
  stopifnot(sum(rows) == 480L)
  ts <- cumsum(c(0, rows[-length(rows)])) / 8   # timestamps consistent with 8 Hz
  b  <- suppressWarnings(make_burst(rows = rows, ts = ts))

  expect_no_warning(out <- reburst(b, 3))
  expect_equal(nrow(out), 20L)
  expect_true(all(n_rows(out) == 24L))
  expect_equal(out$timestamp, seq(0, 57, by = 3))

  # every original row is present exactly once, in order
  expect_equal(unname(do.call(rbind, out$accelerations_raw)[, "x"]), as.numeric(0:479))
})

test_that("jittery timestamps conserve rows and give bins within +/- 1 row", {
  jitter <- c(0, .02, -.01, .03, 0, -.02, .01, 0, .02, -.03, .01, 0)
  b <- make_burst(ts = (0:11) * 5 + jitter)

  expect_no_warning(out <- reburst(b, 3, sampling_rate = 8))
  nr <- n_rows(out)
  expect_equal(nrow(out), 20L)
  expect_true(all(abs(nr - 24L) <= 1L))
  expect_equal(sum(nr), 480L)
})

# gaps and warnings ----

test_that("gaps are reported with their count and never interpolated across", {
  ts <- (0:11) * 5
  ts[7:12] <- ts[7:12] + 10                     # two bursts missing after burst 6
  b <- make_burst(ts = ts)

  w <- testthat::capture_warnings(out <- reburst(b, 3))
  expect_match(w, "1 gap identified", all = FALSE)
  expect_match(w, "bins dropped", all = FALSE)  # partial bins at the gap edges

  expect_true(all(n_rows(out[-nrow(out),]) == 24L))
  # nothing was produced inside the gap (burst 6 ends ~30 s, burst 7 starts at 40 s)
  expect_true(all(out$timestamp <= 27 | out$timestamp >= 42))
  # the bin just before the gap holds consecutive rows only
  before <- out$accelerations_raw[[which(out$timestamp == 27)]][, "x"]
  expect_true(all(diff(before) == 1))
})

test_that("a length that does not divide into whole rows warns", {
  b <- make_burst()
  w <- testthat::capture_warnings(out <- reburst(b, 1 / 3, sampling_rate = 8))
  expect_match(w, "2.667 rows per burst", all = FALSE)
  expect_match(w, "59 of 180 bins dropped", all = FALSE)  # 60 two-row bins, last kept

  nr <- n_rows(out)
  expect_true(all(head(nr, -1L) == 3L))
  expect_equal(tail(nr, 1L), 2L)

  out0 <- suppressWarnings(reburst(b, 1 / 3, sampling_rate = 8, min_fill = 0))
  expect_true(all(n_rows(out0) %in% 2:3))
})

test_that("a single burst needs an explicit sampling_rate", {
  b1 <- make_burst(rows = 40L)
  expect_error(reburst(b1, 3), "at least two bursts")

  out <- reburst(b1, 3, sampling_rate = 8)
  expect_equal(nrow(out), 2L)
  expect_equal(n_rows(out), c(24L, 16L))
})

# input validation ----

test_that("invalid input is rejected", {
  b <- make_burst()
  expect_error(reburst(data.frame(a = 1), 3))
  expect_error(reburst(b, -1))
  expect_error(reburst(b, c(1, 2)))
  expect_error(reburst(b, 3, sampling_rate = c(8, 16)))

  no_ts <- burst(make_df(), "accelerations_raw", "xyz")
  expect_error(reburst(no_ts, 3), "no `ts_col`")
})

# real data ----

test_that("reburst works on the bench data (8 Hz, 40 x 3 matrices every 5 s)", {
  skip_if_not(exists("bench"), "`bench` not available")
  b <- burst(bench, "accelerations_raw", "xyz", "timestamp")

  w <- testthat::capture_warnings({
    out <- reburst(b, 3)
    out <- out[-nrow(out),]
  })
  expect_false(any(grepl("gap", w)))            # regular data: no gaps
  expect_false(any(grepl("rows per burst", w))) # 3 s * 8 Hz is a whole number

  expect_true(all(n_rows(out) == 24L))
  expect_equal(unique(diff(as.numeric(out$timestamp))), 3)
})
