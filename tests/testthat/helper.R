make_test_accel <- function(n = 10, axes = "xyz", sampling_rate = 5,
                            start_time = as.POSIXct("2024-01-01", tz = "UTC")) {
  x <- matrix(seq_len(n * nchar(axes)), nrow = n)
  colnames(x) <- strsplit(axes, "")[[1]]
  new_accelerometry(x, sampling_rate = sampling_rate, start_time = start_time)
}

make_burst_df <- function(n = 3, axes = "xyz", nrow_each = 4) {
  cols <- strsplit(axes, "")[[1]]
  bursts <- replicate(
    n,
    matrix(seq_len(nrow_each * length(cols)), ncol = length(cols),
           dimnames = list(NULL, cols)),
    simplify = FALSE
  )
  data.frame(id = seq_len(n), burst = I(bursts))
}

# attach a timestamp column to a burst data.frame built by make_burst_df()
add_ts_col <- function(df, ts) {
  df$ts <- ts
  df
}

make_df <- function(rows = rep(40L, 12), ts = NULL, dur = 5) {
  m <- length(rows)
  if (is.null(ts)) ts <- (seq_len(m) - 1) * dur
  first <- c(0, cumsum(rows)[-m])
  mats <- unname(Map(function(f, r) {
    k <- f + seq_len(r) - 1
    cbind(x = k, y = 2 * k, z = 3 * k)
  }, first, rows))
  df <- data.frame(timestamp = ts)
  df$accelerations_raw <- mats
  df
}

make_burst <- function(...) {
  burst(make_df(...), "accelerations_raw", "xyz", "timestamp")
}

n_rows <- function(b) vapply(b$accelerations_raw, nrow, integer(1L))
