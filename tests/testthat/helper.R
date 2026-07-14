make_test_accel <- function(n = 10, axes = "xyz", sampling_rate = 5,
                            start_time = as.POSIXct("2024-01-01", tz = "UTC")) {
  x <- matrix(seq_len(n * nchar(axes)), nrow = n)
  colnames(x) <- strsplit(axes, "")[[1]]
  new_accelerometery(x, sampling_rate = sampling_rate, start_time = start_time)
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
