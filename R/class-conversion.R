#' Convert an aclrtm_accelerometery object to an aclrtm_burst object
#'
#' Splits the accelerometery matrix into a list of burst matrices and wraps
#' the result in an \code{aclrtm_burst} object.  Bursts are defined either
#' by duration (\code{burst_length}, in seconds) or by number of samples
#' (\code{burst_size}).  \code{burst_length} takes precedence when both are
#' supplied.
#'
#' @param x           An \code{aclrtm_accelerometery} object.
#' @param burst_length Numeric scalar. Duration of each burst in seconds.
#'   Takes precedence over \code{burst_size} when both are supplied.
#' @param burst_size  Integer scalar. Number of samples per burst. Used only
#'   when \code{burst_length} is \code{NULL}.
#' @param ts_col      Name of the timestamp column in the output.
#'   Defaults to \code{"timestamp"}.
#' @param data_col    Name of the burst data column in the output.
#'   Defaults to \code{"burst"}.
#'
#' @details
#'   A timestamp column is created according to the following rules:
#'   \itemize{
#'     \item Both \code{start_time} and \code{sampling_rate} present: full
#'       timestamp sequence.
#'     \item Only \code{start_time} present: first burst gets \code{start_time},
#'       remaining timestamps are \code{NA}.
#'     \item Only \code{sampling_rate} present: timestamps start at 0.
#'     \item Neither present: no timestamp column.
#'   }
#'   Note that converting an \code{aclrtm_burst} object to
#'   \code{aclrtm_accelerometery} via \code{\link{burst_to_accelerometery}} and
#'   back again with \code{accelerometery_to_burst} is not a lossless round-trip.
#'   An \code{aclrtm_accelerometery} object does not store individual timestamps
#'   — only \code{start_time} and \code{sampling_rate} — and assumes a fixed,
#'   regular sampling interval across the entire recording.  Consequently, if
#'   the original \code{aclrtm_burst} object had irregular timestamps (e.g.
#'   variable gaps between bursts) or variable burst sizes (i.e. ragged
#'   bursts), that information is not recoverable from the accelerometery
#'   object alone.  The reconstructed burst object will have uniformly spaced
#'   timestamps and equal-sized bursts (except possibly the last), regardless
#'   of the original structure.
#'   All additional columns present in the original object will also be lost.
#'
#' @return A \code{aclrtm_burst} object.
#' @export
accelerometery_to_burst <- function(x,
                                    burst_length = NULL,
                                    burst_size   = NULL,
                                    ts_col       = "timestamp",
                                    data_col     = "burst") {

  if (!inherits(x, "aclrtm_accelerometery")) {
    stop("`x` must be an aclrtm_accelerometery object.", call. = FALSE)
  }
  if (!is.character(ts_col) || length(ts_col) != 1L) {
    stop("`ts_col` must be a single character string.", call. = FALSE)
  }
  if (!is.character(data_col) || length(data_col) != 1L) {
    stop("`data_col` must be a single character string.", call. = FALSE)
  }

  axes          <- paste(colnames(x), collapse = "")
  sampling_rate <- attr(x, "sampling_rate")
  start_time    <- attr(x, "start_time")
  n_samples     <- nrow(x)

  # resolve burst_size from burst_length
  if (!is.null(burst_length)) {
    if (!is.numeric(burst_length) || length(burst_length) != 1L ||
        !is.finite(burst_length) || burst_length <= 0) {
      stop("`burst_length` must be a single positive finite number.", call. = FALSE)
    }
    if (is.null(sampling_rate)) {
      stop(
        paste0("`burst_length` requires `sampling_rate` to be defined in `x`. ",
               "Use `burst_size` instead, or reconstruct `x` with a `sampling_rate`."),
        call. = FALSE
      )
    }
    burst_size <- floor(burst_length * sampling_rate)
    if (burst_size < 1L) {
      stop(
        sprintf(
          paste0("`burst_length` (%.4g s) is shorter than one sample at the ",
                 "current `sampling_rate` (%.4g Hz)."),
          burst_length, sampling_rate
        ),
        call. = FALSE
      )
    }
  } else {
    if (is.null(burst_size)) {
      stop("One of `burst_length` or `burst_size` must be supplied.", call. = FALSE)
    }
    if (!is.numeric(burst_size) || length(burst_size) != 1L ||
        burst_size < 1L) {
      stop("`burst_size` must be a single positive integer.", call. = FALSE)
    }
    burst_size <- .as_integer_strict(burst_size)
  }

  # split matrix into burst matrices
  starts     <- seq(1L, n_samples, by = burst_size)
  burst_list <- lapply(starts, function(i) {
    m <- suppressMessages(x[i:min(i + burst_size - 1L, n_samples), , drop = FALSE])
    attributes(m)[c("sampling_rate", "start_time")] <- NULL
    class(m) <- c("matrix", "array")
    m
  })

  # build timestamp column
  has_sr <- !is.null(sampling_rate)
  has_st <- !is.null(start_time)

  if (has_sr && has_st) {

    burst_times <- start_time + (starts - 1L) / sampling_rate
    if (inherits(start_time, "POSIXct")) {
      burst_times <- as.POSIXct(burst_times,
                                tz = attr(start_time, "tzone") %||% "UTC")
    }

  } else if (has_st && !has_sr) {

    message(
      "`sampling_rate` not defined in `x`: only the first timestamp can be ",
      "set. Remaining timestamps will be NA."
    )
    burst_times        <- rep(NA_real_, length(starts))
    burst_times[[1L]]  <- as.numeric(start_time)
    if (inherits(start_time, "POSIXct")) {
      burst_times <- as.POSIXct(burst_times,
                                tz = attr(start_time, "tzone") %||% "UTC")
    }

  } else if (has_sr && !has_st) {

    message(
      "`start_time` not defined in `x`: timestamps will start at 0."
    )
    burst_times <- (starts - 1L) / sampling_rate

  } else {
    burst_times <- NULL
  }

  # assemble output data.frame
  out_ts_col <- if (!is.null(burst_times)) ts_col else NULL

  out_df <- if (!is.null(burst_times)) {
    data.frame(stats::setNames(list(burst_times), ts_col))
  } else {
    data.frame(row.names = seq_along(burst_list))
  }

  out_df[[data_col]] <- burst_list

  new_burst(out_df,
            data_col = data_col,
            axes     = axes,
            ts_col   = out_ts_col)
}


#' @rdname accelerometery_to_burst
#' @export
a2b <- function(x,
                burst_length = NULL,
                burst_size   = NULL,
                ts_col       = "timestamp",
                data_col     = "burst") {
  accelerometery_to_burst(x,
                          burst_length = burst_length,
                          burst_size   = burst_size,
                          ts_col       = ts_col,
                          data_col     = data_col)
}

#' Converts an aclrtm_burst object to an aclrtm_accelerometery object
#'
#' Stacks all burst matrices into a single matrix and
#' constructs an \code{aclrtm_accelerometery} object.  If a timestamp column
#' is present, \code{start_time} is taken from its first value and
#' \code{sampling_rate} is estimated as the median number of samples per
#' second across bursts.  If no timestamp column is defined, \code{start_time}
#' and \code{sampling_rate} can be supplied by the user.
#'
#' @param x An \code{aclrtm_burst} object.
#' @param sampling_rate A single numeric value (Hz), or \code{NULL} (default).
#' @param start_time    A \code{POSIXct} scalar, a positive numeric, or
#'   \code{NULL} (default).
#' @param ...           Additional arguments forwarded to
#'   \code{accelerometery.matrix()}.
#'
#' @details
#'   Converting to \code{aclrtm_accelerometery} and back via
#'   \code{\link{accelerometery_to_burst}} is not a lossless round-trip.
#'   Individual burst timestamps and variable burst sizes are not retained in
#'   the accelerometery representation, which stores only \code{start_time} and
#'   \code{sampling_rate} and assumes a regular sampling interval throughout.
#'   All additional columns present in the original object will also be lost.
#'
#' @return An \code{aclrtm_accelerometery} object.
#' @export
burst_to_accelerometery <- function(x,
                                    sampling_rate = NULL,
                                    start_time    = NULL,
                                    ...) {

  if (!inherits(x, "aclrtm_burst")) {
    stop("`x` must be an aclrtm_burst object.", call. = FALSE)
  }

  if (!is.null(sampling_rate) &&
      (!is.numeric(sampling_rate) || length(sampling_rate) != 1L ||
       !is.finite(sampling_rate) || sampling_rate <= 0)) {
    stop(
      "`sampling_rate` must be a single positive finite number.",
      call. = FALSE
    )
  }

  if (!is.null(start_time) &&
      (length(start_time) != 1L ||
       (!inherits(start_time, "POSIXct") &&
        (!is.numeric(start_time) || start_time < 0)))) {
    stop(
      "`start_time` must be a single POSIXct or non-negative numeric value.",
      call. = FALSE
    )
  }

  data_col <- attr(x, "data_col")
  ts_col   <- attr(x, "ts_col")
  axes     <- attr(x, "axes")

  burst_list <- x[[data_col]]

  # stack bursts
  mat <- do.call(rbind, burst_list)

  # timestamp-derived metadata
  if (!is.null(ts_col)) {

    ts_vals <- x[[ts_col]]
    ts_num  <- as.numeric(ts_vals)

    # start_time from first timestamp, preserving POSIXct if applicable
    if (inherits(ts_vals, c("POSIXct", "POSIXlt"))) {
      start_time <- as.POSIXct(ts_num[[1L]],
                               tz = attr(ts_vals, "tzone") %||% "UTC")
    } else {
      start_time <- ts_num[[1L]]
    }

    # median sampling rate: for each burst, duration = diff to next timestamp;
    # rate = nrow(burst) / duration
    n_bursts    <- length(burst_list)
    check_n     <- min(5L, n_bursts - 1L)   # need at least two timestamps

    if (check_n >= 1L) {
      ts_diffs    <- diff(ts_num)                          # length n_bursts - 1
      burst_nrows <- vapply(burst_list, nrow, integer(1L)) # length n_bursts

      # pair each burst with the gap to the *next* burst
      implied_rates <- burst_nrows[-n_bursts] / ts_diffs   # length n_bursts - 1

      estimated_rate <- stats::median(implied_rates)

      # consistency check against user-supplied sampling_rate
      if (!is.null(sampling_rate)) {
        check_rate <- stats::median(implied_rates[seq_len(check_n)])
        if (abs(sampling_rate - check_rate) > 0.3) {
          warning(
            sprintf(
              paste0("User-supplied `sampling_rate` (%.4g Hz) differs from the ",
                     "rate implied by the first %d burst(s) (%.4g Hz) by more ",
                     "than 0.3 Hz. Using user-supplied value."),
              sampling_rate, check_n, check_rate
            ),
            call. = FALSE
          )
        }
        estimated_rate <- sampling_rate
      }

      sampling_rate <- estimated_rate

    } else {
      warning(
        "Only one burst present; cannot estimate `sampling_rate` from timestamps.",
        call. = FALSE
      )
      # sampling_rate stays as user-supplied or NULL
    }

  }

  # build accelerometery object
  accelerometery.matrix(mat,
                        axes          = axes,
                        sampling_rate = sampling_rate,
                        start_time    = start_time,
                        ...)
}

#' @rdname burst_to_accelerometery
#' @export
b2a <- function(x, sampling_rate = NULL, start_time = NULL, ...) {
  burst_to_accelerometery(x,
                          sampling_rate = sampling_rate,
                          start_time    = start_time,
                          ...)
}
