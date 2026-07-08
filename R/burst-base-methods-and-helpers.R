#' @export
`[.aclrtm_burst` <- function(x, i, j, drop = FALSE) {
  data_col <- attr(x, "data_col")
  ts_col   <- attr(x, "ts_col")
  axes     <- attr(x, "axes")

  # Delegate to data.frame subsetting
  result <- NextMethod()

  # Re-wrap as burstdata only if the result is still a data.frame and the
  # burst column survived
  if (is.data.frame(result) && data_col %in% names(result)) {
    new_ts <- if (!is.null(ts_col) && ts_col %in% names(result)) ts_col else NULL
    result <- new_burst(result,
                        data_col = data_col,
                        axes = axes,
                        ts_col = new_ts)
  }

  result
}

#' @export
print.aclrtm_burst <- function(x, n = 10L, ...) {
  data_col <- attr(x, "data_col")
  ts_col   <- attr(x, "ts_col")
  axes     <- attr(x, "axes")
  nr       <- nrow(x)

  cat(sprintf("<aclrtm_burst>  [%d rows x %d cols]\n", nr, ncol(x)))
  cat(sprintf("  data_col : \"%s\"\n", data_col))
  cat(sprintf("  axes     : \"%s\"\n", axes))
  if (!is.null(ts_col)) {
    cat(sprintf("  ts_col   : \"%s\"\n", ts_col))
  } else {
    cat("  ts_col   : (none)\n")
  }
  cat("\n")

  row_idx <- if (n >= 0L) seq_len(min(n, nr)) else seq_len(max(0L, nr + n))

  if (length(row_idx) == 0L) {
    cat("<0 rows>\n")
    return(invisible(x))
  }

  show_df <- as.data.frame(x[row_idx, , drop = FALSE])
  class(show_df) <- "data.frame"

  show_df[[data_col]] <- vapply(
    show_df[[data_col]],
    function(m) {
      if (is.matrix(m) && nrow(m) == 0L) return("<empty matrix>")
      sprintf("<matrix %dx%d>", nrow(m), ncol(m))
    },
    character(1L)
  )

  print(show_df, ...)

  if (nr > length(row_idx)) {
    cat(sprintf("... with %d more rows\n", nr - length(row_idx)))
  }

  invisible(x)
}

#' Write an aclrtm_burst object to a delimited file
#'
#' Collapses the burst list-column back to a character vector of
#' separator-delimited strings, then delegates to `writer` for the actual I/O.
#' All other columns are passed through unchanged.
#'
#' @param x        An `aclrtm_burst` object.
#' @param filename Path to the output file.
#' @param writer   A writing function that accepts a data.frame as its first
#'                 argument and a file path as its second (e.g. `write.csv`,
#'                 `write.table`, `readr::write_csv`, `data.table::fwrite`).
#'                 Defaults to `write.csv`.
#' @param sep      Separator used when collapsing each burst vector back to a
#'                 string. Should match the `sep` used when the object was
#'                 created so the file round-trips cleanly. Defaults to `","`.
#' @param ...      Additional arguments forwarded to `writer`.
#'
#' @return Invisibly returns `x`.
#' @export
write_burst <- function(x, filename,
                        writer = utils::write.csv,
                        sep = " ", ...) {

  if (!inherits(x, "aclrtm_burst")) {
    stop("`x` must be an aclrtm_burst object.", call. = FALSE)
  }
  if (!is.character(filename) || length(filename) != 1L) {
    stop("`filename` must be a single character string.", call. = FALSE)
  }
  if (!is.function(writer)) {
    stop("`writer` must be a function.", call. = FALSE)
  }

  data_col <- attr(x, "data_col")

  # Flatten burst column: each numeric vector → one delimited string
  out <- as.data.frame(x)
  class(out) <- "data.frame"
  out[[data_col]] <- vapply(
    x[[data_col]],
    function(v) paste(as.vector(t(v)), collapse = sep),
    character(1L)
  )

  writer(out, filename, ...)
  invisible(x)
}

#' Converts an aclrtm_burst object to an aclrtm_accelerometry object
#'
#' Stacks all burst matrices into a single matrix and
#' constructs an \code{aclrtm_accelerometry} object.  If a timestamp column
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
#'   \code{accelerometry.matrix()}.
#'
#' @details
#'   Converting to \code{aclrtm_accelerometry} and back via
#'   \code{\link{accelerometry_to_burst}} is not a lossless round-trip.
#'   Individual burst timestamps and variable burst sizes are not retained in
#'   the accelerometry representation, which stores only \code{start_time} and
#'   \code{sampling_rate} and assumes a regular sampling interval throughout.
#'   All additional columns present in the original object will also be lost.
#'
#' @return An \code{aclrtm_accelerometry} object.
#' @export
burst_to_accelerometry <- function(x,
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

  # build accelerometry object
  accelerometry.matrix(mat,
                       axes          = axes,
                       sampling_rate = sampling_rate,
                       start_time    = start_time,
                       ...)
}

#' @rdname burst_to_accelerometry
#' @export
b2a <- function(x, sampling_rate = NULL, start_time = NULL, ...) {
  burst_to_accelerometry(x,
                         sampling_rate = sampling_rate,
                         start_time    = start_time,
                         ...)
}
