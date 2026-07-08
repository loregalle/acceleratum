#' @export
print.aclrtm_accelerometry <- function(x, ...) {
  sr <- attr(x, "sampling_rate")
  st <- attr(x, "start_time")
  axes <- paste(colnames(x), collapse = "")

  cat(sprintf(
    "<aclrtm_accelerometry>  [%d samples \u00d7 %d axes (%s)]",
    nrow(x), ncol(x), axes
  ))

  if (!is.null(sr)) {
    cat(sprintf("  |  sampling_rate: %.4g Hz", sr))
  }
  if (!is.null(st)) {
    cat(sprintf("  |  start_time: %s", format(st)))
  }
  cat("\n")

  # print the underlying matrix, temporarily stripping the extra attributes
  # so the output stays clean
  m <- unclass(x)
  attr(m, "sampling_rate") <- NULL
  attr(m, "start_time")    <- NULL
  print(m, ...)
  invisible(x)
}


#' @export
summary.aclrtm_accelerometry <- function(object, ...) {
  sr <- attr(object, "sampling_rate")
  st <- attr(object, "start_time")
  axes <- paste(colnames(object), collapse = ", ")

  cat("aclrtm_accelerometry object\n")
  cat("  Samples      :", nrow(object), "\n")
  cat("  Axes         :", axes, "\n")
  if (!is.null(sr)) cat(sprintf("  Sampling rate: %.4g Hz\n", sr))
  if (!is.null(st)) cat("  Start time   :", format(st), "\n")
  cat("\nColumn summaries:\n")
  print(apply(object, 2L, summary))
  invisible(object)
}

#' @export
head.aclrtm_accelerometry <- function(x, n = 6L, ...) {
  sr <- attr(x, "sampling_rate")
  st <- attr(x, "start_time")
  new_accelerometry(NextMethod(), sampling_rate = sr, start_time = st)
}

#' @export
tail.aclrtm_accelerometry <- function(x, n = 6L, ...) {
  sr <- attr(x, "sampling_rate")
  st <- attr(x, "start_time")

  result <- NextMethod()

  # if we have both attributes, shift start_time forward by the number of
  # dropped rows
  n_dropped <- nrow(x) - nrow(result)

  if (n_dropped == 0) {
    return(x)
  }

  if (!is.null(st) && !is.null(sr)) {
    st <- st + n_dropped / sr
  } else if (!is.null(st) && is.null(sr)) {
    message(
      "Note: sampling rate is not defined. Start time cannot be ",
      "estimated."
    )
    st <- NULL
  }

  new_accelerometry(result, sampling_rate = sr, start_time = st)
}

#' @export
`[.aclrtm_accelerometry` <- function(x, i, j, ..., drop = FALSE) {

  sr <- attr(x, "sampling_rate")
  st <- attr(x, "start_time")

  # delegate subsetting to the matrix method
  result <- NextMethod()

  # if result is no longer a matrix (e.g. single row with drop = TRUE)
  # just return it
  if (!is.matrix(result)) return(result)

  # if no row subsetting, sampling_rate and start_time remain unchanged
  if (missing(i)) {
    return(new_accelerometry(result,
                             sampling_rate = sr,
                             start_time    = st))
  }

  # normalise i to positive integer indices to
  # handle logical, negative, and positive integer i
  n <- nrow(x)
  idx <- seq_len(n)[i]

  # if empty selection, early return
  if (length(idx) == 0L) {
    return(new_accelerometry(result,
                             sampling_rate = NULL,
                             start_time    = NULL))
  }

  is_sequential <- length(idx) == 1L || all(diff(idx) == 1L)

  if (is_sequential) {
    n_dropped <- idx[1L] - 1L
    if (!is.null(st) && !is.null(sr) && n_dropped > 0L) {
      st <- st + n_dropped / sr
    } else if (!is.null(st) && is.null(sr) && n_dropped > 0L) {
      message(
        "Note: sampling_rate is not defined. start_time cannot be estimated."
      )
      st <- NULL
    }
  } else {
    message(
      "Note: row subsetting is not sequential. ",
      "start_time and sampling_rate have been dropped."
    )
    st <- NULL
    sr <- NULL
  }

  new_accelerometry(result, sampling_rate = sr, start_time = st)
}

#' @export
plot.aclrtm_accelerometry <- function(x, y, axes = "xyz", ...) {

  # row subsetting via y
  if (!missing(y)) {
    x <- x[y, ]
  }

  # resolve axes to plot
  requested  <- .parse_axes(axes)
  available  <- colnames(x)
  to_plot    <- intersect(requested, available)

  if (length(to_plot) == 0L) {
    stop(
      "None of the requested axes (\"", paste(requested, collapse = ""),
      "\") are present in the object (available: \"",
      paste(available, collapse = ""), "\")."
    )
  }

  # build time axis
  sr <- attr(x, "sampling_rate")
  st <- attr(x, "start_time")
  n  <- nrow(x)

  if (!is.null(sr)) {
    time_axis <- seq(0, by = 1 / sr, length.out = n)
    x_label <- "Time"
  } else {
    time_axis <- seq_len(n)
    x_label   <- "Sample"
  }

  # plot
  n_axes <- length(to_plot)

  op <- graphics::par(
    mfrow = c(n_axes, 1L),
    mar   = c(2, 4, 1, 1),
    oma   = c(3, 0, 2, 0)
  )
  on.exit(graphics::par(op))

  for (i in seq_along(to_plot)) {
    axis_name <- to_plot[i]
    graphics::plot(
      time_axis, x[, axis_name],
      type = "l",
      ylab = paste0("acc. ", axis_name, " (g)"),
      xlab = "",
      xaxt = if (i < n_axes) "n" else "s",   # only draw x axis on bottom panel
      ...
    )
    graphics::mtext(axis_name, side = 3, line = 0, adj = 0.01, cex = 0.8)
  }

  # shared x and title labels
  graphics::mtext(x_label, side = 1, outer = TRUE, line = 1.5)
  graphics::mtext("Accelerometry", side = 3, outer = TRUE, line = 0.5)

  invisible(x)
}

#' Convert an aclrtm_accelerometry object to an aclrtm_burst object
#'
#' Splits the accelerometry matrix into a list of burst matrices and wraps
#' the result in an \code{aclrtm_burst} object.  Bursts are defined either
#' by duration (\code{burst_length}, in seconds) or by number of samples
#' (\code{burst_size}).  \code{burst_length} takes precedence when both are
#' supplied.
#'
#' @param x           An \code{aclrtm_accelerometry} object.
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
#'   \code{aclrtm_accelerometry} via \code{\link{burst_to_accelerometry}} and
#'   back again with \code{accelerometry_to_burst} is not a lossless round-trip.
#'   An \code{aclrtm_accelerometry} object does not store individual timestamps
#'   — only \code{start_time} and \code{sampling_rate} — and assumes a fixed,
#'   regular sampling interval across the entire recording.  Consequently, if
#'   the original \code{aclrtm_burst} object had irregular timestamps (e.g.
#'   variable gaps between bursts) or variable burst sizes (i.e. ragged
#'   bursts), that information is not recoverable from the accelerometry
#'   object alone.  The reconstructed burst object will have uniformly spaced
#'   timestamps and equal-sized bursts (except possibly the last), regardless
#'   of the original structure.
#'   All additional columns present in the original object will also be lost.
#'
#' @return A \code{aclrtm_burst} object.
#' @export
accelerometry_to_burst <- function(x,
                                   burst_length = NULL,
                                   burst_size   = NULL,
                                   ts_col       = "timestamp",
                                   data_col     = "burst") {

  if (!inherits(x, "aclrtm_accelerometry")) {
    stop("`x` must be an aclrtm_accelerometry object.", call. = FALSE)
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
    m <- x[i:min(i + burst_size - 1L, n_samples), , drop = FALSE]
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


#' @rdname accelerometry_to_burst
#' @export
a2b <- function(x,
                burst_length = NULL,
                burst_size   = NULL,
                ts_col       = "timestamp",
                data_col     = "burst") {
  accelerometry_to_burst(x,
                         burst_length = burst_length,
                         burst_size   = burst_size,
                         ts_col       = ts_col,
                         data_col     = data_col)
}
