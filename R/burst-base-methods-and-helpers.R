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

  cat(sprintf("<aclrtm_burst>  [%d rows x %d cols]", nr, ncol(x)))
  cat(sprintf("\n  data_col : \"%s\"", data_col))
  cat(sprintf("\n  axes     : \"%s\"", axes))
  if (!is.null(ts_col)) {
    cat(sprintf("\n  ts_col   : \"%s\"", ts_col))
  } else {
    cat("\n  ts_col   : (none)")
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

#' Plot method for aclrtm_burst objects
#'
#' @param x     An `aclrtm_burst` object.
#' @param y     Optional integer vector of row indices to plot. Must be
#'              *sequential* (e.g. `1:10`, `5:5`), not an arbitrary subset
#'              (e.g. `c(3, 76, 82)` is invalid). If missing, every row/burst
#'              in `x` is plotted.
#' @param axes  A string of axes to plot (e.g. "xyz"). Defaults to the
#'              object's own `axes` attribute.
#' @param ...   Passed on to `graphics::plot()`.
#'
#' @export
plot.aclrtm_burst <- function(x, y, axes = NULL, ...) {

  data_col <- attr(x, "data_col")
  obj_axes <- attr(x, "axes")
  ts_col   <- attr(x, "ts_col")
  n_rows   <- nrow(x)

  # resolve rows to plot
  if (missing(y)) {
    row_idx <- seq_len(n_rows)
  } else {
    row_idx <- as.integer(y)

    if (anyNA(row_idx)) {
      stop("`y` must be an integer/numeric vector of row indices.", call. = FALSE)
    }
    if (any(row_idx < 1L) || any(row_idx > n_rows)) {
      stop(sprintf("`y` contains indices out of range (object has %d row(s)).",
                   n_rows), call. = FALSE)
    }
    if (anyDuplicated(row_idx) != 0L) {
      stop("`y` must not contain duplicate row indices.", call. = FALSE)
    }

    row_idx <- sort(row_idx)
    if (length(row_idx) > 1L && any(diff(row_idx) != 1L)) {
      stop(
        "`y` must select a contiguous, sequential range of rows (e.g. 1:10), ",
        "not an arbitrary subset.",
        call. = FALSE
      )
    }
  }

  n_samples  <- vapply(x[[data_col]][row_idx], nrow, integer(1L))

  # resolve axes to plot
  if (is.null(axes)) axes <- obj_axes
  requested <- .parse_axes(axes)
  available <- colnames(x[[data_col]][[1L]])
  to_plot   <- intersect(requested, available)

  if (length(to_plot) == 0L) {
    stop(
      "None of the requested axes (\"", paste(requested, collapse = ""),
      "\") are present in the burst(s) (available: \"",
      paste(available, collapse = ""), "\")."
    )
  }

  # concatenate the selected bursts into one continuous matrix
  full_mat <- do.call(rbind, lapply(x[[data_col]][row_idx],
                                    function(b) b[, to_plot, drop = FALSE]))

  # build the time axis
  if (!is.null(ts_col)) {
    ts_vals <- x[[ts_col]][row_idx]
    ts_num  <- suppressWarnings(as.numeric(ts_vals))
    ts_num  <- ts_num - min(ts_num, na.rm = TRUE)

    if (anyNA(ts_num)) {
      warning("`ts_col` could not be interpreted numerically; falling back to sample index.",
              call. = FALSE)
      time_axis <- seq_len(nrow(full_mat))
      x_label   <- "Sample"
    } else {
      rates_inv <- diff(ts_num)/n_samples[-length(n_samples)]
      rates_inv <- c(rates_inv, median(rates_inv))
      time_axis <- numeric(0)
      for (i in 1:length(rates_inv)) {
        time_axis <- c(time_axis,
                       seq(ts_num[i], by = rates_inv[i],
                           length.out = n_samples[i]))
      }
      x_label   <- "Time"
    }
  } else {
    time_axis <- seq_len(nrow(full_mat))
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
      time_axis, full_mat[, axis_name],
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

  invisible(x)
}

#' Write an aclrtm_burst object to a delimited file
#'
#' A convenience wrapper for `collapse_burst()` followed by the chosen writer
#' for I/O. Collapses the burst list-column back to a character vector of
#' separator-delimited strings, then delegates to `writer` for the actual I/O.
#' All other columns are passed through unchanged.
#' When working with large objects, for speed purposes, we recommend operating
#' on chunks of `x` in parallel using `collapse_burst()` and call the writer
#' function separately.
#'
#' @param x         An `aclrtm_burst` object.
#' @param filename  Path to the output file.
#' @param writer    A writing function that accepts a data.frame as its first
#'                  argument and a file path as its second (e.g. `write.csv`,
#'                  `write.table`, `readr::write_csv`, `data.table::fwrite`).
#' @param burst_sep Passed to `collapse_burst()`.
#' @param fmt       Passed to `collapse_burst()`.
#' @param ...       Additional arguments forwarded to `writer`.
#'
#' @return Invisibly returns `x` unchanged.
#' @seealso [collapse_burst()]
#' @export
write_burst <- function(x, filename,
                        writer = utils::write.csv,
                        burst_sep = " ",
                        fmt = NULL,
                        ...) {


  if (!is.character(filename) || length(filename) != 1L) {
    stop("`filename` must be a single character string.", call. = FALSE)
  }
  if (!is.function(writer)) {
    stop("`writer` must be a function.", call. = FALSE)
  }

  out <- collapse_burst(x, burst_sep, fmt)

  writer(out, filename, ...)
  invisible(x)
}

#' Collapse burst
#'
#' Collapses the burst list-column back to a character vector of
#' separator-delimited strings.
#' @param x          An `aclrtm_burst` object.
#' @param burst_sep  Separator used when collapsing each burst vector back to a
#'                   string.
#' @param fmt        Optional. Passed to `sprintf()` for values formatting.
#' @returns Returns an object stripped of its `aclrtm_burst` class and with a
#' collapsed burst column.
#' @seealso [write_burst()]
#' @export
collapse_burst <- function(x, burst_sep = " ", fmt = NULL) {
  if (!inherits(x, "aclrtm_burst")) {
    stop("`x` must be an aclrtm_burst object.", call. = FALSE)
  }

  data_col <- attr(x, "data_col")

  # Flatten burst column: each numeric vector → one delimited string
  out <- x
  class(out) <- class(out)[class(out) != "aclrtm_burst"]

  out[[data_col]] <- vapply(
    x[[data_col]],
    function(v) {
      v <- t(v)
      if (!is.null(fmt)) {
        v <- sprintf(fmt, v)
      }
      paste(v, collapse = burst_sep)
    },
    character(1L)
  )

  out
}

#' Re-cut a aclrtm_burst object
#'
#' From a `aclrtm_burst` object, offsets a new object of the same class and/or
#' changes burst length.
#'
#' @param x             A aclrtm_burst object (attribute `ts_col` must be set).
#' @param length        New burst duration in seconds.
#' @param sampling_rate Sampling rate in Hz. If `NULL` (default), it is
#'                      estimated as the median of rows / time-to-next-burst.
#' @param offset        Either a numeric value that offsets, in seconds,
#'                      the start of the first bin. Can be a `POSIXct` value,
#'                      in which case the first bin will be offset to this
#'                      value.
#' @param min_fill      Bins holding fewer than
#'                      `min_fill * length * sampling_rate` rows
#'                      (e.g. at gaps or the end of the data) are dropped.
#' @param gap_tol       A burst whose distance to the next one exceeds
#'                      `gap_tol * median distance` is treated as followed by
#'                      a gap.
#'
#' @details
#' Akin to chaining `b2a()` and `a2b()` with the difference that
#' this function will interpolate timestamps between bursts rather than lose
#' that information entirely. Additionally, columns that have a single value
#' (e.g. a recording ID for the data) will be retained in the output.
#' Columns with more than 1 value cannot be assigned to the output rows
#' unequivocally and will be dropped with a warning.
#'
#' @seealso [burst()], [a2b()], [b2a()]
#'
#' @return A aclrtm_burst object with columns `ts_col` and `data_col` and any
#' additional column with unique value.
#' @export
reburst <- function(x,
                    length,
                    sampling_rate = NULL,
                    offset = 0,
                    min_fill = 0.9,
                    gap_tol = 1.5) {
  stopifnot(inherits(x, "aclrtm_burst"),
            is.numeric(length), length(length) == 1L, length > 0)

  data_col <- attr(x, "data_col")
  ts_col   <- attr(x, "ts_col")
  axes     <- attr(x, "axes")
  if (is.null(ts_col)) stop("`x` has no `ts_col`; cannot reburst.", call. = FALSE)

  # columns other than ts/data: keep those holding a single value, drop the rest
  other_cols <- setdiff(names(x), c(ts_col, data_col))
  is_const   <- vapply(other_cols,
                       function(col) length(unique(x[[col]])) == 1L,
                       logical(1L))
  const_cols <- other_cols[is_const]
  lost_cols  <- other_cols[!is_const]

  if (length(lost_cols) > 0L) {
    warning(sprintf(
      "Column%s %s hold%s more than one value and cannot be assigned to the new bursts; dropped.",
      if (length(lost_cols) == 1L) "" else "s",
      paste0("`", lost_cols, "`", collapse = ", "),
      if (length(lost_cols) == 1L) "s" else ""
    ), call. = FALSE)
  }

  ts_in    <- x[[ts_col]]
  is_posix <- inherits(ts_in, "POSIXct")
  tz       <- if (is_posix) attr(ts_in, "tzone") else NULL
  t        <- as.numeric(ts_in)
  mats     <- x[[data_col]]

  m <- length(t)
  n <- vapply(mats, nrow, integer(1L))

  # distance to next burst (NA for the last one)
  dur <- c(diff(t), NA_real_)

  if (is.null(sampling_rate)) {
    if (m < 2L) {
      stop("Need at least two bursts to estimate `sampling_rate`.",
           call. = FALSE)
    }
    sampling_rate <- median(n[-m] / dur[-m])
  }

  # gaps: only transitions between bursts count (the last burst has no successor)
  typ      <- median(dur, na.rm = TRUE)
  gap      <- is.na(dur) | dur > gap_tol * typ
  n_gaps   <- sum(gap[-m])

  if (n_gaps > 0L) {
    warning(sprintf(
      "%d gap%s identified between bursts (distance to next burst > %.2f x the median of %.3f s).",
      n_gaps, if (n_gaps == 1L) "" else "s", gap_tol, typ
    ), call. = FALSE)
  }

  # a bin can only hold a constant number of rows if length * sampling_rate
  # is a whole number
  rows_per_bin <- length * sampling_rate
  if (abs(rows_per_bin - round(rows_per_bin)) > 1e-8) {
    warning(sprintf(
      "`length` = %g s at `sampling_rate` = %g Hz gives %.3f rows per burst; bursts will not have a regular number of rows and bins below the `min_fill` threshold will be dropped.",
      length, sampling_rate, rows_per_bin
    ), call. = FALSE)
  }

  # row spacing within each burst: interpolate to the next burst,
  # fall back on the nominal rate across gaps / for the last burst
  dt <- ifelse(gap, 1 / sampling_rate, dur / n)

  # interpolated timestamp for every row
  rt <- rep(t, n) + (sequence(n) - 1L) * rep(dt, n)

  # bin offset
  if (inherits(offset, "POSIXct")) {
    if (!is_posix) {
      stop("`offset` is POSIXct but the `ts_col` timestamps are not.",
           call. = FALSE)
    }
    t0 <- as.numeric(offset)
  } else if (is.numeric(offset) && length(offset) == 1L && is.finite(offset)) {
    t0 <- t[1L] + offset
  } else {
    stop("`offset` must be a single number (seconds relative to the first ",
         "burst) or a POSIXct.", call. = FALSE)
  }

  # bin assignment
  bin <- floor((rt - t0) / length + 1e-9)

  big    <- do.call(rbind, mats)
  starts <- c(1L, which(diff(bin) != 0L) + 1L)
  ends   <- c(starts[-1L] - 1L, length(bin))

  keep <- (ends - starts + 1L) >= min_fill * rows_per_bin
  # the last bin ends the data rather than a gap, so a short one is expected
  keep[length(keep)] <- TRUE

  n_dropped <- sum(!keep)
  if (n_dropped > 0L) {
    warning(sprintf(
      "%d of %d bins dropped for holding fewer than %.0f%% of the expected %.2f rows.",
      n_dropped, length(keep), 100 * min_fill, rows_per_bin
    ), call. = FALSE)
  }

  starts <- starts[keep]
  ends   <- ends[keep]
  if (!length(starts)) stop("No complete bins produced.", call. = FALSE)

  new_mats <- Map(function(s, e) big[s:e, , drop = FALSE], starts, ends)
  new_ts   <- t0 + bin[starts] * length
  if (is_posix) new_ts <- .POSIXct(new_ts, tz = tz)

  out <- data.frame(new_ts)
  names(out) <- ts_col
  out[[data_col]] <- new_mats

  # re-attach single-valued columns
  n_out <- nrow(out)
  for (col in const_cols) {
    out[[col]] <- x[[col]][rep(1L, n_out)]
  }

  # restore the input column order (minus any dropped columns)
  out <- out[names(x)[names(x) %in% names(out)]]

  new_burst(out, data_col = data_col, axes = axes, ts_col = ts_col)
}
