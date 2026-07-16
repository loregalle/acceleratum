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
#' Collapses the burst list-column back to a character vector of
#' separator-delimited strings, then delegates to `writer` for the actual I/O.
#' All other columns are passed through unchanged.
#'
#' @param x        An `aclrtm_burst` object.
#' @param filename Path to the output file.
#' @param writer   A writing function that accepts a data.frame as its first
#'                 argument and a file path as its second (e.g. `write.csv`,
#'                 `write.table`, `readr::write_csv`, `data.table::fwrite`).
#' @param sep      Separator used when collapsing each burst vector back to a
#'                 string.
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
