#' @export
print.aclrtm_accelerometery <- function(x, ...) {
  sr <- attr(x, "sampling_rate")
  st <- attr(x, "start_time")
  axes <- paste(colnames(x), collapse = "")

  cat(sprintf(
    "<aclrtm_accelerometery>  [%d samples \u00d7 %d axes (%s)]",
    nrow(x), ncol(x), axes
  ))

  if (!is.null(sr)) {
    cat(sprintf("\n  sampling_rate : %.4g Hz", sr))
  }
  if (!is.null(st)) {
    cat(sprintf("\n  start_time : %s", format(st)))
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
head.aclrtm_accelerometery <- function(x, n = 6L, ...) {
  sr <- attr(x, "sampling_rate")
  st <- attr(x, "start_time")
  new_accelerometery(NextMethod(), sampling_rate = sr, start_time = st)
}

#' @export
tail.aclrtm_accelerometery <- function(x, n = 6L, ...) {
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
    st <- NULL
  }

  new_accelerometery(result, sampling_rate = sr, start_time = st)
}

#' @export
`[.aclrtm_accelerometery` <- function(x, i, j, ..., drop = FALSE) {

  sr <- attr(x, "sampling_rate")
  st <- attr(x, "start_time")

  # delegate subsetting to the matrix method
  result <- NextMethod(drop = drop)

  # if result is no longer a matrix (e.g. single row with drop = TRUE)
  # just return it
  if (!is.matrix(result)) return(result)

  # if no row subsetting, sampling_rate and start_time remain unchanged
  if (missing(i)) {
    return(new_accelerometery(result,
                              sampling_rate = sr,
                              start_time    = st))
  }

  # normalise i to positive integer indices to
  # handle logical, negative, and positive integer i
  n <- nrow(x)
  idx <- seq_len(n)[i]

  # if empty selection, early return
  if (length(idx) == 0L) {
    return(new_accelerometery(result,
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

  new_accelerometery(result, sampling_rate = sr, start_time = st)
}

#' @export
plot.aclrtm_accelerometery <- function(x, y, axes = "xyz", ...) {

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

  invisible(x)
}
