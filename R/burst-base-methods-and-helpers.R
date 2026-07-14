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
