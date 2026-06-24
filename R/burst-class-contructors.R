#' Low-level constructor for a aclrtm_burst object
#'
#' @param df        A data.frame whose `data_col` column is already a list of
#'                  numeric vectors.
#' @param data_col  Character scalar — name of the burst-data column.
#' @param ts_col    Character scalar or NULL — name of the timestamp column.
#'
#' @return An object of class c("aclrtm_burst", "data.frame").
new_burst <- function(df, data_col, ts_col = NULL) {
  stopifnot(is.data.frame(df))
  stopifnot(is.character(data_col), length(data_col) == 1L)
  stopifnot(is.null(ts_col) || (is.character(ts_col) && length(ts_col) == 1L))

  structure(
    df,
    class    = c("aclrtm_burst", "data.frame"),
    data_col = data_col,
    ts_col   = ts_col
  )
}

#' Validate a aclrtm_burst object
#'
#' @param x  An object to validate.
#'
#' @return `x` invisibly if valid.
#' @export
validate_burst <- function(x) {

  # class
  if (!inherits(x, "aclrtm_burst")) {
    stop("`x` must inherit from \"aclrtm_burst\".", call. = FALSE)
  }
  if (!is.data.frame(x)) {
    stop("`x` must also be a data.frame.", call. = FALSE)
  }

  # required attributes
  data_col <- attr(x, "data_col")
  ts_col   <- attr(x, "ts_col")

  if (is.null(data_col) || !is.character(data_col) || length(data_col) != 1L) {
    stop("Attribute `data_col` must be a single character string.", call. = FALSE)
  }

  if (!is.null(ts_col) &&
      (!is.character(ts_col) || length(ts_col) != 1L)) {
    stop("Attribute `ts_col` must be NULL or a single character string.",
         call. = FALSE)
  }

  # column presence
  if (!data_col %in% names(x)) {
    stop(
      sprintf("Column \"%s\" (data_col) not found in the data.frame.", data_col),
      call. = FALSE
    )
  }

  if (!is.null(ts_col) && !ts_col %in% names(x)) {
    stop(
      sprintf("Column \"%s\" (ts_col) not found in the data.frame.", ts_col),
      call. = FALSE
    )
  }

  # burst column type
  burst_col <- x[[data_col]]

  if (!is.list(burst_col)) {
    stop(
      sprintf("Column \"%s\" must be a list, not \"%s\".",
              data_col, class(burst_col)[[1L]]),
      call. = FALSE
    )
  }

  non_numeric <- which(!vapply(burst_col, is.numeric, logical(1L)))
  if (length(non_numeric) > 0L) {
    stop(
      sprintf(
        "All elements of \"%s\" must be numeric vectors. Non-numeric at row(s): %s.",
        data_col,
        paste(non_numeric, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # ts_col uniqueness
  if (!is.null(ts_col)) {
    ts_vals <- x[[ts_col]]
    if (anyDuplicated(ts_vals) != 0L) {
      warning(
        sprintf("Timestamp column \"%s\" contains duplicate values.", ts_col),
        call. = FALSE
      )
    }
  }

  x
}

#' Create a aclrtm_burst object
#'
#' @param df        A data.frame.
#' @param data_col  Name of the column that holds burst data.  The column may
#'                  be a list of numeric vectors (used as-is) or a character
#'                  vector of `sep`-separated numbers (parsed automatically).
#' @param ts_col    Optional name of a timestamp column.  Stored as the `ts_col`
#'                  attribute; the column itself is not modified.
#' @param sep       Separator used when parsing character burst data
#'                  (default `" "`).
#'
#' @return An `aclrtm_burst` object.
#' @export
burst <- function(df, data_col, ts_col = NULL, sep = " ") {

  # input sanity
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.", call. = FALSE)
  }
  if (!is.character(data_col) || length(data_col) != 1L) {
    stop("`data_col` must be a single character string.", call. = FALSE)
  }
  if (!data_col %in% names(df)) {
    stop(sprintf("Column \"%s\" not found in `df`.", data_col), call. = FALSE)
  }
  if (!is.null(ts_col)) {
    if (!is.character(ts_col) || length(ts_col) != 1L) {
      stop("`ts_col` must be NULL or a single character string.", call. = FALSE)
    }
    if (!ts_col %in% names(df)) {
      stop(sprintf("Column \"%s\" (ts_col) not found in `df`.", ts_col),
           call. = FALSE)
    }
  }

  # parse burst column if needed
  raw <- df[[data_col]]

  if (is.list(raw)) {
    # Already a list — coerce every element to numeric just in case
    parsed <- lapply(raw, function(el) {
      if (is.numeric(el)) el else as.numeric(el)
    })
  } else if (is.character(raw)) {
    # CSV-style strings → list of numeric vectors
    parsed <- lapply(
      strsplit(raw, sep, fixed = TRUE),
      function(parts) as.numeric(trimws(parts))
    )
  } else if (is.numeric(raw)) {
    # A plain numeric column — wrap each scalar in a length-1 list
    parsed <- as.list(raw)
  } else {
    stop(
      sprintf(
        "Column \"%s\" must be a list, character, or numeric vector.",
        data_col
      ),
      call. = FALSE
    )
  }

  df[[data_col]] <- parsed

  # build & validate
  obj <- new_burst(df, data_col = data_col, ts_col = ts_col)
  validate_burst(obj)
}
