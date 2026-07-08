#' Low-level constructor for a aclrtm_burst object
#'
#' @param df        A data.frame whose `data_col` column is already a list of
#'                  numeric vectors.
#' @param axes A string composed of any non-repeating combination of "x", "y",
#'   "z" (e.g. "xyz", "xz", "y"). If \code{NULL} (default), an attempt
#'   to estimate axes from the data structure will be performed.
#' @param data_col  Character scalar — name of the burst-data column.
#' @param ts_col    Character scalar or NULL — name of the timestamp column.
#'
#' @return An object of class c("aclrtm_burst", "data.frame").
new_burst <- function(df, data_col, axes, ts_col = NULL) {
  stopifnot(is.data.frame(df))
  stopifnot(is.character(data_col), length(data_col) == 1L)
  stopifnot(is.character(axes),     length(axes)     == 1L)
  stopifnot(is.null(ts_col) || (is.character(ts_col) && length(ts_col) == 1L))

  structure(
    df,
    class    = c("aclrtm_burst", "data.frame"),
    data_col = data_col,
    axes     = axes,
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
    stop("`x` must be a data.frame.", call. = FALSE)
  }

  # required attributes
  data_col <- attr(x, "data_col")
  axes     <- attr(x, "axes")
  ts_col   <- attr(x, "ts_col")

  if (is.null(data_col) || !is.character(data_col) || length(data_col) != 1L) {
    stop("Attribute `data_col` must be a single character string.", call. = FALSE)
  }
  if (is.null(axes) || !is.character(axes) || length(axes) != 1L) {
    stop("Attribute `axes` must be a single character string (e.g. \"xyz\").",
         call. = FALSE)
  }
  if (!is.null(ts_col) && (!is.character(ts_col) || length(ts_col) != 1L)) {
    stop("Attribute `ts_col` must be NULL or a single character string.",
         call. = FALSE)
  }

  # column presence
  if (!data_col %in% names(x)) {
    stop(
      sprintf("Column \"%s\" (data_col) not found in `x`", data_col),
      call. = FALSE
    )
  }

  if (!is.null(ts_col) && !ts_col %in% names(x)) {
    stop(
      sprintf("Column \"%s\" (ts_col) not found in `x`", ts_col),
      call. = FALSE
    )
  }

  # axes validity
  col_names <- .parse_axes(axes)
  n_axes    <- length(col_names)

  # burst column
  burst_col <- x[[data_col]]

  if (!is.list(burst_col)) {
    stop(
      sprintf("Column \"%s\" must be a list, not \"%s\".",
              data_col, class(burst_col)[[1L]]),
      call. = FALSE
    )
  }

  n_bursts      <- length(burst_col)
  empty_idx     <- integer(0L)
  non_matrix    <- integer(0L)
  non_numeric   <- integer(0L)
  wrong_ncol    <- integer(0L)
  ragged_nrow   <- integer(0L)   # bursts whose nrow differs from the first non-empty
  ref_nrow      <- NULL          # nrow of the first non-empty burst

  for (i in seq_len(n_bursts)) {
    b <- burst_col[[i]]

    # empty burst
    if ((is.matrix(b) && nrow(b) == 0L) || length(b) == 0L) {
      empty_idx <- c(empty_idx, i)
      next
    }

    # must be a matrix
    if (!is.matrix(b)) {
      non_matrix <- c(non_matrix, i)
      next          # can't check further without a matrix
    }

    # must be numeric
    if (!is.numeric(b)) {
      non_numeric <- c(non_numeric, i)
    }

    # ncol must match number of axes — hard error collected below
    if (ncol(b) != n_axes) {
      wrong_ncol <- c(wrong_ncol, i)
    }

    # ragged nrow check
    if (is.null(ref_nrow)) {
      ref_nrow <- nrow(b)
    } else if (nrow(b) != ref_nrow) {
      ragged_nrow <- c(ragged_nrow, i)
    }
  }

  # hard errors first
  if (length(non_matrix) > 0L) {
    stop(
      sprintf(
        "All elements of \"%s\" must be matrices. Non-matrix at row(s): %s.",
        data_col, paste(non_matrix, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (length(non_numeric) > 0L) {
    stop(
      sprintf(
        "All burst matrices in \"%s\" must be numeric. Non-numeric at row(s): %s.",
        data_col, paste(non_numeric, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (length(wrong_ncol) > 0L) {
    stop(
      sprintf(
        paste0("Burst matrices in \"%s\" must have %d column(s) to match ",
               "axes = \"%s\". Wrong ncol at row(s): %s."),
        data_col, n_axes, axes, paste(wrong_ncol, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # warnings
  if (length(empty_idx) > 0L) {
    warning(
      sprintf(
        "Empty burst(s) found in \"%s\" at row(s): %s.",
        data_col, paste(empty_idx, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (length(ragged_nrow) > 0L) {
    warning(
      sprintf(
        paste0("Ragged bursts detected in \"%s\": row(s) %s have a different ",
               "number of samples than the first non-empty burst (%d samples)."),
        data_col, paste(ragged_nrow, collapse = ", "), ref_nrow
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
#' @param axes A string composed of any non-repeating combination of "x", "y",
#'   "z" (e.g. "xyz", "xz", "y"). If \code{NULL} (default), an attempt
#'   to estimate axes from the data structure will be performed.
#' @param ts_col    Optional name of a timestamp column.  Stored as the `ts_col`
#'                  attribute; the column itself is not modified.
#' @param sep       Separator used when parsing character burst data
#'                  (default `" "`).
#'
#' @return An `aclrtm_burst` object.
#' @export
burst <- function(df, data_col, axes = NULL, ts_col = NULL, sep = " ") {

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
    # could be matrices already or plain vectors — handle below
    vec_list <- lapply(raw, function(el) {
      if (is.matrix(el)) el          # pass through; checked later
      else if (is.numeric(el)) el
      else as.numeric(el)
    })
  } else if (is.character(raw)) {
    vec_list <- lapply(
      strsplit(raw, sep, fixed = TRUE),
      function(parts) as.numeric(trimws(parts))
    )
  } else if (is.numeric(raw)) {
    # scalar column — each value becomes a 1×1 matrix after reshape
    vec_list <- as.list(raw)
  } else {
    stop(
      sprintf("Column \"%s\" must be a list, character, or numeric vector.",
              data_col),
      call. = FALSE
    )
  }

  # resolve axes
  if (is.null(axes)) {

    vec_lengths <- vapply(vec_list, function(el) {
      if (is.matrix(el)) NA_integer_
      else                length(el)
    }, integer(1L))

    non_empty_lengths <- vec_lengths[!is.na(vec_lengths) & vec_lengths > 0L]

    if (length(non_empty_lengths) > 0L) {
      g         <- Reduce(.gcd2, non_empty_lengths)
      col_names <- .infer_axes_from_length(g)
    } else {
      # all elements are already matrices or all empty — derive from ncol
      mat_ncols <- vapply(vec_list, function(el) {
        if (is.matrix(el) && nrow(el) > 0L) ncol(el) else NA_integer_
      }, integer(1L))
      mat_ncols <- mat_ncols[!is.na(mat_ncols)]

      if (length(mat_ncols) > 0L) {
        col_names <- .infer_axes(mat_ncols[[1L]])
      } else {
        col_names <- "x"
        warning("All bursts are empty; defaulting to axes = \"x\".", call. = FALSE)
      }
    }

    axes <- paste(col_names, collapse = "")
    message(sprintf("Note: `axes` not supplied; inferred axes = \"%s\".", axes))

  } else {
    col_names <- .parse_axes(axes)
  }

  n_axes <- length(col_names)

  # vecs to matrices
  parsed <- lapply(vec_list, function(el) {
    if (is.matrix(el)) {
      # already a matrix: just (re)set column names if needed
      if (!identical(colnames(el), col_names)) colnames(el) <- col_names
      return(el)
    }
    # plain numeric vector
    n <- length(el)
    if (n == 0L) {
      # empty burst → 0-row matrix with correct columns
      return(matrix(numeric(0L), nrow = 0L, ncol = n_axes,
                    dimnames = list(NULL, col_names)))
    }
    if (n %% n_axes != 0L) {
      stop(
        sprintf(
          paste0("Burst vector length (%d) is not divisible by the number of ",
                 "axes (%d for axes = \"%s\")."),
          n, n_axes, axes
        ),
        call. = FALSE
      )
    }
    matrix(el, ncol = n_axes, byrow = TRUE, dimnames = list(NULL, col_names))
  })

  df[[data_col]] <- parsed

  # build & validate
  obj <- new_burst(df, data_col = data_col, axes = axes, ts_col = ts_col)
  validate_burst(obj)
}
