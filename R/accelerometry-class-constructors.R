#' Low-level constructor for accelerometry objects
#'
#' Wraps a numeric matrix in the "aclrtm_accelerometry" class and attaches
#' optional metadata attributes.  Intended for internal use or for developers;
#' for interactive use, call [accelerometry()].
#'
#' @param x A numeric matrix with 1-3 columns whose names are a subset of
#'   \code{c("x","y","z")}.
#' @param sampling_rate \code{NULL} or a single positive numeric (Hz).
#' @param start_time \code{NULL} or a length-1 \code{POSIXct} vector.
#' @returns An object of class \code{c("aclrtm_accelerometry","matrix","array")}.
new_accelerometry <- function(x,
                              sampling_rate = NULL,
                              start_time    = NULL) {
  stopifnot(is.matrix(x), is.numeric(x))

  structure(
    x,
    sampling_rate = sampling_rate,
    start_time    = start_time,
    class         = c("aclrtm_accelerometry", "matrix", "array")
  )
}

#' Validate an accelerometry object
#'
#' Checks that \code{x} satisfies all structural invariants of the
#' \code{aclrtm_accelerometry} class and throws informative errors if
#' it does not.
#'
#' @param x Object to validate.
#' @returns \code{x}
#' @export

validate_accelerometry <- function(x) {

  # class & type
  if (!inherits(x, "aclrtm_accelerometry")) {
    stop("`x` must be an 'aclrtm_accelerometry' object.")
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("An 'aclrtm_accelerometry' object must be a numeric matrix.")
  }

  # dimensions
  nc <- ncol(x)
  if (is.null(nc) || nc < 1L || nc > 3L) {
    stop(
      "An 'aclrtm_accelerometry' matrix must have 1, 2, or 3 columns ",
      "(found ", nc, ")."
    )
  }

  # column names
  cn <- colnames(x)
  if (!is.null(cn)) {
    bad <- setdiff(cn, c("x", "y", "z"))
    if (length(bad) > 0L) {
      stop(
        "Column names must be a subset of c('x','y','z'). ",
        "Invalid name(s): ", paste(bad, collapse = ", "), "."
      )
    }
    if (anyDuplicated(cn)) {
      stop("Column names must be unique; found duplicates: ",
           paste(cn[duplicated(cn)], collapse = ", "), ".")
    }
  }

  # sampling_rate
  sr <- attr(x, "sampling_rate")
  if (!is.null(sr)) {
    if (!is.numeric(sr) || length(sr) != 1L || !is.finite(sr) || sr <= 0) {
      stop(
        "`sampling_rate` must be a single positive finite number; ",
        "got: ", deparse(sr), "."
      )
    }
  }

  # start_time
  st <- attr(x, "start_time")
  if (!is.null(st)) {
    if (length(st) != 1L) {
      stop(
        "`start_time` must be a single POSIXct or positive numeric value; ",
        "got: ", deparse(st), "."
      )
    }
    if (!inherits(st, "POSIXct") && (!is.numeric(st) || st < 0)) {
      stop(
        "`start_time` must be a single POSIXct or positive numeric value; ",
        "got: ", deparse(st), "."
      )
    }
  }

  x
}

#' Create an accelerometry object
#'
#' Constructor for \code{aclrtm_accelerometry} objects.
#' An \code{aclrtm_accelerometry} object
#' represents a time series of acceleration measurements. It is effectively
#' a matrix with rows being observations and each column being a spatial axis
#' (a subset of x, y, z). Two optional metadata attributes can be attached:
#' sampling_rate, a positive numeric giving the acquisition frequency in Hz,
#' and start_time, a POSIXct or numeric marking the beginning of the recording.
#'
#' @param x A numeric vector, matrix, or data frame of acceleration samples.
#' @param axes A string composed of any non-repeating combination of "x", "y",
#'   "z" (e.g. "xyz", "xz", "y"). If \code{NULL} (default), an attempt
#'   to estimate axes from the data structure will be performed.
#' @param sampling_rate Numeric scalar (Hz), or \code{NULL} (default).
#' @param start_time Either a positive real number, a \code{POSIXct} scalar,
#' or \code{NULL} (default).
#' @param ... Additional arguments passed to the appropriate method.
#' @param time_col Optional character vector of length 1 specifying the
#'   time column. It can also be an integer specifying the column index.
#' @details \code{aclrtm_accelerometry} is a utility class that can be used to
#'   work with (tri)axial accelerometry data in a consistent fashion within the
#'   \code{acceleratum} package. The constructor function
#'   \code{accelerometry()} can accept input data (argument \code{x})
#'   as either a numeric vector, a numeric matrix, or a data frame.
#'
#'   If \code{x} is a vector and the argument \code{axes} is provided,
#'   \code{x} is assumed to be ordered in multiplets according to the number
#'   of axes specified. E.g. given a vector of length \eqn{n \cdot 3} and
#'   \code{axes = "xyz"}, the input vector will be assumed to be arranged as
#'   \eqn{\left(x_1, y_1, z_1, x_2, y_2, z_2, \ldots, x_n, y_n, z_n\right)}.
#'
#'   If \code{x} is a matrix, each row is assumed to be an observation and
#'   columns are assumed to be the axes. Effectively, when the input is
#'   a matrix, \code{accelerometry()} will run a few sanity checks and
#'   return the same data structure with attached attributes (if provided).
#'   If the argument \code{axes} was provided and column names were
#'   already present, \code{axes} takes precedence and throws a warning.
#'
#'   If \code{x} is a data frame and argument \code{time_col} has been provided,
#'   \code{start_time} and \code{sampling_rate} will be estimated from
#'   \code{time_col}. \code{time_col} also takes precedence over
#'   \code{start_time} (if provided), but not over \code{sampling_rate}
#'   (if provided). If more than 3 numeric columns (excluding \code{time_col})
#'   are present in \code{x}, the argument \code{axes} must be provided and
#'   those columns must be present in \code{x}.
#'   With 3 or fewer numeric columns, column names will be retained if they are
#'   a subset of \code{c("x", "y", "z")} or overridden (with warning).
#'
#'   As accelerometry data can often be large, an object of class
#'   \code{"aclrtm_accelerometry"} does not store the vector of timestamps.
#'   Instead, attributes \code{start_time} and \code{sampling_rate} can be
#'   provided and timestamps will be estimated as needed within the
#'   \code{acceleratum} package. We feel that this is a good enough
#'   approach for the purposes of this package. That is: calibration and
#'   standardization of observations, tasks that can be performed
#'   without knowledge of timestamps. However, this comes at the cost of
#'   ignoring possible timing drift for the logger device. If that is of
#'   concern to the user, we recommend two approaches:
#'   \enumerate{
#'     \item If memory allows, most functions within the \code{acceleratum}
#'       package work with data frames as input and a time_col specification
#'       (as needed).
#'     \item It might be preferred to store the timestamp vector elsewhere
#'       to have it re-added to the data once most operations are performed.
#'   }
#'
#' @returns An object of class \code{c("aclrtm_accelerometry","matrix","array")}.
#' @examples
#' # vector as input
#' x <- 1:12
#' accelerometry(x,                           # 3 axes inferred
#'               sampling_rate = 4,
#'               start_time = as.POSIXct(0))
#' accelerometry(x, "xyz", 4, 1)              # 3 axes specified
#' accelerometry(x, "yz", 8, as.POSIXct(1e7)) # 2 axes specified
#'
#' ## discrepancy between vector length and axes argument throws an error
#' tryCatch(
#'   accelerometry(1:10, axes = "xyz"),
#'   error = function(e) cat(conditionMessage(e), "\n")
#' )
#'
#' # Matrix as input
#' x <- matrix(1:12, 4, 3)
#' accelerometry(x)                           # column names inferred
#' accelerometry(x, "yzx")                    # column names specified
#' colnames(x) <- letters[1:3]
#' accelerometry(x)                           # warning thrown
#'
#' # Data frame as input
#' x <- data.frame(x = 1:4,
#'                 y = 5:8,
#'                 z = 9:12,
#'                 timestamp = as.POSIXct(seq(.5,2,.5), tz = "UTC"))
#'
#' ## start_time and sampling_rate inferred from the timestamp column
#' accelerometry(x, time_col = "timestamp")
#'
#' ## if sampling rate is known from device settings,
#' ## providing it manually overrides inference.
#' accelerometry(x, sampling_rate = 5, time_col = "timestamp")
#'
#' @export

accelerometry <- function(x,
                          axes          = NULL,
                          sampling_rate = NULL,
                          start_time    = NULL,
                          ...) {
  UseMethod("accelerometry")
}

#' @export
#' @rdname accelerometry

accelerometry.numeric <- function(x,
                                  axes          = NULL,
                                  sampling_rate = NULL,
                                  start_time    = NULL,
                                  ...) {

  # resolve axes
  if (is.null(axes)) {
    col_names <- .infer_axes_from_length(length(x))
    message(
      "Note: `axes` not supplied; inferred axes = \"",
      paste(col_names, collapse = ""), "\"."
    )
  } else {
    col_names <- .parse_axes(axes)
  }

  n_cols <- length(col_names)

  # check divisibility
  if (length(x) %% n_cols != 0L) {
    stop(
      "Vector length (", length(x), ") is not divisible by the number of ",
      "axes (", n_cols, " for axes = \"",
      paste(col_names, collapse = ""), "\")."
    )
  }

  # reshape to matrix
  mat <- matrix(x, ncol = n_cols, byrow = TRUE,
                dimnames = list(NULL, col_names))

  validate_accelerometry(
    new_accelerometry(mat,
                      sampling_rate = sampling_rate,
                      start_time    = start_time)
  )
}

#' @export
#' @rdname accelerometry

accelerometry.matrix <- function(x,
                                 axes          = NULL,
                                 sampling_rate = NULL,
                                 start_time    = NULL,
                                 ...) {

  if (!is.numeric(x)) stop("`x` must be a numeric matrix.")

  nc <- ncol(x)
  if (is.null(nc) || nc < 1L || nc > 3L) {
    stop("`x` must have 1, 2, or 3 columns (found ", nc, ").")
  }

  # resolve column names
  existing <- colnames(x)
  if (!is.null(axes)) {
    # user supplied axes - use them, check length matches
    col_names <- .parse_axes(axes)
    if (length(col_names) != nc) {
      stop(
        "`axes` specifies ", length(col_names), " axis/axes but `x` has ",
        nc, " column(s)."
      )
    }
    if (!is.null(existing) && !all(col_names == existing)) {
      warning(
        "Existing column names have been overridden by `axes`. ",
        "If that is undesired, rerun with `axes = NULL`"
      )
    }
  } else {
    # no axes supplied: try existing column names first
    if (!is.null(existing) && all(existing %in% c("x", "y", "z"))) {
      col_names <- existing
    } else {
      if (!is.null(existing) && length(existing) > 0L) {
        warning(
          "Existing column names (", paste(existing, collapse = ", "),
          ") are not valid axis labels and have been replaced."
        )
      }
      col_names <- .infer_axes(nc)
      message(
        "Note: `axes` not supplied; inferred axes = \"",
        paste(col_names, collapse = ""), "\"."
      )
    }
  }

  colnames(x) <- col_names

  validate_accelerometry(
    new_accelerometry(x,
                      sampling_rate = sampling_rate,
                      start_time    = start_time)
  )
}

#' @export
#' @rdname accelerometry

accelerometry.data.frame <- function(x,
                                     axes          = NULL,
                                     sampling_rate = NULL,
                                     start_time    = NULL,
                                     time_col      = NULL,
                                     ...) {

  if (!is.data.frame(x)) stop("`x` must be a data frame.")

  # handle time column
  if (!is.null(time_col)) {

    # resolve time_col to a column name
    if (is.numeric(time_col)) {
      if (time_col < 1L || time_col > ncol(x)) {
        stop("`time_col` index (", time_col, ") is out of range [1, ",
             ncol(x), "].")
      }
      time_col <- names(x)[time_col]
    } else if (!is.character(time_col) || length(time_col) != 1L) {
      stop("`time_col` must be a single column name or integer index.")
    } else if (!time_col %in% names(x)) {
      stop("Column '", time_col, "' not found in `x`.")
    }

    t_vals <- x[[time_col]]

    t_sec <- as.numeric(t_vals)
    diffs <- diff(t_sec)

    if (any(diffs < 0)) {
      message(
        "Note: input data is not ordered by '", time_col, "' and it will ",
        "be rearranged."
      )
      x <- x[order(t_sec),]
      t_vals <- x[[time_col]]
    }

    # validate time column type
    if (!is.numeric(t_vals) && !inherits(t_vals, "POSIXct") &&
        !inherits(t_vals, "POSIXlt")) {
      stop(
        "The time column ('", time_col, "') must be numeric or POSIXct/POSIXlt."
      )
    }

    # derive start_time if not explicitly provided
    if (!is.null(start_time)) {
      message(
        "Note: the argument `start_time` is ignored when `time_col` is ",
        "defined and it will be set to the lowest value of `time_col`."
      )
    }
    st_raw <- min(t_vals, na.rm = TRUE)

    if (inherits(t_vals, c("POSIXct", "POSIXlt"))) {
      start_time <- as.POSIXct(st_raw, tz = attr(t_vals, "tzone") %||% "UTC")
    } else {
      start_time <- st_raw
    }

    # estimate sampling_rate from time differences if not explicitly provided
    if (is.null(sampling_rate)) {

      t_sec <- as.numeric(t_vals)          # works for both numeric and POSIXct
      diffs <- diff(t_sec)
      if (all(diffs == 0L)) {
        warning(
          "Cannot estimate `sampling_rate`: all time differences are zero."
        )
      } else {
        sampling_rate <- 1 / stats::median(diffs)
        message(
          sprintf("Note: `sampling_rate` estimated as %.4g Hz.", sampling_rate)
        )
      }
    }

    # drop the time column from the data frame before continuing
    x <- x[, setdiff(names(x), time_col), drop = FALSE]
  }

  # select acceleration columns
  if (!is.null(axes)) {
    col_names <- .parse_axes(axes)
    missing_cols <- setdiff(col_names, names(x))
    if (length(missing_cols) > 0L) {
      stop(
        "Column(s) specified in `axes` not found in `x`: ",
        paste(missing_cols, collapse = ", "), "."
      )
    }
    acc_df <- x[, col_names, drop = FALSE]
  } else {
    # fall back to all remaining numeric columns
    num_idx <- vapply(x, is.numeric, logical(1L))
    acc_df  <- x[, num_idx, drop = FALSE]

    if (ncol(acc_df) == 0L) {
      stop("No numeric columns remain after removing the time column.")
    }
    if (ncol(acc_df) > 3L) {
      stop(
        "Found ", ncol(acc_df), " numeric columns but an `aclrtm_accelerometry`",
        "object can have at most 3.  Use `axes` to select the relevant columns."
      )
    }
  }

  # check all selected columns are numeric
  non_num <- names(acc_df)[!vapply(acc_df, is.numeric, logical(1L))]
  if (length(non_num) > 0L) {
    stop(
      "The following column(s) are not numeric: ",
      paste(non_num, collapse = ", "), "."
    )
  }

  mat <- as.matrix(acc_df)

  # delegate to matrix method for column-name inference / validation
  accelerometry.matrix(mat,
                       axes          = axes,
                       sampling_rate = sampling_rate,
                       start_time    = start_time)
}
