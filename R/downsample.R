#' Downsample accelerometry data
#'
#' Description
#'
#' @param acc Accelerometry data input in matrix-like or data.frame format
#' @param ... passed to methods
#' @param from_sr Sampling rate of data input. Can be NULL.
#' @param to_sr Sampling rate of output. Must be lower than, and be a
#'   divisor of, the sampling rate of the input.
#' @param ds_factor Downsampling factor. Alternative to using
#'   \code{to_sr} and \code{from_sr}.
#' @param FUN An optional function to apply to aggregated data rather than
#'   simple downsampling.
#' @param ... Arguments passed to \code{FUN}
#'
#' @returns An object of class \code{c("aclrtm_accelerometry","matrix","array")}.
#' @export
downsample <- function(acc, ...) {
  UseMethod("downsample")
}

#' @rdname downsample
#' @export
downsample.matrix <- function(acc, to_sr = NULL, from_sr = NULL,
                              ds_factor = NULL, FUN = NULL, ...) {
  acc_sr <- attr(acc, "sampling_rate")
  acc_st <- attr(acc, "start_time")

  # check input
  if (!is.null(ds_factor) && (!is.null(to_sr) || !is.null(from_sr))) {
    stop("If argument `ds_factor` is defined, `to_sr` and `from_sr` must be ",
         "NULL",
         call. = FALSE)
  }
  if (!is.null(ds_factor)) {
    if (!is.numeric(ds_factor) || length(ds_factor) != 1L ||
        ds_factor != round(ds_factor) || ds_factor < 1L) {
      stop("`ds_factor` must be a single positive integer.", call. = FALSE)
    }
  }
  if (is.null(ds_factor)) {
    if (is.null(to_sr)) {
      stop("Please provide at least one between `to_sr` and `ds_factor`",
           call. = FALSE)
    }
    if (is.null(acc_sr) && is.null(from_sr)) {
      stop("Sampling rate input of data is neither embedded nor provided. ",
           "Please provide `from_sr`",
           call. = FALSE)
    }
    if (!is.null(acc_sr) && !is.null(from_sr) && acc_sr != from_sr) {
      stop("Input data sampling rate (",
           acc_sr, " Hz",
           ") differs from input `from_sr` (",
           from_sr, " Hz", ")",
           call. = FALSE)
    }
    if (!is.null(acc_sr)) {
      from_sr <- acc_sr
    }
    divisor_safe_check <- (abs(round(from_sr / to_sr) - from_sr / to_sr) >
      sqrt(.Machine$double.eps))
    if (to_sr >= from_sr || divisor_safe_check) {
      stop("`to_sr` must be lower than, and be a divisor of ",
           "input sampling rate\n",
           "Input sampling rate: ", from_sr, " Hz\n",
           "Requested output rate: ", to_sr, " Hz",
           call. = FALSE)
    }

    ds_factor <- trunc(from_sr / to_sr)
  } else if (!is.null(acc_sr)) {
    ds_factor <- trunc(ds_factor)
    to_sr <- acc_sr / ds_factor
  }

  if (is.null(FUN)) {
    # Subsample and return
    out_i <- seq(1, nrow(acc), by = ds_factor)

    accelerometry(suppressMessages(acc[out_i,]),
                  sampling_rate = to_sr,
                  start_time = acc_st)
  } else {
    FUN <- match.fun(FUN)
    nend <- floor(nrow(acc)/ds_factor) * ds_factor

    cn <- colnames(acc)

    out_m <- apply(array(c(acc[1:nend, ]),
                         dim = c(ds_factor, nend/ds_factor, ncol(acc))),
                   MARGIN = c(2L, 3L), FUN = FUN, ...)

    colnames(out_m) <- cn

    accelerometry(out_m,
                  sampling_rate = to_sr,
                  start_time = acc_st)
  }

}
