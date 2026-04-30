#' Apply calibration
#'
#' Apply calibration to raw data
#'
#' @param A_sens Input accelerometry data with axes as columns and samples as
#'   rows. See details.
#' @param S The scale factor diagonal matrix (can be provided as a vector)
#' @param M The misalignment matrix
#' @param R A rotation matrix to account for a sensor mounted "flipped"
#'   around one or two axes (e.g. a collar mounted backwards).
#'   Can be provided as diagonal matrix or numeric vector of 1s and -1s,
#'   or a character vector or character string of axes to flip around.
#' @param b The bias vector
#' @return An \code{accelerometry} object containing the calibrated data
#'   \eqn{A_{true}^T}, with the same dimensions and - if available -
#'   sampling rate and start time as the input \code{A_sens}.
#' @details
#' The sensor data model used is:
#'
#' \deqn{A_{sens} = S M R A_{true} + b}
#'
#' where \eqn{A_{sens}} is the column vector of raw signal recorded by the
#' sensor, \eqn{A_{true}} is the true acceleration in the orthogonal reference
#' frame, \eqn{S} is a diagonal matrix of scale factors, \eqn{M} is the
#' misalignment matrix, \eqn{R} is a rotation matrix accounting for 180 degrees
#' rotations of the sensor around one or two axes,
#' and \eqn{b} is the bias vector. Solving for \eqn{A_{true}}
#' gives the calibration equation:
#'
#' \deqn{A_{true} = R^{-1} M^{-1} S^{-1} (A_{sens} - b)}
#' In the triaxial case and with \eqn{n} samples, \eqn{A_{true}} and
#' \eqn{A_{sens}} are matrices of dimensions \eqn{3 \times n}, where rows
#' represent axes (x, y, z) and columns represent samples. However,
#' as is customary, data is expected with samples as rows and axes as
#' columns (i.e. \eqn{A_{sens}^T}). The function therefore solves for
#' \eqn{A_{true}} in its transposed form:
#'
#' \deqn{A_{true}^T = (A_{sens}^T - b^T) S^{-1} (M^{-1})^T R}
#'
#' @export
apply_cal <- function(A_sens,
                      S = diag(ncol(A_sens)),
                      M = diag(ncol(A_sens)),
                      R = diag(ncol(A_sens)),
                      b = rep(0, ncol(A_sens))) {

  # Checks for A_sens
  if (!is.matrix(A_sens) || !is.numeric(A_sens)) {
    stop("`A_sens` must be a numeric matrix.", call. = FALSE)
  }
  nc <- ncol(A_sens)
  if (nc < 1L || nc > 3L) {
    stop("`A_sens` must have 1, 2, or 3 columns (found ", nc, ").",
         call. = FALSE)
  }

  # Checks for S
  if (is.numeric(S) && is.vector(S)) {
    if (length(S) != nc) {
      stop("`S` as a vector must have length equal to the number of columns ",
           "in `A_sens` (", nc, "); got length ", length(S), ".",
           call. = FALSE)
    }
    S <- diag(S)
  }

  if (!is.matrix(S) || !is.numeric(S)) {
    stop("`S` must be a numeric matrix or numeric vector.", call. = FALSE)
  }

  if (any(dim(S) != nc)) {
    stop("`S` must be a ", nc, "\u00d7", nc, " matrix to match `A_sens`; ",
         "got ", nrow(S), "\u00d7", ncol(S), ".", call. = FALSE)
  }

  if (abs(det(S)) < sqrt(.Machine$double.eps)) {
    stop("`S` is singular or near-singular and cannot be inverted.",
         call. = FALSE)
  }

  # Checks for M
  if (!is.matrix(M) || !is.numeric(M)) {
    stop("`M` must be a numeric matrix.", call. = FALSE)
  }
  if (any(dim(M) != nc)) {
    stop("`M` must be a ", nc, "\u00d7", nc, " matrix to match `A_sens`; ",
         "got ", nrow(M), "\u00d7", ncol(M), ".", call. = FALSE)
  }
  if (abs(det(M)) < sqrt(.Machine$double.eps)) {
    stop("`M` is singular or near-singular and cannot be inverted.",
         call. = FALSE)
  }

  # Checks for R
  if (is.character(R)) {
    if(nc < 3L) {
      stop("Rotation correction via character input is only supported for ",
           "triaxial data. For biaxial or uniaxial data, provide `R` as a ",
           "rotation matrix directly.", call. = FALSE)
    }
    if (length(R) == 1L && nchar(R) > 1L) {
      flip_around <- .parse_axes(R)
    } else {
      # character vector of axis names e.g. c("x", "y")
      bad <- R[!R %in% c("x", "y", "z")]
      if (length(bad) > 0L) {
        stop("`R` contains invalid axis label(s): ",
             paste(bad, collapse = ", "), ".", call. = FALSE)
      }
      flip_around <- R
    }
    # check axes are present in A_sens
    missing_axes <- setdiff(flip_around, colnames(A_sens))
    if (length(missing_axes) > 0L) {
      stop("Axis/axes to flip around not found in `A_sens`: ",
           paste(missing_axes, collapse = ", "), ".", call. = FALSE)
    }

    if (length(flip_around) == 3L) {
      message("Input for R attempts to rotate around all axes. Note that ",
              "this is the same as performing no rotation at all.")
      R <- diag(nc)
    } else {
        if (length(flip_around) == 2L) {
        axes_to_flip <- flip_around
      } else {
        axes_to_flip <- setdiff(colnames(A_sens), flip_around)
      }

      diag_vals <- rep(1, nc)
      names(diag_vals) <- colnames(A_sens)
      diag_vals[axes_to_flip] <- -1L
      R <- diag(diag_vals)
    }
  }

  if (is.numeric(R) && is.vector(R)) {
    # numeric vector: treat as diagonal entries
    if (length(R) != nc) {
      stop("`R` as a vector must have length equal to the number of columns ",
           "in `A_sens` (", nc, "); got length ", length(R), ".", call. = FALSE)
    }
    if (!all(abs(R) == 1L)) {
      stop("`R` as a diagonal vector must contain only 1 or -1.", call. = FALSE)
    }
    R <- diag(R)
  }

  if (is.matrix(R)) {
    if (!is.numeric(R)) {
      stop("`R` must be a numeric matrix.", call. = FALSE)
    }
    if (any(dim(R) != nc)) {
      stop("`R` must be a ", nc, "\u00d7", nc, " matrix; got ",
           nrow(R), "\u00d7", ncol(R), ".", call. = FALSE)
    }
    # check it is a proper rotation matrix
    if (det(R) != 1) {
      stop("`R` is not a proper rotation matrix (det(R) must equal 1).",
           call. = FALSE)
    }
    if (!all((R %*% t(R)) == diag(nc))) {
      stop("`R` is not a proper rotation matrix as it is not orthogonal ",
           "(R %*% t(R) must equal the identity matrix).",
           call. = FALSE)
    }
  } else {
    stop("`R` must be a character string, character vector, numeric vector, ",
         "or numeric matrix.", call. = FALSE)
  }

  # Checks for b
  if (!is.numeric(b)) {
    stop("`b` must be a numeric vector or 1-row matrix.", call. = FALSE)
  }
  b <- as.vector(b)
  if (length(b) != nc) {
    stop("`b` must have length equal to the number of columns in `A_sens` (",
         nc, "); got length ", length(b), ".", call. = FALSE)
  }

  sr <- attr(A_sens, "sampling_rate")
  st <- attr(A_sens, "start_time")

  A_true <- sweep(A_sens, 2L, b, "-") %*% solve(S) %*% t(solve(M)) %*% R
  colnames(A_true) <- colnames(A_sens)

  new_accelerometry(A_true,
                    sampling_rate = sr,
                    start_time    = st)
}
