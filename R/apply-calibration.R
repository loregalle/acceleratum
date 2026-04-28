#' Apply calibration
#'
#' Apply calibration to raw data
#'
#' @param A_sens Input accelerometry data
#' @param S The scale factor matrix (can be provided as a vector)
#' @param M The misalignment matrix
#' @param b The bias vector
#' @return An \code{accelerometry} object containing the calibrated data
#'   \eqn{A_{true}^T}, with the same dimensions and - if available -
#'   sampling rate and start time as the input \code{A_sens}.
#' @details
#' The accelerometer data model used follows:
#'
#' \deqn{A_{sens} = S M A_{true} + b}
#'
#' where \eqn{A_{sens}} is the column vector of raw signal recorded by the
#' sensor, \eqn{A_{true}} is the true acceleration in the orthogonal reference
#' frame, \eqn{S} is a diagonal matrix of scale factors, \eqn{M} is the
#' misalignment matrix, and \eqn{b} is the bias vector. Solving for
#' \eqn{A_{true}} gives the calibration equation:
#'
#' \deqn{A_{true} = M^{-1} S^{-1} (A_{sens} - b)}
#' In the triaxial case and with \eqn{n} samples, \eqn{A_{true}} and
#' \eqn{A_{sens}} are matrices of dimensions \eqn{3 \times n}, where rows
#' represent axes (x, y, z) and columns represent samples. However,
#' as is customary, data is expected with samples as rows and axes as
#' columns (i.e. \eqn{A_{sens}^T}). The function therefore solves for
#' \eqn{A_{true}} in its transposed form:
#'
#' \deqn{A_{true}^T = (A_{sens}^T - b^T) S^{-1} (M^{-1})^T}
#'
#' @export
apply_cal <- function(A_sens,
                      S = diag(ncol(A_sens)),
                      M = diag(ncol(A_sens)),
                      b = rep(0, ncol(A_sens))) {

  if (!is.matrix(A_sens) || !is.numeric(A_sens)) {
    stop("`A_sens` must be a numeric matrix.", call. = FALSE)
  }
  nc <- ncol(A_sens)
  if (nc < 1L || nc > 3L) {
    stop("`A_sens` must have 1, 2, or 3 columns (found ", nc, ").",
         call. = FALSE)
  }

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

  A_true <- sweep(A_sens, 2L, b, "-") %*% solve(S) %*% t(solve(M))
  colnames(A_true) <- colnames(A_sens)

  new_accelerometry(A_true,
                    sampling_rate = sr,
                    start_time    = st)
}
