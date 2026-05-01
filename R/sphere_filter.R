#' Simple sphere filter
#'
#' Identify samples of a (triaxial) sensor data that describe a vector
#' higher than a set length.
#'
#' @param acc (Calibrated) sensor data with samples as rows.
#' @param r Radium of the sphere filter (max distance from the origin).
#' @return A Boolean vector that identifies all samples with a vector length
#'   lower than \code{r}
#' @export

sphere_filter <- function(acc, r = 1.05) {
  row_ss <- acc[, 1]^2
  for (i in seq_len(ncol(acc))[-1]) row_ss <- row_ss + acc[, i]^2
  sqrt(row_ss) < r
}
