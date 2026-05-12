#' Simple sphere filter
#'
#' Identifies samples of a (triaxial) sensor data matrix that describe a vector
#' longer than a set length.
#'
#' @param acc (Calibrated) sensor data with samples as rows.
#' @param r Radium of the sphere filter (max distance from the origin).
#' @returns A Boolean vector that identifies all samples with a vector length
#'   lower than \code{r}
#' @export

sphere_filter <- function(acc, r = 1.05) {
  row_ss <- rowSums(acc^2)
  sqrt(row_ss) < r
}
