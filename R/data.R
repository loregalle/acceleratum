#' Accelerometery Data for Bench Calibration
#'
#' Data from a triaxial accelerometer device mounted on a cubic frame and
#' rotated through 6 known orientations. Sampling rate is fixed at 8 Hz.
#'
#' @format A data frame with 29 rows and 3 columns:
#' \describe{
#'   \item{tag_local_identifier}{ID of the tag used}
#'   \item{individual_local_identifier}{fixed value to "cal" for "calibration"}
#'   \item{timestamp}{A \code{POSIXct} column of time stamps}
#'   \item{acceleration_axes}{A specification for the axes collected}
#'   \item{acceleration_sampling_frequency_per_axis}{Sampling frequency}
#'   \item{accelerations_raw}{Character vector of raw acceleration data in burst
#'     format (see details).}
#' }
#'
#' @details The `accelerations_raw` column stores the data in bursts as a
#' single space-separated character string, with axis measurements
#' interleaved in repeating x, y, z order:
#'
#' ```
#' "x1 y1 z1 x2 y2 z2 ... xn yn zn"
#' ```
#'
#' where `n` is the number of samples in that burst.
#'
#' @source Data collected by the package authors
"bench"

#' Annotation Example for Bench Calibration
#'
#' Segment annotations example for [bench] data.
#'
#' @format A data frame with 6 rows and 5 columns:
#' \describe{
#'   \item{label}{A label specifying the face orientation}
#'   \item{xmin}{The start timestamp (in seconds) of each segment}
#'   \item{xmax}{The end timestamp (in seconds) of each segment}
#'   \item{imin}{Row index of the start of each segment}
#'   \item{imax}{Row index of the end of each segment}
#' }
#'
#' @source Collected by the package authors interactively using
#'   [segment_select()]
"bench_annotations"

#' Muskox Accelerometery Data
#'
#' About 4 days of triaxial accelerometery data from a device mounted on a
#' collar and deployed on a muskox in Zackenberg Research Station, Greenland.
#' Sampling rate is fixed at 8 Hz.
#'
#' @format A data frame with 69103 rows and 6 columns:
#' \describe{
#'   \item{tag_local_identifier}{ID of the tag used}
#'   \item{individual_local_identifier}{ID of the animal}
#'   \item{timestamp}{A \code{POSIXct} column of time stamps}
#'   \item{acceleration_axes}{A specification for the axes collected}
#'   \item{acceleration_sampling_frequency_per_axis}{Sampling frequency}
#'   \item{accelerations_raw}{Character vector of raw acceleration data in burst
#'     format (see details).}
#' }
#'
#' @details The `accelerations_raw` column stores the data in bursts as a
#' single space-separated character string, with axis measurements
#' interleaved in repeating x, y, z order:
#'
#' ```
#' "x1 y1 z1 x2 y2 z2 ... xn yn zn"
#' ```
#'
#' where `n` is the number of samples in that burst.
#'
#' @source Data collected by the package authors
"muskox"
