# Generals ----
#' Parse a axes string into a character vector of axis labels
#'
#' @param axes A string composed of any non-repeating combination of "x", "y",
#'   "z" (e.g. "xyz", "xz", "y").
#' @return A character vector such as c("x","y","z").
#' @noRd

.parse_axes <- function(axes) {
  if (!is.character(axes) || length(axes) != 1L) {
    stop("`axes` must be a single character string (e.g. \"xyz\", \"xz\").")
  }
  chars <- strsplit(axes, "")[[1L]]
  if (length(chars) == 0L) {
    stop("`axes` must not be an empty string.")
  }
  if (!all(chars %in% c("x", "y", "z"))) {
    stop("`axes` must only contain the characters 'x', 'y', and/or 'z'.")
  }
  if (anyDuplicated(chars)) {
    stop("`axes` must not contain duplicate axis labels.")
  }
  chars
}


#' Infer axes from the number of columns / vector length
#'
#' @param n Integer. Length of vector or number of columns.
#' @return A character vector of axis labels.
#' @noRd

.infer_axes <- function(n) {
  switch(as.character(n),
         "3" = c("x", "y", "z"),
         "2" = c("x", "y"),
         "1" = "x",
         stop(
           "Cannot auto-infer `axes`: length/columns (", n, ") is not 1, 2, or 3. ",
           "Please supply `axes` explicitly."
         )
  )
}


#' Infer axes from vector length with divisibility fallback
#'
#' Tries 3, then 2, then 1 column(s).
#'
#' @param len Integer. Length of the numeric vector.
#' @return A character vector of axis labels.
#' @noRd

.infer_axes_from_length <- function(len) {
  if (len %% 3L == 0L) return(c("x", "y", "z"))
  if (len %% 2L == 0L) return(c("x", "y"))
  "x"
}

#' Convert to integer, but stricter
#'
#' Convert to integer checking against tolerance.
#'
#' @param x a numeric vector
#' @return an integer vector
#' @noRd
.as_integer_strict <- function(x, tol = sqrt(.Machine$double.eps)) {
  r <- round(x)
  if (abs(x - r) > tol) {
    stop(deparse(substitute(x)), " (", x, ") is not an integer value.",
         call. = FALSE)
  }
  as.integer(r)
}

# Null-coalescing operator (base R < 4.4 compatibility)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Calibration ----
#' Annotation to label
#'
#' Converts various possible way to annotate acceleration region
#' to a common character string.
#'
#' @param label A character string.
#' @return A character string.
#' @noRd
.canon_face <- function(label) {
  label <- tolower(trimws(label))
  lookup <- c(
    "+x" = "+x", "px" = "+x", "posx" = "+x", "upx"  = "+x",
    "-x" = "-x", "nx" = "-x", "negx" = "-x", "downx"= "-x",
    "+y" = "+y", "py" = "+y", "posy" = "+y", "upy"  = "+y",
    "-y" = "-y", "ny" = "-y", "negy" = "-y", "downy"= "-y",
    "+z" = "+z", "pz" = "+z", "posz" = "+z", "upz"  = "+z",
    "-z" = "-z", "nz" = "-z", "negz" = "-z", "downz"= "-z",

    "x+" = "+x", "xp" = "+x", "xpos" = "+x", "xup"  = "+x",
    "x-" = "-x", "xn" = "-x", "xneg" = "-x", "xdown"= "-x",
    "y+" = "+y", "yp" = "+y", "ypos" = "+y", "yup"  = "+y",
    "y-" = "-y", "yn" = "-y", "yneg" = "-y", "ydown"= "-y",
    "z+" = "+z", "zp" = "+z", "zpos" = "+z", "zup"  = "+z",
    "z-" = "-z", "zn" = "-z", "zneg" = "-z", "zdown"= "-z"
  )
  unname(lookup[label])  # returns NA if not found
}

#' Trims column means
#'
#' Trims column means (10% each side) as a robust estimator
#' to a common character string.
#'
#' @param m A matrix.
#' @return trimmed column means (10% each side)
#' @noRd
.robust_mean <- function(m) {
  apply(m, 2L, mean, trim = 0.1)
}

#' List to vector of means
#'
#' Average a list of mean vectors into a single vector
#'
#' @param lst A list.
#' @return A vector of means
#' @noRd
.avg_face <- function(lst) {
  colMeans(do.call(rbind, lst))
}

#' Safe asin
#'
#' Clamp scalar to [-1, 1] before asin to avoid NaN due to floating point errors
#'
#' @param x A scalar.
#' @return asin
#' @noRd
.safe_asin <- function(x, tol = sqrt(.Machine$double.eps)) {
  if (abs(x) > 1 + tol) {
    stop("asin argument (", x, ") is outside [-1, 1] by more than floating ",
         "point tolerance, which may indicate a problem with the input data ",
         "or calibration procedure.", call. = FALSE)
  }
  asin(max(-1, min(1, x)))
}
