# Generals ----
#' Parse a axes string into a character vector of axis labels
#'
#' @param axes A string composed of any non-repeating combination of "x", "y",
#'   "z" (e.g. "xyz", "xz", "y").
#' @returns A character vector such as c("x","y","z").
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
#' @returns A character vector of axis labels.
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
#' @returns A character vector of axis labels.
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
#' @returns an integer vector
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
#' @returns A character string.
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
#' @returns trimmed column means (10% each side)
#' @noRd
.robust_mean <- function(m) {
  apply(m, 2L, mean, trim = 0.1)
}

#' List to vector of means
#'
#' Average a list of mean vectors into a single vector
#'
#' @param lst A list.
#' @returns A vector of means
#' @noRd
.avg_face <- function(lst) {
  colMeans(do.call(rbind, lst))
}

#' Safe asin
#'
#' Clamp scalar to \[-1, 1\] before asin to avoid NaN due to floating point errors
#'
#' @param x A scalar.
#' @returns asin
#' @noRd
.safe_asin <- function(x, tol = sqrt(.Machine$double.eps)) {
  if (abs(x) > 1 + tol) {
    stop("asin argument (", x, ") is outside [-1, 1] by more than floating ",
         "point tolerance, which may indicate a problem with the input data ",
         "or calibration procedure.", call. = FALSE)
  }
  asin(max(-1, min(1, x)))
}

# Rotation ----
#' Rodrigues' rotation formula
#'
#' Uses Rodrigues' rotation formula to find the rotation matrix that rotates
#' the first input vector to be aligned with the second input vector.
#'
#' @param from A vector.
#' @param to A vector.
#' @param tol Tolerance value to test for parallel vectors
#' @returns A rotation matrix \code{R} such as \code{R %*% from}
#'   equals \code{to}
#' @noRd
.rodrigues_rotation <- function(from, to, tol = 1e-10) {
  # normalise
  from <- from / sqrt(sum(from^2))
  to   <- to   / sqrt(sum(to^2))

  # cross product: rotation vector
  k <- .cross(from, to)

  # sin theta via cross product
  sin_theta <- sqrt(sum(k^2))

  # cos theta via dot product
  cos_theta <- as.vector(from %*% to)

  if (sin_theta < tol) {        # if parallel
    if (cos_theta > 0) {        # if codirectional
      return(diag(3L))          # R is the identity
    }
    # if not codirectional, 180 degree rotation
    seed_vec <- ifelse(        # use some simple logic to pick any
      abs(from[1L]) < 0.9,     # vector different enough from "from"
      c(1,0,0),
      c(0,1,0)
    )
    k <- .gram_schmidt(seed_vec, from)      # orthogonalisation
    k <- k / sqrt(sum(k^2))                 # normalisation
    return(2 * outer(k, k) - diag(3))       # solving Rodrigues with theta = pi
  }

  # skew-symmetric cross product matrix
  k <- k/sin_theta
  K <- diag(0,3)
  K_offdiag <- k[c(3,2,3,1,2,1)] * c(1,-1,-1,1,1,-1)
  K[row(K) != col(K)] <- K_offdiag

  # Rodrigues formula
  diag(3L) + sin_theta * K + (1-cos_theta) * K %*% K
}

#' Givens' rotation
#'
#' Uses Givens' rotation to find the rotation matrix that rotates
#' the first input vector to be aligned with the second input vector.
#'
#' @param from A vector.
#' @param to A vector.
#' @param rot_axes the axes that define the plane being rotated
#' @param axes the names of all the axes in the rotation matrix to be returned
#' @returns A rotation matrix
#' @noRd
.givens_rotation <- function(from, to, rot_axes, axes = c("x", "y", "z")) {
  from <- from / sqrt(sum(from^2))
  to   <- to   / sqrt(sum(to^2))

  R <- diag(length(axes))
  rownames(R) <- colnames(R) <- axes
  theta <- atan2(to[2L], to[1L]) - atan2(from[2L], from[1L])

  R[rot_axes, rot_axes] <- c(
    cos(theta), sin(theta), -sin(theta), cos(theta)
  )
  R
}

#' Gram-Schmidt orthogonalisation
#'
#' Uses Gram-Schmidt orthogonalisation to return the projection of vector
#' \eqn{\mathbf{v}_1} on the plane that is perpendicular to \eqn{\mathbf{v}_2}
#'
#' @param v1 a 3d vector or matrix
#' @param v2 a 3d vector or matrix
#' @returns a 3d vector
#' @noRd
.gram_schmidt <- function(v1,v2) {
  if (inherits(v1, "matrix") != inherits(v2, "matrix")) {
    stop("'v1' and 'v2' must be either both matrices or both vectors")
  }
  if (inherits(v1, "matrix")) {
    v1 - rowSums(v1 * v2) * v2
  } else {
    v1 - sum(v1 * v2) * v2
  }
}

#' Geometrical cross product
#'
#' Finds the vector that is perpendicular to the plane defined by the
#' \code{a} and \code{b} vectors
#'
#' @param a a 3d vector
#' @param b another 3d vector
#' @noRd
.cross <- function(a, b) {
  out <- c(
    a[c(2L,3L,1L)] * b[c(3L, 1L, 2L)] -
      a[c(3L, 1L, 2L)] * b[c(2L,3L,1L)]
  )
  names(out) <- names(a)
  out
}
