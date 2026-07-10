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

#' Greatest common divisor of two integers
#'
#' Computes the GCD of \code{a} and \code{b} using the Euclidean algorithm.
#'
#' @param a A non-negative integer value.
#' @param b A non-negative integer value.
#'
#' @returns A length-1 integer giving the greatest common divisor of \code{a}
#'   and \code{b}. Returns \code{a} when \code{b} is zero.
#'
#' @noRd
.gcd2 <- function(a, b) {
  if (b == 0) a else Recall(b, a %% b)
}

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
    v1 - rowSums(v1 * v2)/rowSums(v2 * v2) * v2
  } else {
    v1 - sum(v1 * v2)/sum(v2 * v2) * v2
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

# Screening and selecting ----

#' Parse a single burst cell
#'
#' Parse a burst cell (space- or sep-separated numeric string) into a
#' samples x n_axes matrix, row-major interleaved convention
#' (x1 y1 z1 x2 y2 z2 ...)
#'
#' @param s A string of sep-separated values
#' @param n_axes Number of axes in the data
#' @param sep Character separator
#' @noRd
.parse_burst_string <- function(s, n_axes, sep = " ") {
  s <- trimws(s)
  if (!nzchar(s)) return(matrix(numeric(0), ncol = n_axes))
  if (identical(sep, " ")) {
    parts <- strsplit(s, "\\s+")[[1]]
  } else {
    parts <- strsplit(s, sep, fixed = TRUE)[[1]]
  }
  vals <- suppressWarnings(as.numeric(parts))
  if (anyNA(vals)) {
    stop("Non-numeric value encountered while parsing a burst cell.")
  }
  if (length(vals) %% n_axes != 0) {
    stop(sprintf(
      "Burst cell length (%d) is not a multiple of n_axes (%d).",
      length(vals), n_axes
    ))
  }
  matrix(vals, ncol = n_axes, byrow = TRUE)
}

#' Record size in bytes
#'
#' Fixed record size (bytes) for one reservoir slot in the scratch file.
#'
#' @param slot_size self_explanatory
#' @param n_axes number of axes in data matrix
#'
#' @details In order (int32 = 4 bytes, float64 = 8 bytes):
#' \itemize{
#'  \item{actual_n : 1 x int32}
#'  \item{timestamps : slot_size x float64}
#'  \item{raw samples : slot_size x n_axes x float64}
#'  \item{source burst-row id per sample : slot_size x int32}
#' }
#' @noRd
.record_size_bytes <- function(slot_size, n_axes) {
  4L + slot_size * 8L + slot_size * n_axes * 8L + slot_size * 4L
}

#' Write slot in reservoir
#'
#' Write one (possibly replacing) reservoir slot to the scratch file.
#'
#' @param state The state environment
#' @param slot Slot in the reservoir
#' @param ts_w Timestamp values
#' @param raw_w Matrix of data values
#' @param row_w The sample burst-row ID
#' @noRd
.write_slot <- function(state, slot, ts_w, raw_w, row_w) {
  n <- length(ts_w)
  if (n > state$slot_size) {
    stop(sprintf(
      paste0(
        "A stationary window required %d samples but `slot_size` is only %d. ",
        "Re-run mine_reservoir() with a larger slot_size (e.g. slot_size = %d) ",
        "or adjust window_sec / sr_tol."
      ),
      n, state$slot_size, n
    ))
  }

  # size of padding needed to fill slot size
  pad_n   <- state$slot_size - n
  # pad timestamp
  ts_pad  <- c(ts_w, numeric(pad_n))
  # pad raw value matrix
  raw_pad <- rbind(raw_w, matrix(0, nrow = pad_n, ncol = ncol(raw_w)))
  # pad row ID
  row_pad <- c(row_w, integer(pad_n))

  # Offset to slot initial position
  offset <- (slot - 1) * state$record_size_bytes
  # reposition connection in scratch file
  seek(state$scratch_con, where = offset, origin = "start")
  # write data to scratch file
  writeBin(as.integer(n), state$scratch_con, size = 4)
  writeBin(as.double(ts_pad), state$scratch_con, size = 8)
  writeBin(as.double(t(raw_pad)), state$scratch_con, size = 8)
  writeBin(as.integer(row_pad), state$scratch_con, size = 4)
}

#' Pairwise angular distance
#'
#' Compute the pairwise angular distance
#' between the rows of two unit-vector matrices. Used both for
#' vector-vs-single-vector distance (Farthest Point Sampling candidate
#' selection) and for full pairwise distance within a selected set
#' (mean pairwise angular distance scoring).
#'
#' @param A A numeric matrix of unit row vectors (n x d).
#' @param B A numeric matrix of unit row vectors (m x d).
#'
#' @return An n x m matrix of angular distances (radians), where entry
#'   \code{[i, j]} is the angular distance between \code{A[i, ]} and
#'   \code{B[j, ]}.
#'
#' @details Rows of \code{A} and \code{B} are assumed to already be unit
#'   vectors; this function does not
#'   normalise them.
#' @noRd
.angular_distance <- function(A, B) {
  dots <- A %*% t(B)
  dots <- pmin(pmax(dots, -1), 1)
  acos(dots)
}

#' Evaluate the current window buffer and admit it to the reservoir if stationary
#'
#' For internal use only!
#' Test the sliding-window buffer's current live contents (\code{state$buf_start}
#' to \code{state$buf_end}) for stationarity using three criteria -- VeDBA
#' (mean per-sample L2 norm of the dynamic acceleration), mean per-axis
#' variance, and minimum window mean-vector magnitude -- all of which must
#' pass for the window to be admitted. On a pass, runs one step of Vitter's
#' reservoir sampling algorithm to determine the target slot \code{j} (the
#' next empty slot while the reservoir is filling, or a randomly drawn slot
#' once full, discarded if it falls outside the reservoir), then stores the
#' window's mean vector at \code{j} and persists its raw samples to the
#' scratch file via a single \code{.write_slot()} call. Assumes the caller
#' has already established that the buffer spans a full window (see
#' \code{.buf_append_and_evict()}); does not itself check window length or
#' sample count.
#'
#' @param state The shared mutable state environment (see \code{mine_reservoir()}).
#'   Updated in place: increments \code{windows_total} on every call, and on
#'   a pass increments \code{passing_seen} and updates \code{reservoir_means}
#'   / \code{reservoir_filled} accordingly.
#'
#' @noRd
.evaluate_and_admit_window <- function(state) {
  idx   <- state$buf_start:state$buf_end
  ts_w  <- state$buf_ts[idx]
  raw_w <- state$buf_raw[idx, , drop = FALSE]
  row_w <- state$buf_row[idx]

  mean_vec <- colMeans(raw_w)
  dyn        <- sweep(raw_w, 2, mean_vec, "-")
  vedba_mean <- mean(sqrt(rowSums(dyn^2))) # VeDBA: L2 norm per sample
  var_mean   <- mean(apply(raw_w, 2, stats::var))
  mag_mean   <- sqrt(sum(mean_vec^2))

  state$windows_total <- state$windows_total + 1L

  pass_now <- (vedba_mean < state$vedba_thresh) &&
    (var_mean   < state$var_thresh)   &&
    (mag_mean  >= state$mag_thresh)

  if (!pass_now) return(invisible(NULL))

  state$passing_seen <- state$passing_seen + 1L

  if (state$reservoir_filled < state$reservoir_size) {
    j <- state$reservoir_filled + 1L
    state$reservoir_filled <- j
  } else {
    j <- sample.int(state$passing_seen, 1L)
    if (j > state$reservoir_size) return(invisible(NULL))
  }

  state$reservoir_means[j, ] <- mean_vec
  .write_slot(state, j, ts_w, raw_w, row_w)
  invisible(NULL)
}

#' Append a sample to the sliding-window buffer and evict expired samples
#'
#' For internal use only!
#' Append one reconstructed sample (timestamp, raw vector, source row id) to
#' the growable sliding-window buffer stored on \code{state}, growing or
#' compacting the underlying storage as needed, then evict samples from the
#' front of the buffer whose timestamp is more than \code{state$window_sec}
#' behind the newly appended sample. Does not evaluate the resulting window
#' for stationarity or admit it to the reservoir; call
#' \code{.evaluate_and_admit_window()} separately once the caller determines
#' the buffer spans a full window.
#'
#' @param state The shared mutable state environment (see \code{mine_reservoir()}).
#' @param ts Timestamp of the sample being appended.
#' @param raw_row A single sample's raw values, as a numeric vector of length
#'   \code{n_axes}.
#' @param row_id Source burst-row identifier for the sample.
#'
#' @noRd
.buf_append_and_evict <- function(state, ts, raw_row, row_id) {
  if (state$buf_end >= state$buf_cap) {
    live_len <- state$buf_end - state$buf_start + 1L
    if (state$buf_start > 1L) {
      idx <- state$buf_start:state$buf_end
      state$buf_ts[seq_len(live_len)]      <- state$buf_ts[idx]
      state$buf_raw[seq_len(live_len), ]   <- state$buf_raw[idx, , drop = FALSE]
      state$buf_row[seq_len(live_len)]     <- state$buf_row[idx]
      state$buf_start <- 1L
      state$buf_end   <- live_len
    } else {
      new_cap <- state$buf_cap * 2L
      new_ts  <- numeric(new_cap)
      new_ts[seq_len(state$buf_cap)]  <- state$buf_ts
      new_raw <- matrix(0, nrow = new_cap, ncol = ncol(state$buf_raw))
      new_raw[seq_len(state$buf_cap), ] <- state$buf_raw
      new_row <- integer(new_cap)
      new_row[seq_len(state$buf_cap)] <- state$buf_row
      state$buf_ts <- new_ts
      state$buf_raw <- new_raw
      state$buf_row <- new_row
      state$buf_cap <- new_cap
    }
  }

  state$buf_end <- state$buf_end + 1L
  state$buf_ts[state$buf_end]    <- ts
  state$buf_raw[state$buf_end, ] <- raw_row
  state$buf_row[state$buf_end]   <- row_id

  while (
    (state$buf_end > state$buf_start) &&
    (state$buf_ts[state$buf_end] - state$buf_ts[state$buf_start] > state$window_sec)
  ) {
    state$buf_start <- state$buf_start + 1L
  }
}

#' Resolve a pending row and feed its samples into the sliding window
#'
#' For internal use only!
#' Derive (or reuse) a per-sample sampling rate for the pending row
#' \code{P}, reconstruct its per-sample timestamps, and push each
#' resulting sample into the sliding-window buffer. If the derived rate is
#' implausible (deviates from the nominal rate by more than
#' \code{state$sr_tol}), the row's samples are discarded and the window
#' buffer is cleared, mirroring a gap-detection reset.
#'
#' @param state The shared mutable state environment (see \code{mine_reservoir()}).
#' @param P A pending row, as a list with elements \code{ts} (row start
#'   timestamp), \code{mat} (samples x n_axes matrix), and \code{row_id}
#'   (source burst-row identifier).
#' @param ts_C Timestamp of the row immediately following \code{P}, used to
#'   derive \code{P}'s per-sample rate as \code{nrow(P$mat) / (ts_C - P$ts)}.
#'   Ignored (and may be omitted) when \code{sr_override = TRUE}.
#' @param sr_override If \code{TRUE}, skip rate derivation and plausibility
#'   checking entirely and use \code{state$nominal_sr} directly. Intended for
#'   the file's final row, which has no successor row to derive a rate from.
#'
#' @noRd
.resolve_and_process_row <- function(state, P,
                                     ts_C = NULL, sr_override = FALSE) {
  n_P <- nrow(P$mat)
  if (n_P == 0L) return(invisible(NULL))  # nothing to contribute either way

  if (sr_override) {
    sr_P <- state$nominal_sr
  } else {
    sr_P <- n_P / (ts_C - P$ts)
    if (abs(sr_P - state$nominal_sr) > state$sr_tol) {
      # cannot reliably reconstruct timestamps across this row, clear buffer
      state$buf_start <- 1L
      state$buf_end   <- 0L
      return(invisible(NULL))
    }
  }

  ts_samples <- P$ts + (seq_len(n_P) - 1L) / sr_P
  for (j in seq_len(n_P)) {
    .buf_append_and_evict(state, ts_samples[j], P$mat[j, ], P$row_id)

    n_buf <- state$buf_end - state$buf_start + 1L
    if (
      n_buf >= state$min_samples_per_window &&
      (state$buf_ts[state$buf_end] - state$buf_ts[state$buf_start]) >= state$window_sec
    ) {
      .evaluate_and_admit_window(state)
    }
  }
  invisible(NULL)
}
