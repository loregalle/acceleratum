#' Process annotations into a segment list for calibration
#'
#' Extracts data segments from an accelerometry object according to
#' \code{annotations}, computes a robust mean per segment, and resolves the
#' target gravity vector for each segment following this priority:
#' \enumerate{
#'   \item If \code{label} resolves to a canonical face (via
#'     \code{.canon_face()}), \code{gvec} is derived from that canonical
#'     axis (scaled to magnitude \code{g}). Any user-supplied \code{gvec}
#'     for that row is ignored, with a message.
#'   \item Else, if the row supplies a \code{gvec} directly, it is
#'     normalised to magnitude \code{g} and used as-is.
#'   \item Else, \code{gvec} is left \code{NULL}, to be estimated from the
#'     data itself during calibration.
#' }
#'
#' @param acc An \code{aclrtm_accelerometry} object with axes x, y, and z.
#' @param annotations A data frame with at least columns \code{xmin} and
#'   \code{xmax}. May optionally include a \code{label} column (character)
#'   and/or a \code{gvec} list-column (each entry a length-3 numeric
#'   vector, or \code{NULL}).
#' @param g Gravitational acceleration magnitude used to scale resolved
#'   gvecs. Default 1.
#'
#' @returns A list, one element per valid annotation row, each a list with:
#'   \describe{
#'     \item{data}{The raw segment data (rows of \code{acc} within range).}
#'     \item{mean}{10%-trimmed mean of \code{data}, per axis.}
#'     \item{glabel}{Canonical face string (\code{"+x"}, \code{"-z"}, etc.)
#'       if resolvable, else \code{NA}.}
#'     \item{gvec}{Resolved target gravity vector (magnitude \code{g}), or
#'       \code{NULL} if it must be estimated from data during calibration.}
#'   }
#'
#' @export
process_annotations <- function(acc, annotations, g = 1) {

  axes <- c("x", "y", "z")
  missing <- setdiff(axes, colnames(acc))
  if (length(missing) != 0) {
    stop("process_annotations requires triaxial data. Columns x, y, and z ",
         "must be present in `acc`. Missing: ",
         paste0(missing, collapse = ", "),
         ".", call. = FALSE)
  }
  if (!all(c("xmin", "xmax") %in% names(annotations))) {
    stop("`annotations` must have at least columns xmin and xmax.",
         call. = FALSE)
  }

  has_label <- "label" %in% names(annotations)
  has_gvec  <- "gvec"  %in% names(annotations)

  sr <- attr(acc, "sampling_rate")
  st <- attr(acc, "start_time")
  time_stamp <- if (is.null(sr)) {
    seq_len(nrow(acc))
  } else {
    seq(0, by = 1 / sr, length.out = nrow(acc))
  }

  canon_axis <- list(
    "+x" = c(1, 0, 0), "-x" = c(-1, 0, 0),
    "+y" = c(0, 1, 0), "-y" = c(0, -1, 0),
    "+z" = c(0, 0, 1), "-z" = c(0, 0, -1)
  )

  n   <- nrow(annotations)
  out <- vector("list", n)

  message("Processing ", n, " annotations.")

  for (i in seq_len(n)) {
    r <- annotations[i, ]

    if (!is.numeric(r$xmin) || !is.numeric(r$xmax) || r$xmin > r$xmax) {
      warning("Invalid xmin/xmax for annotation ", i, ". Skipping.",
              call. = FALSE)
      next
    }

    idx <- which(time_stamp >= r$xmin & time_stamp <= r$xmax)
    if (length(idx) == 0L) {
      warning("No data found for annotation ", i, ". Skipping.",
              call. = FALSE)
      next
    }

    seg <- acc[idx, axes, drop = FALSE]
    m   <- apply(seg, 2L, mean, trim = 0.1)

    glabel <- NULL
    if (has_label && !is.na(r$label)) {
      canon <- .canon_face(r$label)
      if (!is.na(canon)) glabel <- canon
    }

    user_gvec <- if (has_gvec) annotations$gvec[[i]] else NULL
    gvec <- NULL

    if (!is.null(glabel)) {
      gvec <- canon_axis[[glabel]] * g
      if (!is.null(user_gvec)) {
        message("Annotation ", i, ": canonical label \"", glabel,
                "\" present; supplied gvec ignored, using label-derived gvec.")
      }
    } else if (!is.null(user_gvec)) {
      nrm <- sqrt(sum(user_gvec^2))
      if (nrm < 1e-12) {
        warning("Annotation ", i, ": supplied gvec has ~zero norm; ignoring.",
                call. = FALSE)
      } else {
        gvec <- user_gvec / nrm * g
      }
    }
    # else: gvec stays NULL -> resolved from data during calibration

    out[[i]] <- list(data = seg, mean = m, glabel = glabel, gvec = gvec)
  }

  out[!vapply(out, is.null, logical(1))]
}


#' Build the least-squares design matrix and observation vector
#'
#' Generalises the diagonal-only assembler to optionally include
#' off-diagonal (misalignment) terms and/or a bias term.
#'
#' @param rows A list of \code{list(gvec = c(gx,gy,gz), m = c(ax,ay,az))}.
#' @param diagonal_only If \code{TRUE}, fit only \code{diag(sx,sy,sz)}.
#'   If \code{FALSE}, fit a full 3x3 matrix.
#' @param estimate_bias If \code{TRUE}, include a bias term.
#'
#' @returns A list with \code{A} (design matrix) and \code{y} (observations).
#' @noRd
.assemble_ls <- function(rows, diagonal_only, estimate_bias) {
  n <- length(rows)

  ncol_A <- if (diagonal_only) 3L else 9L

  if (estimate_bias) {
    ncol_A <- ncol_A + 3
    idx_bias <- ncol_A - (2:0)
  }

  A <- matrix(0, nrow = n * 3L, ncol = ncol_A)
  y <- numeric(n * 3L)

  for (i in seq_len(n)) {
    row <- rows[[i]]
    ridx <- (i - 1L) * 3L

    if (diagonal_only) {
      A[ridx + 1:3, 1:3] <- diag(row$gvec)
    } else {
      A[ridx + 1:3, 1:9] <- diag(3) %x% t(row$gvec)
    }

    if (estimate_bias) {
      A[ridx + 1:3, idx_bias] <- diag(3)
    }

    y[ridx + 1:3] <- row$m
  }

  list(A = A, y = y)
}



#' Generic iterative affine calibration (diagonal or full misalignment)
#'
#' @param seg_list Output of \code{process_annotations()} (or equivalent
#'   list with \code{$mean} and optionally \code{$gvec} per element).
#' @param g Gravitational acceleration magnitude.
#' @param diagonal_only Logical. Fit \code{M = diag(sx,sy,sz)} only, or a
#'   full 3x3 misalignment matrix.
#' @param estimate_bias Logical. Whether to fit a bias vector.
#' @param max_iter Maximum number of iterations.
#' @param tol Convergence tolerance on max absolute parameter change.
#' @param verbose If true, prints information for each iteration.
#'
#' @returns A named list with \code{M}, \code{Minv}, \code{bias},
#'   \code{seg_means}, \code{iters}, \code{dP_hist}.
#' @noRd
.cal_affine <- function(seg_list, g, diagonal_only, estimate_bias,
                        max_iter, tol, verbose) {

  n <- length(seg_list)
  if (n < 3L) {
    stop("Insufficient segments for calibration (need at least 3, got ",
         n, ").", call. = FALSE)
  }

  M_combined <- diag(3)
  b <- rep(0, 3)

  fixed_gvec <- lapply(seg_list, function(s) s$gvec)  # NULL where unknown

  it_done <- 0L

  rows <- vector("list", n)
  for (i in seq_len(n)) {
    rows[[i]] <- list(gvec = fixed_gvec[[i]], m = seg_list[[i]]$mean)
  }
  free_idx <- which(vapply(fixed_gvec, is.null, logical(1)))

  for (it in seq_len(max_iter)) {
    for (i in free_idx) {
      mcorr <- solve(M_combined, rows[[i]]$m - b)
      nrm   <- sqrt(sum(mcorr^2)) + 1e-12
      rows[[i]]$gvec <- mcorr / nrm * g
    }

    asmb <- .assemble_ls(rows, diagonal_only, estimate_bias)
    fit  <- qr.solve(asmb$A, asmb$y)

    if (diagonal_only) {
      scales_new <- fit[1:3]
      b_new      <- if (estimate_bias) fit[4:6] else rep(0, 3)
      M_new      <- diag(scales_new)
    } else {
      M_new <- matrix(fit[1:9], nrow = 3L, byrow = TRUE)
      b_new <- if (estimate_bias) fit[10:12] else rep(0, 3)
    }

    dP      <- max(abs(c(as.vector(M_new - M_combined), b_new - b)))
    it_done <- it

    if (verbose) {
      message(sprintf("  [Affine] iter %4d: dP = %.4e", it, dP))
    }

    M_combined <- M_new
    b <- b_new
    if (dP < tol) break
  }

  message("  [Affine] converged in ", it_done, " iterations")

  S <- diag(apply(M_combined, 1, \(.x) sqrt(sum(.x^2))))
  M <- solve(S) %*% M_combined


  list(
    S         = S,
    M         = M,
    bias      = b,
    iters     = it_done,
    deltaP    = dP,
    beta_deg  = NA
  )
}


#' Xu et al. (2021) six-position closed-form calibration branch
#'
#' Requires all six canonical faces (\code{+x,-x,+y,-y,+z,-z}) to be
#' present via \code{$glabel} in \code{seg_list}. Estimates scale factors,
#' full misalignment, bias, and platform tilt \code{beta}.
#'
#' @inheritParams .cal_affine
#' @returns A named list matching the original \code{cal_xu()} output shape.
#' @noRd
.cal_xu <- function(seg_list, g, max_iter, tol, verbose) {

  required_faces <- c("+x", "-x", "+y", "-y", "+z", "-z")

  face_means <- list()
  for (s in seg_list) {
    if (!is.null(s$glabel) && s$glabel %in% required_faces) {
      face_means[[s$glabel]] <- c(face_means[[s$glabel]], list(s$mean))
    }
  }

  missing_faces <- setdiff(required_faces, names(face_means))
  if (length(missing_faces) > 0L) {
    stop(
      "estimate_tilt = TRUE requires all six canonical faces. Missing: ",
      paste(missing_faces, collapse = ", "), "\n",
      "Found faces: ", paste(sort(names(face_means)), collapse = ", "),
      call. = FALSE
    )
  }

  Vxu <- .avg_face(face_means[["+x"]])
  Vxd <- .avg_face(face_means[["-x"]])
  Vyu <- .avg_face(face_means[["+y"]])
  Vyd <- .avg_face(face_means[["-y"]])
  Vzu <- .avg_face(face_means[["+z"]])
  Vzd <- .avg_face(face_means[["-z"]])

  Vx1 <- (Vxu - Vxd) / 2;  Vx2 <- (Vxu + Vxd) / 2
  Vy1 <- (Vyu - Vyd) / 2;  Vy2 <- (Vyu + Vyd) / 2
  Vz1 <- (Vzu - Vzd) / 2;  Vz2 <- (Vzu + Vzd) / 2

  gL <- as.double(g)

  beta <- t_xy <- t_xz <- t_yx <- t_yz <- t_zx <- t_zy <- 0
  SFx  <- SFy  <- SFz  <- 1
  bx   <- by   <- bz   <- 0

  X_prev  <- rep(0, 13)
  iters   <- 0L

  for (it in seq_len(max_iter)) {
    iters <- it
    cb <- cos(beta)
    sb <- sin(beta)

    denom_x <- gL * cb * cos(t_xy) * cos(t_xz)
    denom_y <- gL * cb * cos(t_yx) * cos(t_yz)
    denom_z <- gL * cb * cos(t_zx) * cos(t_zy)

    if (abs(denom_x) < 1e-15 || abs(denom_y) < 1e-15 || abs(denom_z) < 1e-15) {
      warning("Near-zero denominator at iteration ", it, ". Aborting.",
              call. = FALSE)
      break
    }

    SFx <- Vx1[1L] / denom_x
    SFy <- Vy1[2L] / denom_y
    SFz <- Vz1[3L] / denom_z

    bx <- Vx2[1L] - gL * SFx * sb * sin(t_xz)
    by <- Vy2[2L] - gL * SFy * sb * cos(t_yx) * sin(t_yz)
    bz <- Vz2[3L] + gL * SFz * sb * cos(t_zy) * sin(t_zx)

    t_yx <- asin(Vz1[2L] / (gL * SFy * cb))
    t_zy <- asin(Vx1[3L] / (gL * SFz * cb))
    t_xz <- asin(Vy1[1L] / (gL * SFx * cb))

    d_xy <- -gL * SFx * cb * cos(t_xz)
    d_yz <- -gL * SFy * cb * cos(t_yx)
    d_zx <- -gL * SFz * cb * cos(t_zy)

    if (abs(d_xy) > 1e-15) t_xy <- asin(Vz1[1L] / d_xy)
    if (abs(d_yz) > 1e-15) t_yz <- asin(Vx1[2L] / d_yz)
    if (abs(d_zx) > 1e-15) t_zx <- asin(Vy1[3L] / d_zx)

    d_sv1 <-  gL * SFy * cos(t_yx) * cos(t_yz)
    d_sv2 <- -gL * SFx * cos(t_xy) * cos(t_xz)

    sv <- 0
    if (abs(d_sv1) > 1e-15) sv <- sv + (Vx2[2L] - by) / d_sv1
    if (abs(d_sv2) > 1e-15) sv <- sv + (Vy2[1L] - bx) / d_sv2
    if (abs(d_sv1) > 1e-15) sv <- sv + (Vz2[2L] - by) / d_sv1

    beta <- asin(sv / 3)

    X_curr <- c(SFx, SFy, SFz, bx, by, bz,
                t_xy, t_xz, t_yx, t_yz, t_zx, t_zy, beta)

    if (it > 1L) {
      e0      <- max(abs(X_curr - X_prev))
      if (verbose) {
        message(sprintf("  [Xu] iter %4d: deltaP = %.4e", it, e0))
      }
      if (e0 < tol) break
    }

    X_prev <- X_curr
  }

  message("  [Xu] converged in ", iters, " iterations")

  A_sb <- matrix(c(
    cos(t_xy) * cos(t_xz),  sin(t_xz), -cos(t_xz) * sin(t_xy),
    -cos(t_yx) * sin(t_yz),  cos(t_yx) * cos(t_yz),  sin(t_yx),
    sin(t_zy), -cos(t_zy) * sin(t_zx),  cos(t_zx) * cos(t_zy)
  ), nrow = 3L, byrow = TRUE)

  S_xu <- diag(c(SFx, SFy, SFz)) %*% A_sb
  M    <- A_sb
  S    <- c(SFx, SFy, SFz)
  bias <- c(bx, by, bz)

  .deg <- function(r) r * 180 / pi

  list(
    S         = S,               # scaling factor matrix
    M         = M,               # misalignment matrix
    bias      = bias,            # bias vector
    iters     = iters,
    deltaP    = if (length(e0) > 0) e0 else 0,
    beta_deg  = .deg(beta)
  )
}


#' Unified affine accelerometer calibration
#'
#' Consolidates diagonal, full-misalignment, and Xu et al. (2021)
#' six-position tilt calibration into a single entry point.
#'
#' @param seg_list A list as returned by \code{process_annotations()}: each
#'   element a list with \code{$data}, \code{$mean}, optionally
#'   \code{$glabel} and/or \code{$gvec}.
#' @param g Gravitational acceleration magnitude.
#' @param diagonal_only Logical, default \code{FALSE}. If \code{TRUE}, fit
#'   \code{M = diag(sx,sy,sz)} only (no cross-axis/misalignment terms).
#'   Ignored if \code{estimate_tilt = TRUE}.
#' @param estimate_bias Logical, default \code{TRUE}. Whether to fit a bias
#'   vector. Ignored if \code{estimate_tilt = TRUE}.
#' @param estimate_tilt Logical, default \code{FALSE}. If \code{TRUE},
#'   runs the Xu et al. (2021) six-position closed-form calibration.
#' @param max_iter Maximum number of iterations.
#' @param tol Convergence tolerance.
#' @param verbose If true, prints information for each iteration.
#'
#' @returns A named list. For the generic branches: \code{M}, \code{Minv},
#'   \code{bias}, \code{seg_means}, \code{iters}, \code{dP_hist}. For the
#'   tilt branch: \code{S}, \code{M}, \code{bias}, \code{seg_means},
#'   \code{iters}, \code{e0_hist}, \code{params_xu}.
#'
#' @references Xu, T., Xu, X., Xu, D., Zhao, H., 2021. A Novel Calibration
#' Method Using Six Positions for MEMS Triaxial Accelerometer. IEEE
#' Transactions on Instrumentation and Measurement 70, 1-11.
#' https://doi.org/10.1109/TIM.2020.3026024
#'
#' @export
cal_accel <- function(seg_list,
                      g = 1,
                      diagonal_only = FALSE,
                      estimate_bias = TRUE,
                      estimate_tilt = FALSE,
                      max_iter = 500L,
                      tol = 1e-8,
                      verbose = TRUE) {

  if (estimate_tilt) {
    if (diagonal_only || !estimate_bias) {
      warning(
        "estimate_tilt = TRUE runs the Xu six-position calibration, ",
        "which always estimates scale, full misalignment, bias, and tilt ",
        "jointly. `diagonal_only` and `estimate_bias` are ignored.",
        call. = FALSE
      )
    }
    return(.cal_xu(seg_list, g = g, max_iter = max_iter,
                   tol = tol, verbose = verbose))
  }

  .cal_affine(seg_list, g = g,
              diagonal_only = diagonal_only,
              estimate_bias = estimate_bias,
              max_iter = max_iter, tol = tol,
              verbose = verbose)
}
