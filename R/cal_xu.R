#' Xu et al. (2021) six-position iterative calibration
#'
#' @param acc An \code{aclrtm_accelerometry} object with axes x, y, and z.
#' @param annotations A data frame as returned by \code{range_select()}.
#'   See details.
#' @param g Gravitational acceleration.
#' @param max_iter Maximum number of iterations.
#' @param tol Convergence tolerance.
#'
#' @details
#' The \code{annotations} data frame should include at least three columns: imin,
#' imax, and label. For each row, imin and imax define the row index range of
#' \code{acc}, while label must describe the face and the direction of the
#' axis being sampled. Each label must be a string formatted as \code{"ab"} or
#' \code{"ba"}, where \code{"a"} represents direction and can be any of \{"p",
#' "pos", "+", "up"\} for positive or \{"n", "neg", "-", "down"\} for negative,
#' and \code{"b"} is one of the axes \{"x", "y", "z"\}.
#' Examples of valid labels: \{"posx", "ypos", "x-", "+z", "nx"\}.
#'
#' @returns A named list with elements:
#' \describe{
#'   \item{S}{Vector of scaling factors}
#'   \item{M}{Misalignment matrix}
#'   \item{bias}{Vector of biases}
#'   \item{seg_means}{List of segment means per annotation}
#'   \item{iters}{Number of iterations performed}
#'   \item{e0_hist}{Numeric vector of max parameter change per iteration}
#'   \item{params_xu}{List of estimated parameters as named in
#'     Xu et al. (2021)}
#' }
#' @export
cal_xu <- function(acc, annotations, g = 1,
                   max_iter = 200L, tol = 1e-6) {

  # input checks
  axes <- c("x", "y", "z")
  missing <- setdiff(axes, colnames(acc))

  if (length(missing) != 0) {
    stop("cal_xu requires triaxial data. Columns x, y, and z must be present",
         " in `acc`. Missing: ", paste0(missing, collapse = ", "), ".",
         .call = FALSE)
  }
  if (!all(c("imin", "imax", "label") %in% names(annotations))) {
    stop("`annotations` must have columns imin, imax, and label. See ?cal_xu",
         " for details.",
         call. = FALSE)
  }

  # reconstruct time axis in plot units
  sr <- attr(acc, "sampling_rate")
  st <- attr(acc, "start_time")

  # 1. Collect segment means
  face_means <- list()

  oj_labels         <- annotations$label
  annotations$label <- .canon_face(annotations$label)
  if (any(is.na(annotations$label))) {
    idl <- which(is.na(annotations$label))
    warning(
      "Labels ", paste(unique(oj_labels[idl]), collapse = ", "),
      " could not be mapped to a canonical face and will be ignored.",
      call. = FALSE
    )
  }

  required_faces <- as.vector(
    outer(c("+", "-"),
          axes,
          paste0,
          "")
  )

  annotations <- annotations[annotations$label %in% required_faces,]
  missing_faces  <- setdiff(required_faces, annotations$label)
  if (length(missing_faces) > 0L) {
    stop(
      "Missing faces: ", paste(missing_faces, collapse = ", "), "\n",
      "Found faces: ", paste(sort(unique(annotations$label)), collapse = ", "),
      call. = FALSE
    )
  }

  n_ann      <- nrow(annotations)
  seg_means  <- vector("list", nrow(annotations))


  message("Collecting segment means from ", n_ann, " valid annotations.")

  for (i in seq_len(n_ann)) {
    r   <- annotations[i, ]

    if (!is.numeric(r$imin) || !is.numeric(r$imax)) {
      warning("At least one of indices imin and imax for annotation ",
      i, " (\"", r$label, "\") is not valid. Skipping.", call. = FALSE)
      next
    }

    idx_limits <- as.integer(round(c(r$imin, r$imax)))
    if (any(abs(idx_limits - c(r$imin, r$imax)) > sqrt(.Machine$double.eps))) {
      warning("At least one of indices imin and imax for annotation ", i,
              " (\"", r$label, "\") was not ",
              "an exact integers and has been rounded.", call. = FALSE)
    }
    if (idx_limits[1] > idx_limits[2]) {
      warning("Index imin is larger than imax for annotation ",
              i, " (\"", r$label, "\"). Skipping.", call. = FALSE)
      next
    }

    idx <- idx_limits[1]:idx_limits[2]

    seg <- acc[idx, axes, drop = FALSE]
    m   <- apply(seg, 2L, mean, trim = 0.1)

    seg_means[[i]] <- list(imin  = r$imin,
                           imax  = r$imax,
                           label = r$label,
                           mean  = m)

    canon <- r$label
    face_means[[canon]] <- c(face_means[[canon]], list(m))
  }

  # 2. Half-differences and half-sums

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

  # 3. Initialise parameters (Step 1)
  beta <- t_xy <- t_xz <- t_yx <- t_yz <- t_zx <- t_zy <- 0
  SFx  <- SFy  <- SFz  <- 1
  bx   <- by   <- bz   <- 0

  X_prev  <- rep(0, 13)
  e0_hist <- double(0)
  iters   <- 0L

  # 4. Iterative loop
  for (it in seq_len(max_iter)) {
    iters <- it
    cb <- cos(beta)
    sb <- sin(beta)

    # Step 1: scale factors (eqs. 19-21)
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

    # Step 2: biases
    bx <- Vx2[1L] - gL * SFx * sb * sin(t_xz)
    by <- Vy2[2L] - gL * SFy * sb * cos(t_yx) * sin(t_yz)
    bz <- Vz2[3L] + gL * SFz * sb * cos(t_zy) * sin(t_zx)

    # Step 3: angles theta_yx, theta_zy, theta_xz
    t_yx <- asin(Vz1[2L] / (gL * SFy * cb))
    t_zy <- asin(Vx1[3L] / (gL * SFz * cb))
    t_xz <- asin(Vy1[1L] / (gL * SFx * cb))

    # Step 4: angles theta_xy, theta_yz, theta_zx
    d_xy <- -gL * SFx * cb * cos(t_xz)
    d_yz <- -gL * SFy * cb * cos(t_yx)
    d_zx <- -gL * SFz * cb * cos(t_zy)

    if (abs(d_xy) > 1e-15) t_xy <- asin(Vz1[1L] / d_xy)
    if (abs(d_yz) > 1e-15) t_yz <- asin(Vx1[2L] / d_yz)
    if (abs(d_zx) > 1e-15) t_zx <- asin(Vy1[3L] / d_zx)

    # Step 5: platform tilt beta
    d_sv1 <-  gL * SFy * cos(t_yx) * cos(t_yz)
    d_sv2 <- -gL * SFx * cos(t_xy) * cos(t_xz)

    sv <- 0
    if (abs(d_sv1) > 1e-15) sv <- sv + (Vx2[2L] - by) / d_sv1
    if (abs(d_sv2) > 1e-15) sv <- sv + (Vy2[1L] - bx) / d_sv2
    if (abs(d_sv1) > 1e-15) sv <- sv + (Vz2[2L] - by) / d_sv1

    beta <- asin(sv / 3)

    # convergence check (skip first pass)
    X_curr <- c(SFx, SFy, SFz, bx, by, bz,
                t_xy, t_xz, t_yx, t_yz, t_zx, t_zy, beta)

    if (it > 1L) {
      e0      <- max(abs(X_curr - X_prev))
      e0_hist <- c(e0_hist, e0)
      message(sprintf("  [Xu] iter %4d: e0 = %.4e", it, e0))
      if (e0 < tol) break
    }

    X_prev <- X_curr
  }

  message("  [Xu] converged in ", iters, " iterations")

  # 5. Build calibration matrix (eq. 2)
  A_sb <- matrix(c(
    cos(t_xy) * cos(t_xz),  sin(t_xz), -cos(t_xz) * sin(t_xy),
    -cos(t_yx) * sin(t_yz),  cos(t_yx) * cos(t_yz),  sin(t_yx),
    sin(t_zy), -cos(t_zy) * sin(t_zx),  cos(t_zx) * cos(t_zy)
  ), nrow = 3L, byrow = TRUE)

  S_xu    <- diag(c(SFx, SFy, SFz)) %*% A_sb
  M       <- A_sb
  S       <- c(SFx, SFy, SFz)
  bias    <- c(bx, by, bz)

  .deg <- function(r) r * 180 / pi

  params_xu <- list(
    S_F        = S,
    A_sb       = A_sb,
    bias       = bias,
    angles_deg = list(
      theta_xy = .deg(t_xy), theta_xz = .deg(t_xz),
      theta_yx = .deg(t_yx), theta_yz = .deg(t_yz),
      theta_zx = .deg(t_zx), theta_zy = .deg(t_zy)
    ),
    beta_deg   = .deg(beta),
    S          = S_xu
  )

  message(sprintf("  Scale factors : %.6f  %.6f  %.6f", SFx, SFy, SFz))
  message(sprintf("  Biases        : %.6f  %.6f  %.6f", bx,  by,  bz))
  message(sprintf("  theta_xy = %.4f deg   theta_xz = %.4f deg", .deg(t_xy), .deg(t_xz)))
  message(sprintf("  theta_yx = %.4f deg   theta_yz = %.4f deg", .deg(t_yx), .deg(t_yz)))
  message(sprintf("  theta_zx = %.4f deg   theta_zy = %.4f deg", .deg(t_zx), .deg(t_zy)))
  message(sprintf("  Platform tilt beta = %.4f deg", .deg(beta)))

  list(
    S         = S,               # scaling factor matrix
    M         = M,               # misalignment matrix
    bias      = bias,            # bias vector
    seg_means = seg_means,
    iters     = iters,
    e0_hist   = if (length(e0_hist) > 0L) e0_hist else 0,
    params_xu    = params_xu
  )
}
