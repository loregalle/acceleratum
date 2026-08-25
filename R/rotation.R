#' Von Mises-Fisher kernel density estimates
#'
#' Calculates kernel density estimates on a Fibonacci lattice grid (3D case)
#' or on a regular angular grid (2D case) using a von Mises-Fisher kernel
#' distribution.
#'
#' @param x A numeric matrix with either 2 or 3 columns and
#'   rows as observations.
#' @param n_grid Number of grid points.
#' @param kappa Concentration parameter for the kernel distribution.
#'   Higher values return sharper kernels.
#' @param normalise Logical. If `TRUE` (default), the density is
#'   normalised by multiplying it by \eqn{\frac{C_p(\mathbf{\kappa})}{n}}.
#'   Else, the sum is returned. See details.
#' @param weights Observation weights. Either a single number or a
#'   vector, recycled to `nrow(x)` length.
#' @param norm_filter Minimum vector length. Observations that define
#'   lower vector lengths are filtered out.
#' @details
#' The probability density function of the von Mises–Fisher distribution for
#' the vector \eqn{\mathbf{x}} in \eqn{p} dimensions is:
#' \deqn{f(\mathbf{x} \mid \boldsymbol{\mu}, \mathbf{\kappa}) = C_p(\mathbf{\kappa})
#' \exp(\mathbf{\kappa}\boldsymbol{\mu}^T\mathbf{x})}
#' where \eqn{\boldsymbol{\mu}} and \eqn{\mathbf{\kappa}} are the mean and
#' concentration parameters (respectively) of the distribution, and
#' \eqn{C_p(\mathbf{\kappa})} is a constant.
#' The kernel density estimator over \eqn{n} observations \eqn{\mathbf{x}_i}
#' at a candidate direction \eqn{\mathbf{v}} is thus:
#' \deqn{\hat{f}(\mathbf{v}) = \frac{C_p(\mathbf{\kappa})}{n}
#' \sum\limits_{i=1}^{n}\exp(\mathbf{\kappa}\mathbf{v}^T\mathbf{x}_i)}
#' This function picks `n_grid` candidate directions on a sphere using a
#' Fibonacci lattice (3D case) or a regular angular grid on a circle (2D case),
#' and calculates the kernel density estimate at each one of them.
#'
#' @returns A matrix with dimensions `n_grid, ncol(x)+1`, the added column
#'   being the density estimate at each candidate direction.
#'
#' @export
vmf_kde <- function(x,
                    n_grid = 1440L,
                    kappa = 10,
                    normalise = TRUE,
                    weights = 1,
                    norm_filter = 1e-10) {

  if (!inherits(x, "matrix") || !ncol(x) %in% c(2, 3)) {
    stop("Input x must be a matrix-like object and have either 2 or 3 columns.",
         call. = FALSE)
  }
  if (!all(colnames(x) %in% c("x", "y", "z"))) {
    stop("Input matrix x must have named columns x, y, and/or z.",
         call. = FALSE)
  }
  if (!is.numeric(weights) || !is.vector(weights)) {
    stop("'weights' can only be a numeric vector")
  }
  if (nrow(x) %% length(weights) != 0) {
    stop("'weights' cannot be recycled to nrow(x) length.",
         call. = FALSE)
  }

  p <- ncol(x)
  norms <- sqrt(rowSums(x^2))
  unit_vec <- x / norms
  unit_vec <- unit_vec[norms >= norm_filter,]
  if (length(weights) > 1) {
    weights <- rep(weights, length.out = nrow(x))
    weights <- weights[norms >= norm_filter]
  }

  if (p == 3) {
    golden <- (1 + sqrt(5)) / 2
    i      <- seq_len(n_grid)
    theta  <- acos(1 - 2 * i / n_grid)
    phi    <- 2 * pi * i / golden

    grid <- cbind(
      sin(theta) * cos(phi),
      sin(theta) * sin(phi),
      cos(theta)
    )
  } else {
    theta <- seq(-pi, pi, length.out = n_grid+1)[-1]
    grid <- cbind(
      sin(theta),
      cos(theta)
    )
  }

  density <- colSums(weights * exp(kappa * (unit_vec %*% t(grid))))

  if (normalise) {
    numer <- kappa^(p/2 - 1)
    denom <- ((2*pi)^(p/2)) * besselI(kappa, p/2 - 1)
    Cpk <- numer/denom
    density <- Cpk * density / nrow(unit_vec)
  }

  out <- cbind(grid, density)
  colnames(out) <- c(colnames(x),"density")
  out
}

#' Rotation matrix
#'
#' Return a rotation matrix that aligns the density peak in the sample to
#' the a specified direction.
#'
#' @param x matrix
#' @param align_to Target direction to rotate the density peak to; either a
#'   character string giving a signed axis (e.g. `"+x"`, `"-z"`)
#'   or a numeric vector of length equal to the number of rotation axes.
#' @param align_secondary Optional target direction for a secondary
#'   alignment (same format as `align_to`), used to additionally
#'   orient the rotation around the primary axis in the 3D case.
#' @param secondary_policy Policy for the alignment of the secondary axis.
#'   Ignored if `align_secondary = NULL`. See details.
#' @param fixed_ax Optional character string (`"x"`, `"y"`, `"z"`)
#'   naming an axis to hold fixed, restricting the rotation to
#'   the remaining two axes. Ignored for 2-column input.
#' @param ... other arguments passed to [vmf_kde()]
#' @details
#' The rotation is derived from the empirical directional density of
#' `x`, estimated via [vmf_kde()] on the sphere (3D) or
#' circle (2D). The direction of peak density is identified and a rotation
#' matrix is constructed that maps this peak to `align_to`: via
#' Rodrigues' rotation formula in the 3D case, or a Givens rotation in the
#' 2D case.
#'
#' If `align_secondary` is supplied (3D only), a second density is
#' estimated on the plane perpendicular to the primary peak direction,
#' and the direction of max (or min, depending on \code{secondary_policy})
#' density within that plane
#' is aligned to \code{align_secondary}. This additionally constrains
#' rotation about the primary axis, which is otherwise left unspecified
#' when only \code{align_to} is provided.
#'
#' @returns A rotation matrix
#' @seealso [vmf_kde()] [apply_rotation()]
#' @export
rotation_to_align <- function(x,
                              align_to,
                              align_secondary = NULL,
                              secondary_policy = c("max", "min"),
                              fixed_ax = NULL,
                              ...) {
  secondary_policy <- match.arg(secondary_policy)
  if (!inherits(x, "matrix") || !ncol(x) %in% c(2, 3)) {
    stop("Input x must be a matrix-like object and have either 2 or 3 columns.",
         call. = FALSE)
  }
  if (is.null(colnames(x)) || !all(colnames(x) %in% c("x", "y", "z"))) {
    stop("Input matrix x must have named columns x, y, and/or z.")
  }
  cn <- colnames(x)
  if (!is.null(fixed_ax)) {
    if (ncol(x) == 2) {
      message("'fixed_ax' is ignored when x is a 2-column matrix",
              call. = F)
      fixed_ax <- NULL
      rot_ax <- cn
    } else {
      if (!is.character(fixed_ax) ||
          length(fixed_ax) != 1 ||
          !fixed_ax %in% c("x", "y", "z") ||
          !fixed_ax %in% cn) {
        stop(
          "'fixed_ax' must be either NULL or a character string of length ",
          "one: 'x', 'y', or 'z'. The axis must be present as a column in ",
          "'x'.",
          call. = F
        )
      }
      rot_ax <- setdiff(cn, fixed_ax)
    }
  } else {
    rot_ax <- cn
  }

  n_ax <- length(rot_ax)

  if (!is.null(align_secondary) && n_ax == 2) {
    message("'align_secondary' is ignored in the two-dimensional case.")
    align_secondary <- NULL
  }

  # align_to processing ----
  if (is.character(align_to)) {
    if(length(align_to) != 1) {
      stop("'align_to' must be of length 1 when provided as a character string",
           call. = FALSE)
    }
    align_to <- .canon_face(align_to)
    if (is.na(align_to)) {
      stop("Invalid 'align_to' argument.", call. = FALSE)
    }
    align_to_ax <- substr(align_to, 2,2)
    align_to_sy <- substr(align_to, 1,1)

    if (!is.null(fixed_ax) && align_to_ax == fixed_ax) {
      stop("'fixed_ax' can't be the same axis defined in 'align_to'",
           call. = FALSE)
    }

    align_to <- rep(0,n_ax)
    names(align_to) <- rot_ax
    align_to[align_to_ax] <- ifelse(align_to_sy == "+", 1, -1)
  }
  if (is.numeric(align_to)) {
    if (length(align_to) != length(rot_ax)) {
      stop("When provided as a numeric vector, the length of 'align_to' (",
           length(align_to),
           ") must match the dimensionality of the rotation (",
           length(rot_ax),
           ").",
           call. = FALSE)
    }
    if (is.null(names(align_to))) {
      names(align_to) <- rot_ax
    }
    align_to <- align_to/sqrt(sum(align_to^2))
  } else {
    stop("'align_to' must be either a numeric vector of length equal to ",
         "the number of columns of 'x' or a character string of length 1",
         call. = FALSE)
  }

  # align_secondary processing ----
  if (!is.null(align_secondary)) {
    if (is.character(align_secondary)) {
      if(length(align_secondary) != 1) {
        stop("'align_secondary' must be of length 1 when provided as a ",
             "character string",
             call. = FALSE)
      }
      align_secondary <- .canon_face(align_secondary)
      if (is.na(align_secondary)) {
        stop("Invalid 'align_secondary' argument.", call. = FALSE)
      }
      align_secondary_ax <- substr(align_secondary, 2,2)
      align_secondary_sy <- substr(align_secondary, 1,1)

      if (exists("align_to_ax") && align_to_ax == align_secondary_ax) {
        stop("'align_to' and 'align_secondary' cannot be parallel.",
             call. = F)
      }

      align_secondary <- rep(0,n_ax)
      names(align_secondary) <- rot_ax
      align_secondary[align_secondary_ax] <- ifelse(
        align_secondary_sy == "+", 1, -1
      )
    }
    if (is.numeric(align_secondary)) {
      if (length(align_secondary) != length(rot_ax)) {
        stop("When provided as a numeric vector, the length of ",
             "'align_secondary' (",
             length(align_secondary),
             ") must match the dimensionality of the rotation (",
             length(rot_ax),
             ").",
             call. = FALSE)
      }
      if (is.null(names(align_secondary))) {
        names(align_secondary) <- rot_ax
      }

      align_secondary <- align_secondary / sqrt(sum(align_secondary^2))
      proj_secondary <- .gram_schmidt(align_secondary, align_to)
      sin_theta <- sqrt(sum(proj_secondary^2))

      if (sin_theta < 1e-10) {
        stop("'align_to' and 'align_secondary' cannot be parallel.",
             call. = F)
      }

      if (sin_theta < sin(pi/4)) {
        warning(
          "'align_secondary' is within 45 degrees of 'align_to', ",
          "the secondary alignment may be unreliable.",
          call. = FALSE
        )
      }
      proj_secondary <- proj_secondary / sin_theta
    } else {
      stop("'align_secondary' must be either a numeric vector of length equal ",
           "to the number of columns of 'x' or a character string of length 1",
           call. = FALSE)
    }
  }

  # find directional density. On the sphere if 3D, on the circle if 2D.
  dens <- vmf_kde(x[,rot_ax], ...)
  dens_peak <- dens[which.max(dens[,"density"]), rot_ax, drop = FALSE][1,]

  if (!is.null(align_secondary)) {
    # This branch goes only in the 3D case and align_secondary defined
    # Finds the direction to align the secondary axis to.
    # first, project the density sphere grid onto the plane perpendicular
    # to the axis defined by dens_peak
    dens_proj <- .gram_schmidt(dens[, rot_ax, drop = FALSE],
                               matrix(dens_peak, nrow(dens),
                                      3,
                                      byrow = TRUE))
    colnames(dens_proj) <- rot_ax

    # calculate the vector lengths
    proj_lengths <- sqrt(rowSums(dens_proj^2))

    # To use vmf_kde I need to restrict to 2D, which means I need to
    # work on the plane perpendicular to dens_peak.
    # Conveniently, dens_proj are all vectors orthogonal to dens_peak.
    # so I can use any of these and then find the 2nd axis of the plane
    # by cross product. For safety, I will pick the longest vector
    u <- dens_proj[which.max(proj_lengths), ]
    u <- u/sqrt(sum(u^2))
    v <- .cross(dens_peak, u)
    # now that the new 2D reference system is defined, get all vectors
    # in dens_proj into this new reference
    proj_2d <- dens_proj %*% unname(cbind(u,v))
    colnames(proj_2d) <- c("x", "y") # add names just because vmf_kde requires
    # compute densities. Use vector lengths as weights
    dens2 <- vmf_kde(proj_2d,
                     weights = proj_lengths,
                     norm_filter = 1e-3)
    # based on policy, peak the maximum or minimum density
    if (secondary_policy == "max") {
      sec_selected <- unname(
        dens2[which.max(dens2[,"density"]),1:2, drop = FALSE][1,]
      )
    } else {
      sec_selected <- unname(
        dens2[which.min(dens2[,"density"]),1:2, drop = FALSE][1,]
      )
    }
    # from the 2D reference system to the 3D
    sec_axis <- unname(sec_selected[1]) * u + unname(sec_selected[2]) * v

    # normalisation should not be needed here, but
    # to avoid dragging floating point errors...
    sec_axis <- sec_axis/sqrt(sum(sec_axis^2))

    # find the third axis of both the data frame and the target frame
    a3 <- .cross(dens_peak, sec_axis)
    b3 <- .cross(align_to, proj_secondary)

    # A is the frame defined by axis aligned to dens_peak, sec_axis, a3
    A <- cbind(dens_peak, sec_axis, a3)
    # B is the target frame
    B <- cbind(align_to, proj_secondary, b3)

    # rotation matrix
    R <- B %*% t(A)

  } else if (n_ax == 3) {
    R <- .rodrigues_rotation(dens_peak, align_to)
  } else {
    R <- .givens_rotation(dens_peak, align_to, rot_ax, cn)
  }
  R
}

#' Rotate data
#'
#' Rotate the data according to a rotation matrix
#'
#' @param x Input data matrix, rows are observations.
#' @param R Rotation matrix such as output of [rotation_to_align()].
#'   Transposed internally. See details.
#' @returns Rotated matrix
#' @details In mathematical convention, a vector \eqn{\mathbf{v}}
#'   is a column vector and a
#'   rotation matrix \eqn{R} transforms it to its
#'   rotated equivalent \eqn{\mathbf{v'}} via:
#'   \deqn{\mathbf{v'} = R\mathbf{v}}
#'   However, data matrices conventionally store observations as rows rather
#'   than columns. This function handles the transposition internally, so that
#'   \eqn{R} can be provided in its standard mathematical form and applied
#'   correctly to row-vector data:
#'   \deqn{X' = X R^T}
#'   where \eqn{X} is the input data matrix.
#' @export
apply_rotation <- function(x, R) {
  out <- x %*% t(R)
  if (inherits(x, "aclrtm_accelerometery")) {
    new_accelerometery(
      out,
      sampling_rate = attr(x, "sampling_rate"),
      start_time = attr(x, "start_time")
    )
  } else {
    out
  }
}
