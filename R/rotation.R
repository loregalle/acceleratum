#' Von Mises-Fisher kernel density estimates
#'
#' Calculates kernel density estimates on a Fibonacci lattice grid (3D case)
#' or on a regular angular grid (2D case) using a von Mises-Fisher kernel
#' distribution.
#'
#' @param m A numeric matrix with either 2 or 3 columns and
#'   rows as observations.
#' @param n_grid Number of grid points. Default is 1000.
#' @param kappa Concentration parameter for the kernel distribution.
#'   Higher values return sharper kernels.
#' @param normalise Logical. If \code{TRUE} (default), the density is
#'   normalised by multiplying it by \eqn{\frac{C_p(\mathbf{\kappa})}{n}}.
#'   Else, the sum is returned. See details.
#' @return A data frame.
#' @details
#' The probability density function of the von Mises–Fisher distribution for
#' the vector \eqn{\mathbf{x}} in \eqn{p} dimensions is:
#' \deqn{f(\mathbf{x} \mid \boldsymbol{\mu}, \mathbf{\kappa}) = C_p(\mathbf{\kappa})
#' \exp(\mathbf{\kappa}\boldsymbol{\mu}^T\mathbf{x})}
#' where \eqn{\boldsymbol{\mu}} and \eqn{\mathbf{\kappa}} are the mean and
#' concentration parameters (respectively) of the distribution, and
#' \eqn{C_p(\mathbf{\kappa})} is a constant.
#' The kernel density estimator over \eqn{n} observations \eqn{\mathbf{x}_i}
#' at a candidate direction \eqn{\mathbf{g}} is thus:
#' \deqn{\hat{f}(\mathbf{g}) = \frac{C_p(\mathbf{\kappa})}{n}
#' \sum\limits_{i=1}^{n}\exp(\mathbf{\kappa}\mathbf{g}^T\mathbf{x}_i)}
#' This function picks \code{n_grid} candidate directions on a sphere using a
#' Fibonacci lattice (3D case) or a regular angular grid on a circle (2D case),
#' and calculates the kernel density estimate at each one of them.
#'
#' @export
vmf_kde <- function(m, n_grid = 1000L, kappa = 4, normalise = TRUE) {

  if (!inherits(m, "matrix") || !ncol(m) %in% c(2, 3)) {
    stop("Input m must be a matrix and have either 2 or 3 columns.",
         call. = FALSE)
  }
  if (!all(colnames(m) %in% c("x", "y", "z"))) {
    stop("Input matrix m must have named columns x, y, and/or z.")
  }

  p <- ncol(m)
  norms <- sqrt(rowSums(m^2))
  unit_vec  <- m / norms

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

  density <- rowSums(exp(kappa * (grid %*% t(unit_vec))))
  if (normalise) {
    numer <- kappa^(p/2 - 1)
    denom <- ((2*pi)^(p/2)) * besselI(kappa, p/2 - 1)
    Cpk <- numer/denom
    density <- Cpk * density / nrow(unit_vec)
  }

  out <- cbind(grid, density)
  colnames(out) <- c(colnames(m),"density")
  as.data.frame(out)
}

