# vmf_kde ----
test_that("vmf_kde validates input matrix class and columns", {
  x_df <- data.frame(x = 1, y = 1, z = 1)
  expect_error(vmf_kde(x_df), "must be a matrix-like object")

  x_bad_ncol <- matrix(1:4, ncol = 4, dimnames = list(NULL, c("x", "y", "z", "w")))
  expect_error(vmf_kde(x_bad_ncol), "must be a matrix-like object")
})

test_that("vmf_kde validates column names", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("a", "b", "c")))
  expect_error(vmf_kde(x), "named columns x, y, and/or z")
})

test_that("vmf_kde validates weights argument", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  expect_error(vmf_kde(x, weights = "a"), "'weights' can only be a numeric vector")
  expect_error(vmf_kde(x, weights = list(1)), "'weights' can only be a numeric vector")
})

test_that("vmf_kde validates weights recycling length", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  # nrow(x) = 10, weights length 3 does not divide evenly
  expect_error(vmf_kde(x, weights = c(1, 2, 3)), "cannot be recycled")
})

test_that("vmf_kde returns a matrix with correct columns (3D)", {
  set.seed(42)
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  out <- vmf_kde(x, n_grid = 100)
  expect_true(is.matrix(out))
  expect_equal(colnames(out), c("x", "y", "z", "density"))
  expect_equal(nrow(out), 100)
})

test_that("vmf_kde returns a matrix with correct columns (2D)", {
  set.seed(42)
  x <- matrix(rnorm(20), ncol = 2, dimnames = list(NULL, c("x", "y")))
  out <- vmf_kde(x, n_grid = 100)
  expect_true(is.matrix(out))
  expect_equal(colnames(out), c("x", "y", "density"))
  expect_equal(nrow(out), 100)
})

test_that("vmf_kde grid points lie on the unit sphere/circle", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  out <- vmf_kde(x, n_grid = 50)
  norms <- sqrt(rowSums(out[, c("x", "y", "z")]^2))
  expect_equal(norms, rep(1, 50), tolerance = 1e-8)

  x2 <- matrix(rnorm(20), ncol = 2, dimnames = list(NULL, c("x", "y")))
  out2 <- vmf_kde(x2, n_grid = 50)
  norms2 <- sqrt(rowSums(out2[, c("x", "y")]^2))
  expect_equal(norms2, rep(1, 50), tolerance = 1e-8)
})

test_that("vmf_kde density peaks near the concentration of the data (3D)", {
  set.seed(42)
  # simulate points tightly clustered around +z
  n <- 200
  base <- matrix(c(0, 0, 1), nrow = n, ncol = 3, byrow = TRUE)
  noise <- matrix(rnorm(n * 3, sd = 0.05), ncol = 3)
  x <- base + noise
  colnames(x) <- c("x", "y", "z")

  out <- vmf_kde(x, n_grid = 2000, kappa = 20)
  peak <- out[which.max(out[,"density"]), c("x", "y", "z")]

  expect_equal(as.numeric(peak), c(0, 0, 1), tolerance = 0.1)
})

test_that("vmf_kde density peaks near the concentration of the data (2D)", {
  set.seed(42)
  n <- 200
  base <- matrix(c(1, 0), nrow = n, ncol = 2, byrow = TRUE)
  noise <- matrix(rnorm(n * 2, sd = 0.05), ncol = 2)
  x <- base + noise
  colnames(x) <- c("x", "y")

  out <- vmf_kde(x, n_grid = 1440, kappa = 20)
  peak <- out[which.max(out[,"density"]), c("x", "y")]

  expect_equal(as.numeric(peak), c(1, 0), tolerance = 0.1)
})

test_that("vmf_kde normalise = FALSE returns unnormalised sums", {
  set.seed(42)
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  out_raw <- vmf_kde(x, n_grid = 50, normalise = FALSE)
  out_norm <- vmf_kde(x, n_grid = 50, normalise = TRUE)

  # normalised and raw densities should be proportional, not identical
  expect_false(isTRUE(all.equal(out_raw[,"density"], out_norm[,"density"])))
})

test_that("vmf_kde filters out observations below norm_filter", {
  x <- rbind(
    matrix(rnorm(30), ncol = 3),
    c(1e-12, 1e-12, 1e-12) # near-zero row, should be filtered
  )
  colnames(x) <- c("x", "y", "z")
  # should not error and should effectively ignore the near-zero row
  expect_silent(out <- vmf_kde(x, n_grid = 50))
})

test_that("vmf_kde accepts scalar and vector weights", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  expect_silent(out1 <- vmf_kde(x, n_grid = 50, weights = 1))
  expect_silent(out2 <- vmf_kde(x, n_grid = 50, weights = rep(2, 10)))
  # doubling all weights should double the (unnormalised) density
  out1_raw <- vmf_kde(x, n_grid = 50, weights = 1, normalise = FALSE)
  out2_raw <- vmf_kde(x, n_grid = 50, weights = 2, normalise = FALSE)
  expect_equal(out2_raw[,"density"], out1_raw[,"density"] * 2, tolerance = 1e-8)
})

# rotation_to_align ----

test_that("rotation_to_align validates input matrix class and columns", {
  x_df <- data.frame(x = 1, y = 1, z = 1)
  expect_error(rotation_to_align(x_df, align_to = "+z"), "must be a matrix-like object")

  x_bad_names <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("a", "b", "c")))
  expect_error(rotation_to_align(x_bad_names, align_to = "+z"), "named columns x, y, and/or z")
})

test_that("rotation_to_align validates fixed_ax argument", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  expect_error(
    rotation_to_align(x, align_to = "+z", fixed_ax = "w"),
    "must be either NULL or a character string"
  )
  expect_error(
    rotation_to_align(x, align_to = "+z", fixed_ax = c("x", "y")),
    "must be either NULL or a character string"
  )
})

test_that("rotation_to_align messages when fixed_ax given for 2D input", {
  x <- matrix(rnorm(20), ncol = 2, dimnames = list(NULL, c("x", "y")))
  expect_message(
    rotation_to_align(x, align_to = "+x", fixed_ax = "x"),
    "ignored when x is a 2-column matrix"
  )
})

test_that("rotation_to_align validates align_to character input", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  expect_error(rotation_to_align(x, align_to = c("+x", "+y")), "length 1")
  expect_error(rotation_to_align(x, align_to = "+q"), "Invalid 'align_to'")
})

test_that("rotation_to_align rejects align_to matching fixed_ax", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  expect_error(
    rotation_to_align(x, align_to = "+x", fixed_ax = "x"),
    "'fixed_ax' can't be the same axis"
  )
})

test_that("rotation_to_align validates align_to numeric length", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  expect_error(
    rotation_to_align(x, align_to = c(1, 0)),
    "must match the dimensionality"
  )
})

test_that("rotation_to_align validates align_to type", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  expect_error(rotation_to_align(x, align_to = TRUE),
               "must be either a numeric vector")
})

test_that("rotation_to_align validates align_secondary type", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  expect_error(rotation_to_align(x, align_to = "+x",
                                 align_secondary = TRUE),
               "must be either a numeric vector")
})

test_that("rotation_to_align ignores align_secondary in 2D case", {
  x <- matrix(rnorm(20), ncol = 2, dimnames = list(NULL, c("x", "y")))
  expect_message(
    rotation_to_align(x, align_to = "+x", align_secondary = "+y"),
    "ignored in the two-dimensional case"
  )
})

test_that("rotation_to_align rejects character align_secondary with length more than 1", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  expect_error(
    rotation_to_align(x, align_to = "+x", align_secondary = c("+y", "+z")),
    "'align_secondary' must be of length 1"
  )
})

test_that("rotation_to_align rejects invalid character align_secondary", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  expect_error(
    rotation_to_align(x, align_to = "+x", align_secondary = c("+q")),
    "Invalid 'align_secondary' argument."
  )
})

test_that("rotation_to_align checks for dimensionality of align_secondary matching input matrix", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  expect_error(
    rotation_to_align(x, align_to = "+x", align_secondary = c(0,1)),
    "must match the dimensionality of the rotation"
  )
})

test_that("rotation_to_align rejects parallel align_to / align_secondary", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  expect_error(
    rotation_to_align(x, align_to = c(1, 0, 0), align_secondary = c(1, 0, 0)),
    "cannot be parallel"
  )
  expect_error(
    rotation_to_align(x, align_to = "+x", align_secondary = "+x"),
    "cannot be parallel"
  )
})

test_that("rotation_to_align warns when align_secondary is close to align_to", {
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  n <- 100
  # numeric vectors close together (small angle) rather than the canonical
  # +x/+y case, so it passes the parallel check but trips the 45-degree one
  expect_warning(
    rotation_to_align(
      x,
      align_to = c(1, 0, 0),
      align_secondary = c(1, 0.1, 0),
      n_grid = n
    ),
    "within 45 degrees"
  )
})

test_that("rotation_to_align returns a 3x3 rotation matrix for 3D input", {
  set.seed(42)
  n <- 200
  base <- matrix(c(0, 0, 1), nrow = n, ncol = 3, byrow = TRUE)
  noise <- matrix(rnorm(n * 3, sd = 0.05), ncol = 3)
  x <- base + noise
  colnames(x) <- c("x", "y", "z")

  R <- rotation_to_align(x, align_to = "+z", n_grid = 2000, kappa = 20)
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(3, 3))
  # rotation matrices are orthogonal: R %*% t(R) = I
  expect_equal(R %*% t(R), diag(3), tolerance = 1e-6)
  expect_equal(det(R), 1, tolerance = 1e-6)
})

test_that("rotation_to_align returns a 2x2 rotation matrix for 2D input", {
  set.seed(42)
  n <- 200
  base <- matrix(c(1, 0), nrow = n, ncol = 2, byrow = TRUE)
  noise <- matrix(rnorm(n * 2, sd = 0.05), ncol = 2)
  x <- base + noise
  colnames(x) <- c("x", "y")

  R <- rotation_to_align(x, align_to = "+x", n_grid = 1440, kappa = 20)
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(2, 2))
  expect_equal(unname(R %*% t(R)), diag(2), tolerance = 1e-6)
})

test_that("rotation_to_align with fixed_ax restricts rotation to 2 axes", {
  set.seed(42)
  n <- 200
  base <- matrix(c(0, 1, 0), nrow = n, ncol = 3, byrow = TRUE)
  noise <- matrix(rnorm(n * 3, sd = 0.05), ncol = 3)
  x <- base + noise
  colnames(x) <- c("x", "y", "z")

  R <- rotation_to_align(x, align_to = "+y", fixed_ax = "z", n_grid = 2000, kappa = 20)
  expect_equal(dim(R), c(3, 3))
  # the fixed axis row/col should look like an identity contribution:
  # rotating a pure z unit vector should leave it (approximately) as z
  z_vec <- c(0, 0, 1)
  expect_equal(as.numeric(R %*% z_vec), z_vec, tolerance = 1e-6)
})

test_that("rotation_to_align actually aligns the density peak to the target", {
  set.seed(42)
  n <- 300
  base <- matrix(c(0, 1, 0), nrow = n, ncol = 3, byrow = TRUE)
  noise <- matrix(rnorm(n * 3, sd = 0.03), ncol = 3)
  x <- base + noise
  colnames(x) <- c("x", "y", "z")

  R <- rotation_to_align(x, align_to = "+x", n_grid = 20000, kappa = 20)
  rotated <- x %*% t(R)
  # after rotation, the mean direction should point close to +x
  mean_dir <- colMeans(rotated)
  mean_dir <- mean_dir / sqrt(sum(mean_dir^2))
  expect_equal(mean_dir, c(1, 0, 0), tolerance = 1e-2)
})

test_that("rotation_to_align with align_secondary produces a valid rotation matrix", {
  set.seed(42)
  n <- 300
  base <- matrix(c(0, 0, 1), nrow = n, ncol = 3, byrow = TRUE)
  noise <- matrix(rnorm(n * 3, sd = 0.05), ncol = 3)
  x <- base + noise
  colnames(x) <- c("x", "y", "z")

  R <- rotation_to_align(
    x,
    align_to = "+z",
    align_secondary = "+x",
    n_grid = 2000,
    kappa = 20
  )
  expect_equal(dim(R), c(3, 3))
  expect_equal(unname(R %*% t(R)), diag(3), tolerance = 1e-6)
})

test_that("rotation_to_align secondary_policy = 'min' produces a valid rotation matrix", {
  set.seed(42)
  n <- 300
  base <- matrix(c(0, 0, 1), nrow = n, ncol = 3, byrow = TRUE)
  noise <- matrix(rnorm(n * 3, sd = 0.05), ncol = 3)
  x <- base + noise
  colnames(x) <- c("x", "y", "z")

  R <- rotation_to_align(
    x,
    align_to = "+z",
    align_secondary = "+x",
    secondary_policy = "min",
    n_grid = 2000,
    kappa = 20
  )
  expect_equal(dim(R), c(3, 3))
  expect_equal(unname(R %*% t(R)), diag(3), tolerance = 1e-6)
})

# apply_rotation ----

test_that("apply_rotation rotates a plain matrix correctly", {
  x <- matrix(c(1, 0, 0), nrow = 1, dimnames = list(NULL, c("x", "y", "z")))
  # 90 degree rotation about z: x-axis -> y-axis
  R <- matrix(c(0, -1, 0,
                1,  0, 0,
                0,  0, 1), nrow = 3, byrow = TRUE)
  rotated <- apply_rotation(x, R)
  expect_equal(as.numeric(rotated), c(0, 1, 0), tolerance = 1e-8)
})

test_that("apply_rotation with identity matrix returns input unchanged", {
  set.seed(42)
  x <- matrix(rnorm(9), ncol = 3)
  R <- diag(3)
  expect_equal(apply_rotation(x, R), x %*% t(R))
  expect_equal(apply_rotation(x, R), x)
})

test_that("apply_rotation preserves vector norms", {
  set.seed(42)
  x <- matrix(rnorm(30), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  R <- rotation_to_align(x, align_to = "+z", n_grid = 200)
  rotated <- apply_rotation(x, R)
  expect_equal(sqrt(rowSums(x^2)), sqrt(rowSums(rotated^2)), tolerance = 1e-8)
})

test_that("apply_rotation preserves aclrtm_accelerometery class and attributes", {
  skip_if_not(exists("new_accelerometery"), "new_accelerometery constructor not available")

  x <- matrix(rnorm(9), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  x <- new_accelerometery(x, sampling_rate = 100, start_time = 0)

  R <- diag(3)
  rotated <- apply_rotation(x, R)

  expect_s3_class(rotated, "aclrtm_accelerometery")
  expect_equal(attr(rotated, "sampling_rate"), 100)
  expect_equal(attr(rotated, "start_time"), 0)
})
