test_that(".parse_axes converts a valid axes string to a character vector", {
  expect_equal(.parse_axes("xyz"), c("x", "y", "z"))
  expect_equal(.parse_axes("xz"), c("x", "z"))
  expect_equal(.parse_axes("y"), "y")
})

test_that(".parse_axes rejects invalid input", {
  expect_error(.parse_axes(123), "single character string")
  expect_error(.parse_axes(c("x", "y")), "single character string")
  expect_error(.parse_axes(""), "must not be an empty string")
  expect_error(.parse_axes("xw"), "only contain")
  expect_error(.parse_axes("xx"), "duplicate")
})

test_that(".infer_axes returns the correct labels for 1-3 columns", {
  expect_equal(.infer_axes(1), "x")
  expect_equal(.infer_axes(2), c("x", "y"))
  expect_equal(.infer_axes(3), c("x", "y", "z"))
})

test_that(".infer_axes errors outside 1-3", {
  expect_error(.infer_axes(4), "not 1, 2, or 3")
  expect_error(.infer_axes(0), "not 1, 2, or 3")
})

test_that(".infer_axes_from_length picks the largest axis count that divides evenly", {
  expect_equal(.infer_axes_from_length(12), c("x", "y", "z")) # divisible by 3
  expect_equal(.infer_axes_from_length(8), c("x", "y"))       # not by 3, but by 2
  expect_equal(.infer_axes_from_length(7), "x")                # neither
})

test_that(".as_integer_strict accepts values within tolerance and rejects others", {
  expect_equal(.as_integer_strict(4), 4L)
  expect_equal(.as_integer_strict(4 + 1e-10), 4L)  # within tolerance
  expect_error(.as_integer_strict(4.3), "not an integer value")
})

test_that(".gcd2 computes the greatest common divisor", {
  expect_equal(.gcd2(12, 8), 4)
  expect_equal(.gcd2(7, 3), 1)
  expect_equal(.gcd2(5, 0), 5)  # base case
})
