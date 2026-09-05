library(kinference)

find_root <- getFromNamespace("find.root", "kinference")

# Increasing and decreasing functions retain their ordinary root estimates.
expect_equal(
  find_root(identity, start = 0, step = .1, target = .6,
            min.x = 0, max.x = 1, fdirection = "increasing"),
  .6, tolerance = 1e-6)
expect_equal(
  find_root(function(x) -x, start = 0, step = .1, target = -.6,
            min.x = 0, max.x = 1, fdirection = "decreasing"),
  .6, tolerance = 1e-6)

# Return an exact hit without evaluating beyond it.
exact_hit <- function(x) {
  if (x > .5) stop("Search passed the exact root")
  x
}
expect_equal(
  find_root(exact_hit, start = 0, step = .5, target = .5,
            fdirection = "increasing"), .5)

# Unreachable targets and zero steps must terminate with an error.
expect_error(
  find_root(identity, start = 0, step = .1, target = 2,
            min.x = 0, max.x = 1, fdirection = "increasing"),
  pattern = "Cannot bracket target")
expect_error(
  find_root(identity, start = 0, step = 0, target = 1,
            fdirection = "increasing"),
  pattern = "Cannot bracket target")

# Invalid function values during bracketing receive an informative error.
for (bad_value in list(NA_real_, Inf, c(1, 2), numeric(0))) {
  bad_function <- function(x) if (x == 0) 0 else bad_value
  expect_error(
    find_root(bad_function, start = 0, step = 1, target = 2,
              fdirection = "increasing"),
    pattern = "Non-finite or non-scalar function value")
}

# Coordinate overflow also terminates cleanly.
expect_error(
  find_root(function(x) 0, start = 1e308, step = 1e308, target = 1,
            fdirection = "increasing"),
  pattern = "Cannot bracket target")

# A tiny initial step exercises the iteration cap before overflow occurs.
expect_error(
  find_root(function(x) 0, start = 0, step = 1e-300, target = 1,
            fdirection = "increasing"),
  pattern = "Cannot bracket target after 1024 iterations")
