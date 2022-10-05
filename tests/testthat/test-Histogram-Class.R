context("Histogram-class")

test_that("Generating a valid Histogram object", {

  x <- Histogram(
    histogram_data = c(1, 2, 3, 4, 3, 3, 1)
  )

  # Interval starts
  expect_equal(x$interval_start, 1:7)

  # Interval ends
  expect_equal(x$interval_end, 2:8)

  # region_ids
  expect_equal(x$region_id, "1-8")

})

test_that("Generating an invalid Histogram object", {

  # Invalid data
  expect_error(
    Histogram(
      c(-1, 2, 2, -1)
    )
  )

  # Invalid intervals
  expect_error(
    Histogram(
      runif(10),
      interval_start = 1:10,
      interval_end = 1:10 - 1
    )
  )

  # Non ordered intervals
  expect_error(
    Histogram(
      runif(2),
      interval_start = c(2, 1),
      interval_end = c(3, 1)
    )
  )

  # Overlapping intervals
  expect_error(
    Histogram(
      runif(2),
      interval_start = c(1, 1),
      interval_end = c(2, 3)
    )
  )

})

test_that("Histogram-class methods", {

  x <- Histogram(
    histogram_data = c(1, 2, 3, 4, 3, 3, 1)
  )

  # Print
  # TODO: Add the expected print output
  expect_invisible(print(x))

  # Length
  expect_length(x , 7)

  # Extract Indices
  x_subset <- x[2:5]
  expect_equal(
    x_subset$histogram_data,
    c(2, 3, 4, 3)
  )

  # reassign_ids
  x_reassign <- reassign_region_id(x, "TEST")
  expect_equal(x_reassign$region_id, "TEST")

})
