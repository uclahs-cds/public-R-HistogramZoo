context("gamma_flip")

test_that("gamma_flip", {
  
  # rgamma_flip
  expect_error(
    rgamma_flip(100, shape = 10),
    NA
  )
  
  expect_error(
    rgamma_flip(100, shape = 10, offset = 20),
    NA
  )
  
  # dgamma_flip
  dgamma_data <- dgamma(1:20, shape = 10)
  dgamma_flip_data <- dgamma_flip(-c(1:20), shape = 10)
  dgamma_flip_offset_data <- dgamma_flip(rev(0:19), shape = 10, offset =20)

  expect_equal(dgamma_data, dgamma_flip_data)
  expect_equal(dgamma_data, dgamma_flip_offset_data)
  
  # pgamma_flip
  pgamma_data <- pgamma(1:20, shape = 10)
  pgamma_flip_data <- pgamma_flip(-c(1:20), shape = 10, lower.tail = F)
  pgamma_flip_offset_data <- pgamma_flip(rev(0:19), shape = 10, offset = 20, lower.tail = F)
  
  expect_equal(pgamma_data, pgamma_flip_data)
  expect_equal(pgamma_data, pgamma_flip_offset_data)
  
  # qgamma_flip
  expect_equal(
    pgamma(qgamma(seq(0, 0.9, 0.1), shape = 10), shape = 10), 
    seq(0, 0.9, 0.1)
  )
  
  expect_equal(
    pgamma_flip(qgamma_flip(seq(0, 0.9, 0.1), shape = 10), shape = 10),
    seq(0, 0.9, 0.1)
  )
  
  expect_equal(
    pgamma_flip(qgamma_flip(seq(0, 0.9, 0.1), shape = 10, offset = 20), shape = 10, offset = 20),
    seq(0, 0.9, 0.1)
  )
  
})

test_that("tgamma_flip ", {
  
  # rtgamma_flip
  expect_error(
    rtgamma_flip(10000, shape = 10, a = 5, b = 20),
    NA
  )
  
  expect_error(
    rtgamma_flip(10000, shape = 10),
    "Need to supply upper bound."
  )
  
  # dtgamma_flip
  dtgamma_data <- dtgamma(5:20, shape = 10, a = 5, b = 20)
  dtgamma_flip_data <- dtgamma_flip(5:20, shape = 10, a = 5, b = 20)
  
  expect_equal(
    dtgamma_data,
    rev(dtgamma_flip_data)
  )
  
  # ptgamma_flip
  ptgamma_data <- ptgamma(5:20, shape = 10, a = 5, b = 20)
  ptgamma_flip_data <- ptgamma_flip(rev(5:20), shape = 10, lower.tail = F, a = 5, b = 20)
  
  expect_equal(ptgamma_data, ptgamma_flip_data)

  # qtgamma_flip
  expect_equal(
    ptgamma(qtgamma(seq(0, 0.9, 0.1), shape = 10), shape = 10), 
    seq(0, 0.9, 0.1)
  )
  
  expect_equal(
    ptgamma(qtgamma(seq(0, 0.9, 0.1), shape = 10, a = 5, b = 20), shape = 10, a = 5, b = 20),
    seq(0, 0.9, 0.1)
  )
  
  expect_equal(
    ptgamma_flip(qtgamma_flip(seq(0, 0.9, 0.1), shape = 10, a = 5, b = 20), shape = 10, a = 5, b = 20),
    seq(0, 0.9, 0.1)
  )
  
})