test_that("uniform.mle works", {

  expected_zero_mle <- c(0, 1, 50)

  for(x in expected_zero_mle){
    expect_equal(
      uniform.mle(x, a = 0, b = 1, log = T),
      0
      )
    }
})
