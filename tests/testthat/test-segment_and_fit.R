context("segment_and_fit")

test_that("segment_and_fit works without removing low entropy regions", {

  x_norm_mix <- Histogram(
    c(`-3` = 0L, `-2` = 0L, `-1` = 3L, `0` = 18L, `1` = 33L, `2` = 30L,
      `3` = 12L, `4` = 4L, `5` = 0L, `6` = 0L, `7` = 0L, `8` = 0L,
      `9` = 0L, `10` = 0L, `11` = 0L, `12` = 0L, `13` = 0L, `14` = 0L,
      `15` = 0L, `16` = 0L, `17` = 0L, `18` = 0L, `19` = 0L, `20` = 0L,
      `21` = 0L, `22` = 0L, `23` = 0L, `24` = 0L, `25` = 0L, `26` = 7L,
      `27` = 5L, `28` = 2L, `29` = 7L, `30` = 7L, `31` = 4L, `32` = 2L,
      `33` = 9L, `34` = 5L, `35` = 2L, `36` = 0L)
    )
  res <- segment_and_fit(x_norm_mix, remove_low_entropy = F)
  res_trunc <- segment_and_fit(x_norm_mix, truncated_models = TRUE, remove_low_entropy = F)

  # Expect 2 segments
  # segment one should be [3,8]
  # segment two should be [30, 39]

  p_expect <- matrix(c(3,8, 30, 39), nrow = 2, ncol = 2, byrow = TRUE)
  segment_points_matrix <- as.matrix(res$p[, c('start', 'end')])
  segment_points_trunc_matrix <- as.matrix(res_trunc$p[, c('start', 'end')])
  expect_equivalent(segment_points_matrix, p_expect)
  expect_equivalent(segment_points_trunc_matrix, p_expect)
})
