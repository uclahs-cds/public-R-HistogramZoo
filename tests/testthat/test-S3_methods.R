context("S3 Histogram methods")

# Writing a set of tests to test S3 methods for fitting histograms
# Essentially testing that the coordinate-based systems work well

# Coordinate based functions
# 1. fit_distributions
# 2. fit_uniform
# 3. identify_uniform_segment

# Index based functions
# 4. find_local_optima
# 5. FTC
# 6. meaningful_gaps_local, find_all_meaningful_gap

# Setting up for coordinate-based functions
set.seed(314)

# Normal distribution
norm_data <- rnorm(1000, mean = 10, sd = 2)
norm_histogram <- observations_to_histogram(norm_data, histogram_bin_width = 2)
norm_numeric <- norm_histogram$histogram_data

norm_genomichistogram <- GenomicHistogram(
  histogram_data = norm_numeric,
  interval_start = norm_histogram$interval_start+1,
  interval_end = norm_histogram$interval_end,
  chr = "chr1",
  strand = "*",
  consecutive_start = norm_histogram$interval_start+1,
  consecutive_end = norm_histogram$interval_end
)

# Uniform distribution
unif_data <- runif(1000, min = 10, max = 20)
unif_histogram <- observations_to_histogram(unif_data, histogram_bin_width = 2)
unif_numeric <- unif_histogram$histogram_data

unif_genomichistogram <- GenomicHistogram(
  histogram_data = unif_numeric,
  interval_start = unif_histogram$interval_start+1,
  interval_end = unif_histogram$interval_end,
  chr = "chr1",
  strand = "*",
  consecutive_start = unif_histogram$interval_start+1,
  consecutive_end = unif_histogram$interval_end
)


# A segment is uniform
unif_trim <- c(1, 1, rep(2, 7), 1)
unif_trim_histogram <- Histogram(
  histogram_data = unif_trim,
  interval_start = seq(10, 28, 2),
  interval_end = seq(12, 30, 2)
)
unif_trim_genomichistogram <- GenomicHistogram(
  histogram_data = unif_trim,
  interval_start = seq(10, 28, 2) + 1,
  interval_end = seq(12, 30, 2),
  chr = "chr1",
  strand = "*",
  consecutive_start = seq(10, 28, 2) + 1,
  consecutive_end = seq(12, 30, 2)
)

  
test_that("S3: fit_distributions", {
  
  # Numeric
  num_res <- fit_distributions(norm_numeric, dist = "norm", metric = "jaccard")
  num_par <- num_res[[1]]$par
  
  # shifted back to 1
  expect_true(
    num_par['mean'] < 4 & num_par['mean'] > 3
  )
  
  expect_true(
    num_par['sd'] < 1 & num_par['sd'] > 0.9
  )
  
  # Histogram
  histogram_res <- fit_distributions(norm_histogram, dist = "norm", metric = "jaccard")
  histogram_par <- histogram_res[[1]]$par
  
  expect_true(
    histogram_par['mean'] < 10.5 & histogram_par['mean'] > 9.5
  )
  
  expect_true(
    histogram_par['sd'] < 2.1 & histogram_par['sd'] > 1.9
  )
  
  # GenomicHistogram
  genomichistogram_res <- fit_distributions(norm_genomichistogram, dist = "norm", metric = "jaccard")
  genomichistogram_par <- genomichistogram_res[[1]]$par
  
  expect_true(
    genomichistogram_par['mean'] < 10.5 & genomichistogram_par['mean'] > 9.5
  )
  
  expect_true(
    genomichistogram_par['sd'] < 2.1 & genomichistogram_par['sd'] > 1.9
  )
  
})

test_that("S3: fit_uniform", {
  
  # Numeric
  num_res <- fit_uniform(unif_numeric, metric = "jaccard")
  num_dens <- num_res$dens

  expect_true(
    num_res$value > 0.9
  ) 
  
  expect_equal(
    num_dens(3),
    mean(unif_numeric)
  )
  
  expect_equal(
    num_dens(6),
    0
  )
  
  # Histogram
  histogram_res <- fit_uniform(unif_histogram, metric = "jaccard")
  histogram_dens <- histogram_res$dens
  
  expect_true(
    histogram_res$value > 0.9
  ) 
  
  expect_equal(
    histogram_dens(14),
    sum(unif_histogram$histogram_data)/length(unif_histogram)
  )
  
  expect_equal(
    histogram_dens(6),
    0
  )
  
  # GenomicHistogram
  genomichistogram_res <- fit_uniform(unif_genomichistogram, metric = "jaccard")
  genomichistogram_dens <- genomichistogram_res$dens
  
  expect_true(
    genomichistogram_res$value > 0.9
  ) 
  
  expect_equal(
    genomichistogram_dens(14),
    sum(unif_genomichistogram$histogram_data)/length(unif_histogram)
  )
  
  expect_equal(
    genomichistogram_dens(6),
    0
  )
  
})

test_that("S3: identify_uniform_segment", {
  
  # Numeric
  num_res <- identify_uniform_segment(unif_trim, max_sd_size = 0)
  
  expect_equal(num_res$seg_start, 3)
  expect_equal(num_res$seg_end, 9)
  
  # Histogram
  histogram_res <- identify_uniform_segment(unif_trim_histogram, max_sd_size = 0) 
  
  expect_equal(histogram_res$dens(9), 0)
  expect_equal(histogram_res$dens(17), 2)
  
  # GenomicHistogram
  genomichistogram_res <- identify_uniform_segment(unif_trim_genomichistogram, max_sd_size = 0)
  
  expect_equal(genomichistogram_res$dens(9), 0)
  expect_equal(genomichistogram_res$dens(17), 2)
  
})

local_optima <- c(1, 2, 3, 2, 1)
local_optima_histogram <- Histogram(
  histogram_data = local_optima,
  interval_start = seq(10, 18, 2),
  interval_end = seq(12, 20, 2)
)
local_optima_genomichistogram <- GenomicHistogram(
  histogram_data = local_optima,
  interval_start = seq(10, 18, 2) + 1,
  interval_end = seq(12, 20, 2)
)

test_that("S3: find_local_optima", {
  
  num_res <- find_local_optima(local_optima)
  
  expect_equal(num_res$min_ind, c(1,5))
  expect_equal(num_res$max_ind, c(3))
  
  histogram_res <- find_local_optima(local_optima_histogram)
  
  expect_equal(histogram_res$min_ind, c(1,5))
  expect_equal(histogram_res$max_ind, c(3))
  
  genomichistogram_res <- find_local_optima(local_optima_genomichistogram)

  expect_equal(genomichistogram_res$min_ind, c(1,5))
  expect_equal(genomichistogram_res$max_ind, c(3))
  
})

set.seed(314)
ftc_data <- rnorm(1000, mean = 20, sd = 5)
ftc_histogram <- observations_to_histogram(ftc_data, histogram_bin_width = 2)
ftc_numeric <- ftc_histogram$histogram_data
ftc_genomichistogram <- GenomicHistogram(
  ftc_numeric,
  interval_start = ftc_histogram$interval_start + 1,
  interval_end = ftc_histogram$interval_end,
  chr = "chr1",
  strand = "*",
  consecutive_start = ftc_histogram$interval_start + 1,
  consecutive_end = ftc_histogram$interval_end
)
ftc_local_optima <- sort(unlist(find_local_optima(ftc_numeric)))

test_that("S3: FTC", {
  
  num_res <- ftc(ftc_numeric)
  expect_equal(num_res, c(1,17))
  
  histogram_res <- ftc(ftc_histogram)
  expect_equal(histogram_res, c(1,17))
  
  genomic_histogram_res <- ftc(ftc_genomichistogram)
  expect_equal(genomic_histogram_res, c(1,17))
  
})

test_that("S3: meaningful_gap", {
  
  num_res <- find_all_meaningful_gap(ftc_numeric, ftc_local_optima)
  
  num_res_gaps <- meaningful_gaps_local(
    ftc_numeric, 
    seg_points = c(1, length(ftc_numeric)), 
    change_points = ftc_local_optima,
    min_gap = 0
  )
  
  expect_equal(num_res_gaps$start[1], 16)
  expect_equal(num_res_gaps$end[1], 17)
  
  histogram_res <- find_all_meaningful_gap(ftc_histogram, ftc_local_optima)
  
  histogram_res_gaps <- meaningful_gaps_local(
    ftc_histogram, 
    seg_points = c(1, length(ftc_numeric)), 
    change_points = ftc_local_optima,
    min_gap = 0
  )
  
  expect_equal(histogram_res_gaps$start[1], 16)
  expect_equal(histogram_res_gaps$end[1], 17)
  
  genomichistogram_res <- find_all_meaningful_gap(ftc_genomichistogram, ftc_local_optima)
  
  genomichistogram_res_gaps <- meaningful_gaps_local(
    ftc_genomichistogram, 
    seg_points = c(1, length(ftc_numeric)), 
    change_points = ftc_local_optima,
    min_gap = 0
  )
  
  expect_equal(genomichistogram_res_gaps$start[1], 16)
  expect_equal(genomichistogram_res_gaps$end[1], 17)
  
})