context("C++ Optimization Functions")

test_that("fastGetCenters computes correct centroids", {
  skip_if_not_installed("Rcpp")
  skip_if_not_installed("Matrix")
  
  library(Matrix)
  
  # Create test sparse matrix
  x <- sparseMatrix(
    i = c(1, 1, 2, 2, 3, 3),
    j = c(1, 2, 1, 3, 2, 3),
    x = c(1, 2, 3, 4, 5, 6),
    dims = c(3, 3)
  )
  
  clslabels <- c(1, 1, 2)
  
  # Compute centroids
  centers <- fastCenters(x, clslabels)
  
  # Check dimensions
  expect_equal(nrow(centers), 2)  # 2 clusters
  expect_equal(ncol(centers), 3)  # 3 features
  
  # Check values - cluster 1 should be mean of rows 1 and 2
  expect_equal(centers[1, 1], 2)  # (1 + 3) / 2
  expect_equal(centers[1, 2], 1)  # (2 + 0) / 2
  expect_equal(centers[1, 3], 2)  # (0 + 4) / 2
  
  # Cluster 2 should be row 3
  expect_equal(centers[2, 1], 0)
  expect_equal(centers[2, 2], 5)
  expect_equal(centers[2, 3], 6)
})

test_that("fastDetermineWeightsVector calculates positional weights correctly", {
  skip_if_not_installed("Rcpp")
  
  # Test sequence
  seq <- "ACGTACGTACGT"  # 12 characters
  # Thirds: 1-4 (v-area, wt=5), 5-8 (middle, wt=10), 9-12 (j-area, wt=1)
  
  # Test kmer in v-area
  kmers <- c("ACGT")  # at positions 1 and 5
  weights <- fastDetermineWeightsVector(kmers, seq)
  # Position 1 (v-area) = 5, Position 5 (middle) = 10
  expect_equal(weights[1], 15)
  
  # Test multiple kmers
  kmers <- c("ACGT", "CGTA", "GTAC")
  weights <- fastDetermineWeightsVector(kmers, seq)
  expect_length(weights, 3)
  expect_true(all(weights >= 0))
})

test_that("fastDetermineWeightsVector handles missing kmers", {
  skip_if_not_installed("Rcpp")
  
  seq <- "AAAAAAAAAA"
  kmers <- c("TTTT", "GGGG")  # Not present in sequence
  
  weights <- fastDetermineWeightsVector(kmers, seq)
  expect_equal(weights, c(0, 0))
})

test_that("C++ functions handle edge cases", {
  skip_if_not_installed("Rcpp")
  skip_if_not_installed("Matrix")
  
  library(Matrix)
  
  # Empty cluster labels (only zeros)
  x <- sparseMatrix(i = c(1, 2), j = c(1, 1), x = c(1, 2), dims = c(2, 2))
  clslabels <- c(0, 0)
  centers <- fastCenters(x, clslabels)
  expect_equal(nrow(centers), 0)  # No valid clusters
  
  # Single cluster
  clslabels <- c(1, 1)
  centers <- fastCenters(x, clslabels)
  expect_equal(nrow(centers), 1)
  
  # Empty sequence
  weights <- fastDetermineWeightsVector(c("ACGT"), "")
  expect_equal(weights[1], 0)
})
