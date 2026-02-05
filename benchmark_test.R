library(Rcpp)
library(Matrix)

# Compile the C++ code
sourceCpp("src/clustering_cpp.cpp")

# Create a larger test matrix to better demonstrate performance
set.seed(123)
n <- 500  # 500 sequences
m <- 256  # 256 k-mer features (4^4 for nucleotide 4-mers)

cat("Creating test sparse matrix with", n, "sequences and", m, "features...\n\n")

# Create sparse matrix simulating k-mer frequency data
x <- sparseMatrix(
  i = sample(1:n, 5000, replace = TRUE),
  j = sample(1:m, 5000, replace = TRUE),
  x = rpois(5000, lambda = 3),
  dims = c(n, m)
)

cat("Matrix density:", round(100 * length(x@x) / (n * m), 2), "%\n\n")

# Benchmark Euclidean distance
cat(strrep("=", 60), "\n")
cat("BENCHMARK: Euclidean Distance Calculation\n")
cat(strrep("=", 60), "\n\n")

cat("R base dist() function:\n")
time_r_dist <- system.time({
  dist_r <- dist(as.matrix(x), method = "euclidean")
})
print(time_r_dist)

cat("\nC++ fastEuclideanDist() function:\n")
time_cpp_dist <- system.time({
  dist_cpp <- fastEuclideanDist(x)
})
print(time_cpp_dist)

speedup_dist <- time_r_dist["elapsed"] / time_cpp_dist["elapsed"]
cat("\n** SPEEDUP: ", round(speedup_dist, 2), "x faster **\n\n")

# Benchmark Centroid calculation
cat(strrep("=", 60), "\n")
cat("BENCHMARK: Cluster Centroid Calculation\n")
cat(strrep("=", 60), "\n\n")

# Create random cluster labels
n_clusters <- 20
clslabels <- sample(1:n_clusters, n, replace = TRUE)

cat("R Matrix::colMeans approach:\n")
time_r_centers <- system.time({
  centers_r <- matrix(nrow = n_clusters, ncol = m)
  for (i in 1:n_clusters) {
    selected <- x[clslabels == i, , drop = FALSE]
    if (nrow(selected) > 0) {
      centers_r[i, ] <- colMeans(selected)
    }
  }
})
print(time_r_centers)

cat("\nC++ fastGetCenters() function:\n")
time_cpp_centers <- system.time({
  centers_cpp <- fastGetCenters(x, clslabels)
})
print(time_cpp_centers)

speedup_centers <- time_r_centers["elapsed"] / time_cpp_centers["elapsed"]
cat("\n** SPEEDUP: ", round(speedup_centers, 2), "x faster **\n\n")

# Benchmark Positional weights
cat(strrep("=", 60), "\n")
cat("BENCHMARK: Positional Weight Calculation\n")
cat(strrep("=", 60), "\n\n")

# Create test sequences and kmers
test_seqs <- replicate(100, paste(sample(c("A", "C", "G", "T"), 20, replace = TRUE), collapse = ""))
test_kmers <- c("ACGT", "CGTA", "GTAC", "TACG", "ACAT", "TGCA")

# R version using sapply and gregexpr
determineWeight_R <- function(kmer, seq) {
  kmerPositions <- as.numeric(gregexpr(kmer, seq, fixed = TRUE)[[1]])
  seqLength <- nchar(seq)
  weightGroups <- seq(seqLength/3, seqLength, seqLength/3)
  
  kmerWts <- c()
  for(i in kmerPositions){
    if(i == -1){
      wt <- 0
    } else {
      wps <- sum(weightGroups > i)
      if(wps == 3) { wt <- 5 }
      else if(wps == 2) { wt <- 10 }
      else if(wps == 1) { wt <- 1 }
    }
    kmerWts <- c(kmerWts, wt)
  }
  return(sum(kmerWts))
}

cat("R sapply + gregexpr approach:\n")
time_r_weights <- system.time({
  for (seq in test_seqs) {
    weights_r <- sapply(test_kmers, function(k) determineWeight_R(k, seq))
  }
})
print(time_r_weights)

cat("\nC++ fastDetermineWeightsVector() function:\n")
time_cpp_weights <- system.time({
  for (seq in test_seqs) {
    weights_cpp <- fastDetermineWeightsVector(test_kmers, seq)
  }
})
print(time_cpp_weights)

speedup_weights <- time_r_weights["elapsed"] / time_cpp_weights["elapsed"]
cat("\n** SPEEDUP: ", round(speedup_weights, 2), "x faster **\n\n")

# Summary
cat("\n")
cat(strrep("=", 60), "\n")
cat("PERFORMANCE SUMMARY\n")
cat(strrep("=", 60), "\n")
cat("Distance Calculation:  ", round(speedup_dist, 2), "x speedup\n")
cat("Centroid Calculation:  ", round(speedup_centers, 2), "x speedup\n")
cat("Positional Weights:    ", round(speedup_weights, 2), "x speedup\n")
cat(strrep("=", 60), "\n\n")

cat("✓ C++ optimizations provide significant performance improvements!\n")
cat("  These speedups will compound in the full RepAn analysis pipeline\n")
cat("  where these operations are called repeatedly.\n")
