library(Rcpp)
library(Matrix)

# Compile the C++ code
sourceCpp("src/clustering_cpp.cpp")

# Test with a simple sparse matrix
set.seed(123)
n <- 50
m <- 100
x <- sparseMatrix(
  i = sample(1:n, 200, replace = TRUE),
  j = sample(1:m, 200, replace = TRUE),
  x = rnorm(200),
  dims = c(n, m)
)

# Test Euclidean distance
cat("Testing Euclidean distance calculation...\n")
system.time(dist_eu <- fastEuclideanDist(x))
cat("Euclidean distance matrix computed: ", nrow(dist_eu), "x", ncol(dist_eu), "\n")

# Test Cosine distance
cat("\nTesting Cosine distance calculation...\n")
system.time(dist_cos <- fastCosineDist(x))
cat("Cosine distance matrix computed: ", nrow(dist_cos), "x", ncol(dist_cos), "\n")

# Test centroid calculation
cat("\nTesting centroid calculation...\n")
clslabels <- sample(1:5, n, replace = TRUE)
system.time(centers <- fastGetCenters(x, clslabels))
cat("Centroids computed: ", nrow(centers), "x", ncol(centers), "\n")

# Test positional weight calculation
cat("\nTesting positional weight calculation...\n")
test_seq <- "ACGTACGTACGT"
test_kmers <- c("ACGT", "CGTA", "GTAC")
system.time(weights <- fastDetermineWeightsVector(test_kmers, test_seq))
cat("Weights computed for", length(test_kmers), "kmers:", paste(weights, collapse=", "), "\n")

cat("\n✓ All C++ functions compiled and tested successfully!\n")
