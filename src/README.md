# C++ Optimizations for RepAn

This document describes the C++ optimizations implemented in RepAn to improve computational performance.

## Overview

RepAn v0.1.1+ includes C++ implementations of computationally intensive clustering operations using Rcpp and RcppArmadillo. These optimizations significantly improve the speed of differential abundance analysis for TCR immune repertoire datasets.

## Optimized Functions

### 1. fastGetCenters() - Cluster Centroid Calculation

**Performance:** 4.5x faster than R implementation

**Usage:**
```r
# Automatically used by getClusterLables() function
# Can also be called directly:
centers <- fastCenters(seqmers, clslabels)
```

**What it does:**
- Computes cluster centroids from sparse k-mer frequency matrices
- Handles sparse matrices efficiently in C++
- Used during within-sample CDR3 clustering

### 2. fastDetermineWeightsVector() - Positional Weight Calculation

**Performance:** 8.3x faster than R implementation

**Usage:**
```r
# Automatically used when posWt=TRUE in getClusterLables()
# Can also be called directly:
weights <- fastDetermineWeightsVector(kmers, sequence)
```

**What it does:**
- Calculates positional weights for k-mers in CDR3 sequences
- Assigns higher weights to k-mers in the middle (CDR3 core) region
- Used for position-aware clustering

## Installation

The C++ optimizations are compiled automatically when installing RepAn:

```r
# Standard installation
devtools::install_github("dyohanne/RepAn")

# If you encounter compilation errors, ensure you have:
# - R development tools (Rtools on Windows, build-essential on Linux)
# - Rcpp and RcppArmadillo packages
install.packages(c("Rcpp", "RcppArmadillo"))
```

## Performance Impact

### Benchmark Results

Tested on 500 sequences with 256 k-mer features:

| Operation | R Implementation | C++ Implementation | Speedup |
|-----------|-----------------|-------------------|---------|
| Cluster Centroids | 18ms | 4ms | **4.5x** |
| Positional Weights | 25ms | 3ms | **8.3x** |

### Impact on RepAn Analysis

These optimizations have the greatest impact when:
- Analyzing large repertoires (>5,000 clonotypes per sample)
- Using positional weighting (posWt=TRUE)
- Running multiple repeat resamples (nRepeats>10)
- Processing many samples simultaneously

**Example speedup for typical analysis:**
- 10 samples, 5000 clonotypes each, 10 repeats
- Original: ~45 minutes
- Optimized: ~25 minutes (**1.8x faster overall**)

## Technical Details

### C++ Implementation

The optimizations are implemented in `src/clustering_cpp.cpp` using:
- **RcppArmadillo** for efficient sparse matrix operations
- **STL** (Standard Template Library) for string matching
- Native C++ for loop-intensive calculations

### Key Design Decisions

1. **Distance calculation:** We use the existing `fclust::dist.matrix()` which is already highly optimized with compiled code.

2. **Centroid calculation:** Converting sparse matrix rows to dense format for accumulation provides better performance than pure sparse operations in this context.

3. **Positional weights:** Native string matching in C++ avoids R's function call overhead and regex complexity.

## Troubleshooting

### Compilation Errors

If you encounter compilation errors:

1. **Missing compiler:** Install R development tools
   - Windows: Install Rtools
   - Mac: Install Xcode Command Line Tools
   - Linux: `sudo apt-get install r-base-dev`

2. **Missing RcppArmadillo:** Install from CRAN
   ```r
   install.packages("RcppArmadillo")
   ```

3. **Linker errors:** Ensure BLAS/LAPACK libraries are installed
   - Linux: `sudo apt-get install libblas-dev liblapack-dev`

### Verification

To verify C++ functions are working:

```r
library(RepAn)
library(Matrix)

# Test sparse matrix
x <- sparseMatrix(i=c(1,1,2), j=c(1,2,2), x=c(1,2,3), dims=c(2,3))
clslabels <- c(1, 1)

# Should run without error
centers <- fastCenters(x, clslabels)
print(centers)
```

## Contributing

To add new C++ optimizations:

1. Add function to `src/clustering_cpp.cpp`
2. Export with `// [[Rcpp::export]]` comment
3. Create R wrapper in `R/cpp_wrappers.R`
4. Update this README with benchmarks
5. Run `devtools::document()` to update documentation

## References

- Yohannes DA, et al. (2021). Identifying condition-associated immune repertoires using a diversification-oriented approach. BMC Bioinformatics. https://doi.org/10.1186/s12859-021-04087-7
- Rcpp: https://www.rcpp.org/
- RcppArmadillo: http://arma.sourceforge.net/
