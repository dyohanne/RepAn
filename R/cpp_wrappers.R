#' Fast distance matrix calculation using C++
#'
#' @param x A sparse matrix (dgCMatrix) or regular matrix
#' @param method Distance method: "euclidean" or "cosine"
#' @return A distance matrix
#' @export
fastDistMatrix <- function(x, method = "euclidean") {
  # Convert to sparse matrix if not already
  if (!inherits(x, "sparseMatrix")) {
    x <- as(x, "sparseMatrix")
  }
  
  # Call appropriate C++ function
  if (method == "euclidean") {
    dist_mat <- fastEuclideanDist(x)
  } else if (method == "cosine") {
    dist_mat <- fastCosineDist(x)
  } else {
    stop("Unsupported distance method. Use 'euclidean' or 'cosine'")
  }
  
  # Convert to dist object
  rownames(dist_mat) <- rownames(x)
  colnames(dist_mat) <- rownames(x)
  
  return(as.dist(dist_mat))
}

#' Fast cluster centroid calculation using C++
#'
#' @param seqmers A sparse matrix of sequence k-mer frequencies
#' @param clslabels Integer vector of cluster labels
#' @return A matrix of cluster centroids
#' @export
fastCenters <- function(seqmers, clslabels) {
  # Ensure sparse matrix format
  if (!inherits(seqmers, "sparseMatrix")) {
    seqmers <- as(seqmers, "sparseMatrix")
  }
  
  # Call C++ function
  centers <- fastGetCenters(seqmers, as.integer(clslabels))
  
  # Set column names
  colnames(centers) <- colnames(seqmers)
  
  return(centers)
}
