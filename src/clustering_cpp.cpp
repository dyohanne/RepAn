// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Fast Euclidean distance calculation for sparse matrices
// [[Rcpp::export]]
NumericMatrix fastEuclideanDist(const arma::sp_mat& x) {
  int n = x.n_rows;
  NumericMatrix dist(n, n);
  
  // Precompute squared norms for all rows
  NumericVector row_norms_sq(n);
  for (int i = 0; i < n; i++) {
    row_norms_sq[i] = arma::accu(x.row(i) % x.row(i));
  }
  
  // Compute pairwise Euclidean distances using ||a-b||^2 = ||a||^2 + ||b||^2 - 2*a·b
  for (int i = 0; i < n; i++) {
    dist(i, i) = 0.0;
    for (int j = i + 1; j < n; j++) {
      double dot_product = arma::accu(x.row(i) % x.row(j));
      double d_sq = row_norms_sq[i] + row_norms_sq[j] - 2.0 * dot_product;
      // Ensure non-negative due to numerical errors
      double d = sqrt(std::max(0.0, d_sq));
      dist(i, j) = d;
      dist(j, i) = d;
    }
  }
  
  return dist;
}

// Fast Cosine distance calculation for sparse matrices
// [[Rcpp::export]]
NumericMatrix fastCosineDist(const arma::sp_mat& x) {
  int n = x.n_rows;
  NumericMatrix dist(n, n);
  
  // Precompute norms for all rows
  NumericVector norms(n);
  for (int i = 0; i < n; i++) {
    norms[i] = sqrt(arma::accu(x.row(i) % x.row(i)));
  }
  
  // Compute pairwise cosine distances
  for (int i = 0; i < n; i++) {
    dist(i, i) = 0.0;
    for (int j = i + 1; j < n; j++) {
      double dot_product = arma::accu(x.row(i) % x.row(j));
      double similarity = dot_product / (norms[i] * norms[j] + 1e-10);
      double d = 1.0 - similarity;
      dist(i, j) = d;
      dist(j, i) = d;
    }
  }
  
  return dist;
}

// Fast computation of cluster centroids for sparse matrices
// [[Rcpp::export]]
arma::mat fastGetCenters(const arma::sp_mat& seqmers, const IntegerVector& clslabels) {
  int n_features = seqmers.n_cols;
  
  // Find unique cluster labels (excluding 0)
  std::set<int> unique_labels_set;
  for (int i = 0; i < clslabels.size(); i++) {
    if (clslabels[i] > 0) {
      unique_labels_set.insert(clslabels[i]);
    }
  }
  
  std::vector<int> unique_labels(unique_labels_set.begin(), unique_labels_set.end());
  int n_clusters = unique_labels.size();
  
  arma::mat centers(n_clusters, n_features);
  centers.zeros();
  
  // Compute centroid for each cluster
  for (int c = 0; c < n_clusters; c++) {
    int label = unique_labels[c];
    int count = 0;
    
    for (int i = 0; i < clslabels.size(); i++) {
      if (clslabels[i] == label) {
        // Convert sparse row to dense row for accumulation
        arma::rowvec dense_row = arma::conv_to<arma::rowvec>::from(arma::mat(seqmers.row(i)));
        centers.row(c) += dense_row;
        count++;
      }
    }
    
    if (count > 0) {
      centers.row(c) /= count;
    }
  }
  
  return centers;
}

// Fast positional weight calculation
// [[Rcpp::export]]
double fastDetermineWeight(const std::string& kmer, const std::string& seq) {
  double total_weight = 0.0;
  int seq_length = seq.length();
  int kmer_length = kmer.length();
  
  if (kmer_length > seq_length) {
    return 0.0;
  }
  
  // Define weight groups based on sequence thirds
  double third = seq_length / 3.0;
  double group1 = third;
  double group2 = 2.0 * third;
  
  // Find all occurrences of kmer in seq
  for (int i = 0; i <= seq_length - kmer_length; i++) {
    bool match = true;
    for (int j = 0; j < kmer_length; j++) {
      if (seq[i + j] != kmer[j]) {
        match = false;
        break;
      }
    }
    
    if (match) {
      // Position i is 0-indexed, convert to 1-indexed for weight calculation
      int pos = i + 1;
      double weight;
      
      if (pos <= group1) {
        weight = 5.0;  // v-area weight
      } else if (pos <= group2) {
        weight = 10.0;  // middle area weight
      } else {
        weight = 1.0;  // j-area weight
      }
      
      total_weight += weight;
    }
  }
  
  return total_weight;
}

// Vectorized version of positional weight calculation for multiple kmers
// [[Rcpp::export]]
NumericVector fastDetermineWeightsVector(const CharacterVector& kmers, const std::string& seq) {
  int n_kmers = kmers.size();
  NumericVector weights(n_kmers);
  
  for (int i = 0; i < n_kmers; i++) {
    std::string kmer = Rcpp::as<std::string>(kmers[i]);
    weights[i] = fastDetermineWeight(kmer, seq);
  }
  
  return weights;
}

// Compute weighted kmer matrix for a set of sequences
// [[Rcpp::export]]
arma::mat fastWeightedKmerMatrix(const arma::sp_mat& seq_mers, 
                                  const CharacterVector& seqs,
                                  const CharacterVector& kmers) {
  int n_seqs = seqs.size();
  int n_kmers = kmers.size();
  arma::mat weighted_mers(n_seqs, n_kmers);
  
  for (int i = 0; i < n_seqs; i++) {
    std::string seq = Rcpp::as<std::string>(seqs[i]);
    NumericVector weights = fastDetermineWeightsVector(kmers, seq);
    
    for (int j = 0; j < n_kmers; j++) {
      weighted_mers(i, j) = seq_mers(i, j) * weights[j];
    }
  }
  
  return weighted_mers;
}
