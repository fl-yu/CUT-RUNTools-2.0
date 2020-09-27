#` TF-IDF normalization with the binary matrix (sparse matrix): calculation of LSI (latent semantic indexing)
#` @param bmat          (sparse matrix)  A saparse matrix of binarized cell-by-feature matrix
#` @param TF            (logic)          If calculate the term frequencies of the input matrix. Use either raw counts or divide by total counts in each cell (default:T)
#` @param log_TF        (logic)          If the input/TF matrix will be log scaled (default:T)
#` @param scale_factor  (numberic)       Multiply terms in TF matrix by scale_factor prior to log1p. Equivalent to adding small pseudocount but doesn't cast to dense matrix at any point 
#` @return              (sparse matrix)  A sparse matrix of TF-IDF 
tfidf <- function(bmat, TF=TRUE, log_TF=TRUE, scale_factor=100000) {
  if (TF) {
    # "term frequency" method
    tf <- t(t(bmat) / Matrix::colSums(bmat))
  } else {
    # "raw count" method
    tf <- bmat
  }
  # Either TF method can optionally be log scaled
  if (log_TF) {
    if (TF) {
      tf@x <- log1p(tf@x * scale_factor)
    } else {
      tf@x <- log1p(tf@x * 1)
    }
  }
  # IDF "inverse document frequency smooth" method
  idf <- log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  # TF-IDF
  tf_idf_counts <- fast_tfidf(tf, idf)
  rownames(tf_idf_counts) <- rownames(bmat)
  colnames(tf_idf_counts) <- colnames(bmat)
  return(tf_idf_counts)
}
