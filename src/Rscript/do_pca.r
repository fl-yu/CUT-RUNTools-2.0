#` Perform fast PCA (use irlba package) on LSI sparse matrix, retaining observation names (no further scaling or centering done)
#` @param mat  (sparse matrix)  A saparse matrix of TF-IDF matrix (LSI) used for PCA 
#` @param dims (integer)        number of PCs to calculate (default:10)
#` @return     (sparse matrix)  A sparse matrix of PCA
do_pca <- function(mat, dims=10) {
  pca.results <- irlba(t(mat), nv=dims, fastpath=FALSE)
  PCA_result <- pca.results$u %*% diag(pca.results$d)
  rownames(PCA_result) <- colnames(mat)
  colnames(PCA_result) <- paste0('PC_', 1:dims)
  return(PCA_result)
}
