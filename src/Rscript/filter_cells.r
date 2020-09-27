#` Filtering of cells with below a given number of non-zero features
#` @param bmat            (matrix)   A saparse matrix with rows of features and columns of cells
#` @param features_number (interger) Filter the cells if they have less than this number (default:100)
#` @return                (matrix)   A filtered sparse matrix
filter_cells = function(bmat, feature_number=100) {
  bmat = bmat[, Matrix::colSums(bmat) >= feature_number]
  return(bmat)
}
