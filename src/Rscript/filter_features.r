#` Filtering of features detected by a given number of cells
#` @param bmat            (matrix)   A saparse matrix with rows of features and columns of cells
#` @param features_number (interger) Filter features if they were detected less than this number of cells (default:20)
#` @return                (matrix)   A filtered sparse matrix
filter_features = function(bmat, cell_number=20) {
  bmat = bmat[Matrix::rowSums(bmat) >= cell_number, ]
  return(bmat)
}
