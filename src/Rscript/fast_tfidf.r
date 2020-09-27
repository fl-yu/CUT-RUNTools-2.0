#` A fast version of row scaling of sparse TF matrix by IDF vector, this function wiil be invoked by function tfidf
#` @param tf   (sparse matrix)   A saparse matrix of term frequency matrix
#` @param idf  (vector)          A vector of inverse document frequency 
#` @return     (sparse matrix)   A sparse matrix of TF-IDF 
fast_tfidf = function(tf, idf) {
   tf <- t(tf)
   tf@x <- tf@x * rep.int(idf, diff(tf@p))
   tf <- t(tf)
   return(tf)
}
