# Normalise ATAC data by binarisation, followed by TF-IDF normalisation, 
# function with thanks to Ricard Argelaguet:

# mtx is the binary peak matrix with dim (npeaks,ncells)
#'  method=1: The TF-IDF implementation used by Stuart & Butler et al. 2019. This computes \eqn{\log(TF \times IDF)}.
#'  method=2: The TF-IDF implementation used by Cusanovich & Hill et al. 2018. This computes \eqn{TF \times (\log(IDF))}.
#'  method=3: The log-TF method used by Andrew Hill. This computes \eqn{\log(TF) \times \log(IDF)}.
tfidf <- function(mtx, method = 1, scale.factor = 1e4) {
  
  require(Matrix)
  
  npeaks <- colSums(mtx)
  if (any(npeaks == 0)) {
    warning("Some cells contain 0 total counts")
  }
  tf <- tcrossprod(mtx, y = Diagonal(x=1/npeaks))
  rsums <- rowSums(mtx)
  if (any(rsums == 0)) {
    warning("Some features contain 0 total counts")
  }
  idf <- ncol(mtx) / rsums
  if (method == 2) {
    idf <- log(1 + idf)
  } else if (method == 3) {
    tf <- log1p(tf * scale.factor)
    idf <- log(1 + idf)
  }
  mtx.tfidf <- Diagonal(n = length(idf), x = idf) %*% tf
  if (method == 1) {
    mtx.tfidf <- log1p(mtx.tfidf * scale.factor)
  }
  colnames(mtx.tfidf) <- colnames(mtx)
  rownames(mtx.tfidf) <- rownames(mtx)
  # set NA values to 0
  mtx.tfidf[is.na(mtx.tfidf)] <- 0
  return(mtx.tfidf)
}
