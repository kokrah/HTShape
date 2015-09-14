#' Sample L-moments ratios.
#' 
#' This function computes the sample l-moments ratios for a given
#' dataset (vector / matrix / data.frame).
#' @param x data vector / matrix.
#' @param nLmom the number of l-moments to compute (Default = 4).
#' @param samplesByCols do columns represent samples ? (Default = TRUE)
#' @export
fitShape <- function (x, nLmom=4) {

  if (nLmom < 2) {
    stop("nLmom must be 2+ (K. Okrah)")
  }
  
  # Compute L-moments
  lmoms = Lmoments(x, rmax=nLmom)
  colnames(lmoms) = paste0("Lmom-", 1:nLmom)
  rownames(lmoms) = colnames(x)
  
  # Compute L-ratios
  lrats = lmoms / lmoms[, "Lmom-2"]
  colnames(lrats) = paste0("Lrat-", 1:nLmom)     
  
  # Compute L-CV
  lcv <- 1 / lrats[, "Lrat-1"]
  names(lcv) = rownames(lrats)
  
  if (nLmom == 2) {
    lrats = NULL
  }else{
    lrats = t(lrats)[-(1:2),]
  }
  
  # Results
  list(lcv = lcv, lrats = lrats, lmoms = t(lmoms))
  
}
