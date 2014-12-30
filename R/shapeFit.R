#' Sample L-moments ratios.
#' 
#' This function computes the sample l-moments ratios for a given
#' dataset (vector / matrix / data.frame).
#' @param x data vector / matrix.
#' @param nLmom The number of l-moments to compute (Default = 4).
#' @param samplesByCols Do columns represent samples ? (Default = TRUE)
#' @export
shapeFit <- function (x, nLmom = 4, samplesByCols = TRUE) {
  
  # 1. If x is a vector transform it to a single column matrix
  if ( is.vector(x) ) {
    
    x <- as.matrix(x) # the rows are the samples
    
  }else{
    
    if ( samplesByCols ) {
      x <- t(x) # put genes in the colums and samples in the rows  
    }
    
  }
  
  # 2. Check for NAs (Replace NAs with median)
  xNA <- is.na(x)
  
  if ( any(xNA) ) {
    
    warning("Replacing NA values with median!")
    
    # A. How many NAs are in each column
    numberOfNAs <- colSums(xNA)
    
    # B. Grab columns that contain NAs
    naCols <- numberOfNAs > 0
    
    # C. Compute the median of column(s) that contain NAs 
    xSub <- x[, naCols] # xSub is vector if only 1 column contains NAs
    xSub <- as.matrix(xSub) # make xSub a column matrix in case it is a vector
        
    # D. Replace NAs with column medians 
    colMedians <- apply(xSub, 2, median, na.rm=TRUE)
    x[xNA] <- rep(colMedians, numberOfNAs[naCols])
  
  }
  
  # 4. Compute weight matrix
  n <- nrow(x) 
  W <- weightMatrix(nLmom, n)
  
  # 5. Compute L-moments
  xs <- apply(x, 2, sort)
  lmoms <- (t(W) %*% xs) / n
  lmoms <- t(lmoms)
  colnames(lmoms) <- paste0("l", 1:nLmom)
  
  # 6. Modify results for identical columns genes
  id <- lmoms[, 2] < .Machine$double.eps ^ 0.5
  
  if ( any(id) ) {
    
    if ( samplesByCols ) {
      
      warning(paste0(sum(id), " out of ", length(id), " row(s) are identical!"))
      
    }else{
      
      warning(paste0(sum(id), " out of ", length(id), " column(s) are identical!"))
      
    }
    
    lmoms[id, 2:nLmom] <- 0
  }
  
  # 7. Compute L-ratios if nLmom > 2
  if ( nLmom > 2 ) {
    
    lrats <- lmoms[, 3:nLmom] / lmoms[, 2]
        
    if ( is.vector(lrats) ) {
      
      lrats <- matrix(lrats, nrow = 1)
      
    }
    
    colnames(lrats) <- paste0("t", 3:nLmom)

  }else{
    
    lrats <- NULL
    
  }

  
  # 8. Results
  list(lrats = lrats, lmoms = lmoms, W = W)
  
}
