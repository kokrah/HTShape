#' Perform Manova on shape (t3, t4) of samples (Wilk's Exact Manova).
#' 
#' @param sfit output from fitShape()
#' @param groups a character vector indicating sample group membership
#' @param nPerm number of permutations for non-parametric test 
#' @export
shapeManova = function (sfit, groups, nPerm=500) {
  Y = t(sfit$lmoms)
  
  # Make design matrix
  X = model.matrix(~ 0 + factor(groups))
  
  # Compute group means
  B = solve(t(X) %*% X) %*% t(X) %*% Y
  
  # Predicted Y values
  Yhat = X %*% B
  
  # Compute within SS matrix
  resids = Y - Yhat
  SSwithin = t(resids) %*% resids
  
  # Compute total SS matrix
  SStotal = t(Y) - colMeans(Y)
  SStotal = SStotal %*% t(SStotal)
  
  # Compute Wilk's lambda
  WL = det(SSwithin) / det(SStotal)

  # Compute random WLs
  randWL = c()
    
  for (k in 1:nPerm) {
    # Make design matrix
    groups = sample(groups, length(groups), replace=FALSE)
      
    X = model.matrix(~ 0 + factor(groups))
      
    # Compute group means
    B = solve(t(X) %*% X) %*% t(X) %*% Y
      
    # Predicted Y values
    Yhat = X %*% B
      
    # Compute within SS matrix
    resids = Y - Yhat
    SSwithin = t(resids) %*% resids
      
    # Compute total SS matrix
    SStotal = t(Y) - colMeans(Y)
    SStotal = SStotal %*% t(SStotal)
      
    # Compute Wilk's lambda
    randWL = c(det(SSwithin) / det(SStotal), randWL)
  }
    
  # Transform WL into Fstat
  g = ncol(X)
  n = nrow(X)
  Fstat =  ((n - g - 1) / (g - 1)) * ((1 - sqrt(WL)) / sqrt(WL))
  randFstat = ((n - g - 1) / (g - 1)) * ((1 - sqrt(randWL)) / sqrt(randWL))
  
  # Compute p-value
  mean(randFstat > Fstat)
}

