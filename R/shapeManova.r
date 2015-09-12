#' Perform Manova on shape (t3, t4) of samples (Wilk's Exact Manova).
#' 
#' @param sfit sfit
#' @param groups a character vector indicating sample group membership.
#' @param plot plot L-skew and L-kurt estimates
#' @param groupCol col for each sample
#' @param loc.scal loc.scal
#' @export
shapeManova = function (sfit, groups, plot=FALSE, loc.scal=FALSE, groupCol=NULL) {
  
  if (loc.scal) {
    Y = t(sfit$lrats)
  }else{
    Y = t(sfit$lmoms)
  }
  
  # Make design matrix
  uGroups = unique(groups)
  X = model.matrix(~ 0 + factor(groups, levels=uGroups))
  
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

  # Transform WL into Fstat
  g = ncol(X)
  n = nrow(X)
  Fstat =  ((n - g - 1) / (g - 1)) * ((1 - sqrt(WL)) / sqrt(WL))
  
  # Compute p-value
  df1 = 2 * (g - 1)
  df2 = 2 * (n - g - 1)
  pval = df(Fstat, df1=df1, df2=df2)
  
  if (plot) {
    
    oldpar = par(mar=c(4, 4, 1.5, 0.5))
    plot(0, 0)
    par(oldpar)

  }
  
  list(WL=WL, Fstat=Fstat, df1=df1, df2=df2, pval=pval)
}
