#' Perform Manova on shape (t3, t4) of samples (Wilk's Exact Manova).
#' 
#' @param lrats a (2 by number of samples) matrix containing (t3, t3) for each sample.
#' @param groups a character vector indicating sample group membership.
#' @param plot plot L-skew and L-kurt estimates
#' @param groupCol col for each sample
#' @export
shapeManova = function (exprs, groups, plot=TRUE, groupCol=NULL) {
  
  Y = fitShape(t(exprs))$lrats
  
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
    
    if (is.null(groupCol)) {
      
      pch = 0:(length(uGroups) - 1) + 6     
      
      plot(Y, pch=pch, ylab="L-kurt", xlab="L=skew",
           main=paste0("Shape Manova (Wilk's P-value: ", format(round(pval, 4), nsmall=4), ")"))
      
      legend("topleft", legend=uGroups, pch=pch)
      
    }else{
      
      plot(Y, pch=19, ylab="L-kurt", xlab="L-skew", cex=0.8, col=groupCol,
           main=paste0("Shape Manova (Wilk's P-value: ", format(round(pval, 4), nsmall=4), ")"))
      
      legend("topleft", legend=uGroups, fill=unique(groupCol))
      
    }
    
    par(oldpar)

  }
  
  list(WL=WL, Fstat=Fstat, df1=df1, df2=df2, pval=pval)
}
