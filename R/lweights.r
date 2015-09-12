#' Compute weights for each expression based on l-scale
#' 
#' @param exprs normalized expression on log scale
#' @param groups character vector indicating the group membership of each sample
#' @param span span for lowess (default=0.5)
#' @param plot show fit ? (default=FALSE)
#' @export
lweights = function (exprs, groups, span=0.5, plot=FALSE) {
  
  # Compute residuals
  Y = exprs
  X = model.matrix(~ groups)
  B = solve(t(X) %*% X) %*% t(X) %*% t(Y)
  Yhat = t(X %*% B)
  resids = Y - Yhat
  
  x = rowMeans(Y)
  y = fitShape(resids, 2)$lmoms["Lmom-2", ]
  
  l = lowess(x, y, f=span)
  f = approxfun(l, rule=2)
  
  if (plot) {
    
    oldpar = par(mar=c(4, 4, 1.5, 0.5))
    
    plot(x, y, pch=".", col="gray40", xlab="Average expression", ylab="L-scale of residuals",
         main="L-scale weights")
    
    points(x, f(x), pch=".", col="red")
    
    par(oldpar)
    
  }
  
  W = 1 / apply(exprs, 2, f)
  
  colnames(W) = colnames(exprs)
  rownames(W) = rownames(exprs)
  
  W
}