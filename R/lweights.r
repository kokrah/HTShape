#' Compute weights for each expression based on l-scale
#' 
#' @param exprs normalized expression on log scale
#' @param groups character vector indicating the group membership of each sample
#' @param span span for lowess (default=0.5)
#' @param plot show fit ? (default=FALSE)
#' @export
lweights = function (exprs, groups, span=0.5, plot=FALSE) {
  
  uGroups = unique(groups)
  
  l2 = c()
  
  for (g in uGroups) {
    l2 = cbind(l2, fitShape(exprs[,g==groups], 2)$lmoms[, "l2"])
  }
  
  y = rowMeans(l2)
  x = rowMeans(exprs)
  
  l = lowess(x, y, f=span)
  f = approxfun(l, rule=2)
  
  if (plot) {
    
    oldpar = par(mar=c(4, 4, 1.5, 0.5))
    
    plot(x, y, pch=".", col="gray", xlab="Average expression", ylab="L-scale",
         main="L-scale weights")
    
    points(x, f(x), pch=".", col="red")
    
    par(oldpar)
    
  }
  
  W = 1 / apply(exprs, 2, f)
  
  colnames(W) = colnames(exprs)
  rownames(W) = rownames(exprs)
  
  W
}