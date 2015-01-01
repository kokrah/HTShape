#' Find volatile genes.
#' 
#' @param t3 A vector containing L-skew estimates for each gene.
#' @param t4 A vector containing L-kurt estimates for each gene.
#' @param plots Indicating whether to show intermediary plots (Default=FALSE).
#' @return d-values for each gene.
#' @export
computeDvals <- function (t3, t4, plots=FALSE, span=0.5) {
  
  # 1. Adjust L-kurt (t4).
  lowess.fit <- lowess(t3, t4, f=span)
  approx.lowess.fit <- approxfun(lowess.fit, rule=2)
  adj.t4 <- t4 - approx.lowess.fit(t3)
  
  # 2. Compute squared distances (D).  
  shape <- cbind(t3, adj.t4)
    
  # compute the inverse covariance matrix
  A <- cov(shape) 
  A.inv <- solve(A)
    
  # center t3 and adj.t4
  mu <- colMeans(shape)
  shape.centered <- t(t(shape) - mu)
  
  # squared distance
  D <- shape.centered %*% A.inv
  D <- D * shape.centered
  D <- rowSums(D)
  
  # 3. Compute d-values
  dvals <- 1 - pchisq(D, df=2)
  
  # 4. Plot summary 
  if ( plots ) {
    oldpar <- par(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0.5))
    
    Lab.palette <- colorRampPalette(c("gray80", "gray65", "gray30", "gray20"), 
                                    space = "Lab")
    data.pts.col <- densCols(t3, t4, colramp=Lab.palette)
    cex <- 0.4
    
    # Plot 1 (plain SO-plot)
    plot(t3, t4,  
         cex = cex, 
         pch = 16,
         xlab = expression(tau[3]),
         ylab = expression(tau[4]),
         main = expression(Kurt-Skew~Trend),
         col = data.pts.col)
    lines(lowess.fit, col = "gray90", lwd=2)
    
    # Plot 2 (skewness adj. kurt. plot)
    plot(t3, adj.t4,  
         cex = cex, 
         pch = 16,
         xlab = expression(tau[3]),
         ylab = expression(paste("adjusted ", tau[4])),
         main = "Skew adjusted kurtosis",
         col = data.pts.col)
    abline(v=0, col="gray")
    sel1 <- dvals <= 1/10^4 
    sel2 <- dvals <= 1/10^3 & !sel1
    
    # show outliers 
    points(t3[sel1], adj.t4[sel1], pch=4) 
    points(t3[sel2], adj.t4[sel2], pch=6)
    
    # Plot 3 (skewness volcano plot)
    plot(t3, -log10(dvals), cex=0.3)
    points(t3[sel1], -log10(dvals)[sel1], pch=4)
    points(t3[sel2], -log10(dvals)[sel2], pch=6)
    abline(h=c(3, 4), v=0, col="gray")
    
    # Plot 4 (hist. of squared distances)
    hist(D, breaks = 40, freq = F, xlim=c(0, 20),
         main = expression(Distribution~of~sqaured~distance~D))
    curve(dchisq(x, 2), 0.05, 20, add=T, lty=2)  
    
    par(oldpar)
  }
  
  # 5. Return d-values
  dvals              
}
