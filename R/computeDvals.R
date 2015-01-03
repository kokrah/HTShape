#' Compute d-values.
#' 
#' @param t3 a vector containing L-skew estimates for each gene.
#' @param t4 a vector containing L-kurt estimates for each gene.
#' @param plot indicating whether to show intermediary plots (Default=FALSE).
#' @param span the lowess smoother span
#' @export
computeDvals <- function (t3, t4, plot=FALSE, span=0.5) {
  # 0. Remove any NaNs (genes with constant expression level)
  nans <- is.nan(t3)
  
  if ( any(nans) ) {
    
    warning("Removing NaN values!")
    t3 <- t3[!nans]
    t4 <- t4[!nans]
  }
  
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
  if ( plot ) {
    
    oldpar <- par(mgp=c(1.5, 0.5, 0), 
                  mar=c(2.5, 2.5, 1, 0.5), 
                  mfrow=c(2,2))
    cex <- 0.4
    
    # Plot 1 (skewness adj. kurt. plot)
    Lab.palette <- colorRampPalette(c("gray80", "gray65", "gray40", "gray30"), 
                                    space = "Lab")
    data.pts.col <- densCols(shape.centered[, 1],
                             shape.centered[, 2], 
                             colramp=Lab.palette)
    
    plot(shape.centered[, 1], shape.centered[, 2],  
         cex = cex, 
         pch = 16,
         xlab = expression(paste("centered ", tau[3])),
         ylab = expression(paste("centered adjusted ", tau[4])),
         main = "Skew adjusted kurtosis",
         col = data.pts.col,
         cex.main = 0.9)
    abline(v=0, h=0, col="tomato")

    # Plot 2 (skewness volcano plot)
    plot(shape.centered[, 1], 
         -log10(dvals), pch = ".",
         xlab = expression(paste("centered ", tau[3])),
         ylab = expression(-log[10](dvals)),
         col = "gray40",
         cex.main = 0.9)
    abline(v=0, col="tomato")
    
    # Plot 3 (hist. of squared distances)
    hist(D, breaks = 40, freq = FALSE, xlim=c(0, 20),
         main = expression(Distribution~of~squared~distance~D),
         cex.main = 0.9)
    curve(dchisq(x, 2), 0.05, 20, add = TRUE, lty=2)  
    
    # Plot 4 (hist. of d-values)
    hist(dvals, freq = FALSE,
         main = "Histogram of d-values",
         cex.main = 0.9)

    par(oldpar)
  }
  
  # 5. Return d-values
  dvals              
}
