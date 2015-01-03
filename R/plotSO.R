#' The Symmetry Outlier Plot.
#' 
#' @param t3 a vector containing L-skew estimates for each gene.
#' @param t4 a vector containing L-kurt estimates for each gene.
#' @param dataName character vector indicating name of dataset.
#' @param medianCol color of L-skew and L-kurt median lines (Default = "red")
#' @param defaultLook determines whether the default appearance of the (L-skew, L-kurt) points
#' is used (Default=TRUE).
#' @param ... arguments to control the (L-skew, L-kurt) points if defaultLook = FALSE.
#' @export
plotSO <- function(t3, t4,
                   dataName = "Data name",
                   medianCol = "gray90",
                   defaultLook = TRUE, ...){

  # Make a blank plot
  old.par <- par(mar=c(3.5, 2.5, 1, 0.5))
  plot(c(-1, 1.18), c(-0.8, 1),
       pch = "",
       xaxt = "none",
       yaxt = "none",
       xlab = "",
       ylab = "",
       main = "Symmetry-Outlier Plot (SO-Plot)",
       cex.main = 0.8)
  
  if ( defaultLook ) {
    
    # Define colors for points based on density
    Lab.palette <- colorRampPalette(c("gray80", "gray65", "gray40", "gray30"), 
                                    space = "Lab")
    data.pts.col <- densCols(t3, t4, colramp=Lab.palette)

    # Plots L-skew and L-kurt values
    points(t3, t4, 
           col = data.pts.col, 
           cex = 0.3,
           pch = 16) 
    
  }else{
    
    points(t3, t4, ...)
    
  }
  
  abline(v = c(-1, 1), col = "gray60")
  
  # Plot extreme borders
  points(seq(-1, 1, 0.01), 0.25 * (5 * seq(-1, 1, 0.01)^2 - 1),
         type = "l")
  points(c(-1, 1), c(1, 1), type = "l")
  
  # Add t4 boxplot
  boxplot(t4, add = TRUE, horizontal = FALSE, 
          at = 1.14, yaxt = "none", boxwex = 0.3, 
          cex = 0.8)
  
  # Add t3 boxplot
  boxplot(t3, add = TRUE, horizontal = TRUE, 
          at = -0.75, xaxt = "none", boxwex = 0.3, 
          cex = 0.8)
  
  # Add lines and text indicating moderate to extreme t3 estimates
  abline(v = c(-.35, -.2, -.05, .05, .2, .35), 
         col = "gray40", 
         lty = c(1, 4, 2, 2, 4, 1),
         lwd = 0.5)
  
  text(c(-.35, -.2, -.05) + .02, rep(.75, 3),
       rev(c("minor", "moderate", "large")), 
       srt = 90, cex = 0.6, pos = 1)
  text(-.6, .75, "extreme/\nvolatile", cex=0.6)
  
  # Make axis for t3
  axis(1, seq(-1, 1, .05), F, tck = 0.01)
  tmp1 <- c(-1, -.75, -.5, -.35, -.2, -.05, .05, .2, .35, .5, .75, 1) 
  mtext(tmp1, at = tmp1, side = 1, line = 0.2, 
        cex = 0.6, las = 3)
  mtext(expression(L-skew~(tau[3])), 1, line=1.5, cex = 0.6)
  mtext(dataName, 1, line = 2.2, cex = 0.6)
  
  # Make axis for t4
  axis(2, seq(-0.8, 1, 0.05), F, tck = 0.01)
  tmp2 <- seq(-0.8, 1, 0.2)
  mtext(tmp2, at = tmp2, side = 2, line = 0.2, 
        cex = 0.6, las=2)
  mtext(expression(L-kurt~(tau[4])), 2, line = 1.5, cex = 0.6)
  
  # Add t3 and t4 median lines
  medt3 <- median(t3)
  medt4 <- median(t4)
  abline(v = medt3, col = medianCol)
  abline(h = medt4, col = medianCol)
  points(medt3, medt4, col = medianCol)
  
  # plot gh-curves
  quad.h.0 <- function(x) {
    
    t1 <-  1.216105e-01 - 8.662085e-06  * x + 8.111085e-01 * x^2 
    t2 <-  5.229546e-05 * x^3 - 4.484714e-02 * x^4 - 5.125358e-05 * x^5 + 1.105605e-01 * x^6
    t1 + t2
  }
  
  quad.h.25 <- function(x) {
    
    t1 <-  2.754634e-01 + 2.144862e-06 * x + 6.201384e-01 * x^2 
    t2 <- 1.571900e-05 * x^3 -5.646926e-02 * x^4 - 2.357846e-05 * x^5 + 1.597199e-01 * x^6
    t1 + t2
  }
  
  quad.h.5 <- function(x) {
    
    t1 <- 4.626708e-01 + 1.941385e-06 * x + 4.555605e-01 * x^2 
    t2 <- 4.812860e-06 * x^3 -1.255911e-01 * x^4 - 1.260371e-05 * x^5 + 2.064426e-01 * x^6
    t1 + t2
  }
  
  x <- seq(-1, 1, 0.01)
  points(x, quad.h.0(x), type = "l")
  points(x, quad.h.25(x), type = "l")
  points(x, quad.h.5(x), type = "l", lty=2)
  
  # Show Gaussian an Uniform distr. lines
  abline(h = 0.1226, v = 0, col = "gray20")
  abline(h = 0, v = 0, col = "gray20", lty = 2)
  points(0, 0.1226, bg = "white", pch = 21)
  
  # Print l3 summary (0.25, 0.5, 0.75)-quantile
  t3.summary <- round(quantile(t3, probs=c(0.25, 0.5, 0.75)), 2)
  t3.summary <- paste0(dataName, 
                       " L-skew: (25%, 50%, 75%) = (",
                       t3.summary[1], ", ",
                       t3.summary[2], ", ",
                       t3.summary[3], ")")
  
  if ( TRUE ) {
  
    print(t3.summary)
    
  }
  
  par(old.par)
}