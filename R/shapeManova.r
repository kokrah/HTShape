#' Perform Manova on shape (t3, t4) of samples (Wilk's Exact Manova).
#' 
#' @param sfit output from fitShape()
#' @param groups a character vector indicating sample group membership
#' @param nPerm number of permutations for non-parametric test 
#' @param plot plot
#' @param lrats lrats
#' @param groupCol group color
#' @export
shapeManova = function (data, groups, nPerm=500, lrats=FALSE, plot=FALSE, groupCol=NULL) {
  
  sfit = fitShape(data, nLmom=4)
  
  if (lrats) {
    Y = t(sfit$lrats)  
  }else{
    Y = t(sfit$lmoms)
  }
  
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
  pval = mean(randFstat > Fstat)
  
  if (plot) {
  
    pal = colorRampPalette(c("red", "orange", "white", "steelblue3", "navy"))(n=99)
    
    main = paste0("Shape Manova (nPerm=", nPerm,")\n P-value = ", round(pval, 3))
    
    gplots::heatmap.2(t(Y), trace="none", Rowv=F, Colv=F, dendrogram="none", col=pal,
                      ColSideColors=groupCol, density.info="none", scale="row",
                      colsep=cumsum(table(groups)), main=main)
    
    legend("bottomleft", legend=unique(groups), pch=19, col=unique(groupCol), 
           title="Groups")
  }
  
  pvals = apply(Y, 2, function(x) summary(aov(x ~ groups))[[1]][1,"Pr(>F)"])
  pvals = p.adjust(pvals, method="bonferroni")
  
  list(apvals = round(pvals,3), manova.pval = round(pval, 3))
  
}

