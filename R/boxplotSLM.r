#' Boxplot of sample L-moments (boxplotSLM)
#' 
#' @param data data matrix
#' @param groups a character vector indicating sample group membership
#' @param data.name character vector indicating name of dataset.
#' @param show.points show.points
#' @param group.col group.col
#' @param ... arguments passed to boxplots
#' @export
#' 
boxplotSLM = function(data, groups, data.name=NULL, show.points=TRUE, group.col=NULL, ...) {
  fit = fitShape(data)
  LM = fit$lmoms  
  
  op = par(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 2, 0.5), mfrow=c(2, 2))
  
  # L1
  if (is.null(data.name)) {
    main = "L1"    
  }else{
    main = paste0(data.name, ": L1")
  }
  
  boxplot(LM["L1",] ~ groups, main=main, ...)
  
  if (show.points) {
    if (is.null(group.col)) {
      group.col = factor(groups)
    }
    
    points(jitter(as.numeric(factor(groups)), amount=0.25), LM["L1",],
           pch=19, cex=0.7, col=group.col)
  }
  
  # L2
  if (is.null(data.name)) {
    main = "L2"    
  }else{
    main = paste0(data.name, ": L2")
  }
  
  boxplot(LM["L2",] ~ groups, main=main, ...)
  
  if (show.points) {
    if (is.null(group.col)) {
      group.col = factor(groups)
    }
    
    points(jitter(as.numeric(factor(groups)), amount=0.25), LM["L2",],
           pch=19, cex=0.7, col=group.col)
  }
  
  # L3
  if (is.null(data.name)) {
    main = "L3"    
  }else{
    main = paste0(data.name, ": L3")
  }
  
  boxplot(LM["L3",] ~ groups, main=main, ...)
  
  if (show.points) {
    if (is.null(group.col)) {
      group.col = factor(groups)
    }
    
    points(jitter(as.numeric(factor(groups)), amount=0.25), LM["L3",],
           pch=19, cex=0.7, col=group.col)
  }
  
  # L4
  if (is.null(data.name)) {
    main = "L4"    
  }else{
    main = paste0(data.name, ": L4")
  }
  
  boxplot(LM["L4",] ~ groups, main=main, ...)
  
  if (show.points) {
    if (is.null(group.col)) {
      group.col = factor(groups)
    }
    
    points(jitter(as.numeric(factor(groups)), amount=0.25), LM["L4",],
           pch=19, cex=0.7, col=group.col)
  }
  
  par(op)
}