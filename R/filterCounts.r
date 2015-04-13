#' Filter out genes with low counts
#' 
#' @param counts raw counts matrix
#' @param lib.size library size (defualt = NULL)
#' @param thresh cpm threshold (defualt =1)
#' @param minSamples number of samples required to exceed threshold (default=2)
#' @export
filterCounts = function (counts, lib.size = NULL, thresh = 1, minSamples = 2) {
  
  log2CPM = function (qcounts, lib.size = NULL) {
    if (is.null(lib.size)) 
    lib.size <- colSums(qcounts)
    y <- t(log2(t(qcounts + 0.5)/(lib.size + 1) * 1e+06))
    return(list(y = y, lib.size = lib.size))
  }
  
  cpms <- 2^log2CPM(counts, lib.size = lib.size)$y
  keep <- rowSums(cpms > thresh) >= minSamples
  counts <- counts[keep, ]
  counts
}