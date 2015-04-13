#' Filter out genes with low counts
#' 
#' @param counts raw counts matrix
#' @param lib.size library size (defualt = NULL)
#' @param thresh cpm threshold (defualt =1)
#' @param minSamples number of samples required to exceed threshold (default=2)
#' @export
filterCounts = function (counts, lib.size = NULL, thresh = 1, minSamples = 2) {
  cpms <- 2^log2CPM(counts, lib.size = lib.size)$y
  keep <- rowSums(cpms > thresh) >= minSamples
  counts <- counts[keep, ]
  counts
}