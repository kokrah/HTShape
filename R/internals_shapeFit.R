#' Weight Coefficient for Sample L-moment (Internal).
#' 
#' @param r The r-th sample l-moment.
#' @param j The j-th term.
#' @param n The sample size.
#' @return The j-th weight factor of the r-th sample l-moment.
weight.factor <- function (r, j, n) {
  
  upper.limit <- min(j - 1, r - 1)
  i <- 0:upper.limit
  
  factor1 <- (-1)^(r - 1 - i)
  factor2 <- choose(r - 1, i)
  factor3 <- choose(r - 1 + i, i)
  factor4 <- choose(j - 1, i)
  factor5 <- choose(n - 1, i)
  
  sum(factor1 * factor2 * factor3 * factor4 / factor5)
  
}

#' Weight Coefficients for Sample L-moment (Internal).
#' 
#' @param r The r-th sample l-moment.
#' @param n The sample size.
#' @return A vector of weight factors.
weightFactors <- function (r, n) {
  
  sapply(1:n, weight.factor, r = r, n = n)
  
}

#' Weight Coefficients for Sample L-moments.
#' 
#' @param nlmom The number of sample l-moments.
#' @param n The sample size.
#' @return A matrix of weight factors.
#' @export
weightMatrix <- function (nlmom, n) {
  
  if (n < 1) {
    stop("n must be greater than or equal to 1.")
  }
  
  if (nlmom < 1) {
    stop("nlmom must be greater than or equal to 1.")
  }
  
  if (nlmom > n) {
    stop("nlmom must be less than or equal to n.")
  }
  
  res <- sapply(1:nlmom, weightFactors, n = n)
  colnames(res) <- paste0("W", 1:nlmom)
  
  res
  
}
