#' Compute Lmoments (Code lifted from Lmoments package).
#' 
#' @param data matrix or dataframe
#' @param rmax number of lmoments
Lmoments = function (data, rmax) {
  
  na.rm = FALSE 
  returnobject = FALSE 
  trim = c(0, 0)
  
  if (!identical(trim, c(0, 0)) & !identical(trim, c(1, 1))) {
    stop("The current version of Lmoments supports only ordinary L-moments (trim=c(0,0)) and T1L-moments (trim=c(1,1))")
  }
  if (identical(trim, c(1, 1)) & rmax > 4) {
    warning("The current version of t1lmoments uses rmax=4.")
    rmax <- 4
  }
  data <- as.matrix(data)
  p <- dim(data)[2]
  if (!na.rm) {
    if (identical(trim, c(0, 0))) 
      L <- Lmoments_calc(data, rmax)
    if (identical(trim, c(1, 1))) 
      L <- t1lmoments(data, rmax)
  }
  else {
    L <- array(, c(p, rmax))
    for (i in 1:p) {
      xi <- data[, i]
      xi <- xi[!is.na(xi)]
      if (identical(trim, c(0, 0))) 
        L[i, ] <- Lmoments_calc(xi, rmax)
      if (identical(trim, c(1, 1))) 
        L[i, ] <- t1lmoments(xi, rmax)
    }
  }
  colnames(L) <- paste("L", 1:rmax, sep = "")
  if (p > 1) 
    rownames(L) <- names(data)
  if (returnobject) {
    Lobject <- list()
    Lobject$lambdas <- L
    Lobject$trim <- trim
    if (identical(trim, c(0, 0))) 
      Lobject$source <- "Lmoments"
    if (identical(trim, c(1, 1))) 
      Lobject$source <- "t1lmoments"
    if (rmax > 2) {
      ratios <- L[, 3:rmax]/(cbind(L[, 2]) %*% rep(1, (rmax - 2)))
      colnames(ratios) <- paste("tau", 3:rmax, sep = "")
      if (p == 1) 
        ratios <- cbind(rbind(L[, 1:2]), cbind(ratios))
      if (p > 1) 
        ratios <- cbind(L[, 1:2], ratios)
      if (p > 1) 
        rownames(ratios) <- names(data)
      Lobject$ratios <- ratios
    }
    return(Lobject)
  }
  else {
    return(L)
  }
}


#' Compute Lmoments (Code lifted from Lmoments package).
#' 
#' @param data matrix or dataframe
#' @param rmax number of lmoments
Lmoments_calc = function (data, rmax = 4) 
{
  data <- as.matrix(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  x <- array(, c(p, n))
  L <- array(, c(p, rmax))
  for (i in 1:p) {
    x[i, ] <- sort(data[, i])
  }
  if (rmax == 1) 
    return(rowMeans(x))
  bcoef <- array(, c(rmax, n))
  bcoefm <- array(, c(rmax, p, n))
  b <- array(, c(p, rmax))
  bcoef[1, ] <- seq(0, 1, by = (1/(n - 1)))
  bcoefm[1, , ] <- t(array(rep(bcoef[1, ], p), c(n, p)))
  b[, 1] <- rowMeans(x)
  b[, 2] <- rowMeans(bcoefm[1, , ] * x)
  L[, 1] = b[, 1]
  if (rmax > 2) {
    for (r in 2:(rmax - 1)) {
      rr <- r + 1
      bcoef[r, ] <- bcoef[r - 1, ] * seq((-(r - 1)/(n - 
                                                      r)), 1, by = (1/(n - r)))
      bcoefm[r, , ] <- t(array(rep(bcoef[r, ], p), c(n, 
                                                     p)))
      b[, rr] <- rowMeans(bcoefm[r, , ] * x)
    }
  }
  for (r in 1:(rmax - 1)) {
    L[, r + 1] <- 0
    for (k in 0:r) {
      kk <- k + 1
      L[, r + 1] <- L[, r + 1] + (-1)^(r - k) * gamma(r + 
                                                        k + 1)/(gamma(k + 1)^2)/gamma(r - k + 1) * b[, 
                                                                                                     kk]
    }
  }
  return(L)
}
