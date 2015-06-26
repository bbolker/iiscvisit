rmvn.spa <- function(x, y, p, method = "exp",
                     vsqrt=NULL,
                     ret.vsqrt=FALSE,
                     useMASS=FALSE,
                     empirical=FALSE) {
#######################################################################
  ## Function to generate spatially
  ## autocorrelated random normal variates using 
  ## the eigendecomposition method. Spatial covariance can follow either an 
  ## exponential or a Gaussian model.

  ## REQUIRED ARGUMENTS
  ##x         vector of length n representing the x coordinates
  ##y         vector of length n representing the y coordinates
  ##p         the range of the spatial models
  ##method    either "exp" (exponential) or "gaus" (gaussian)
  ##vsqrt     pre-computed SVD of covariance matrix
  ##ret.vsqrt return SVD of covariance matrix for re-use?
  ##
  ##VALUE
  ## a vector of spatially correlated random normal variates
  ## with zero mean and unit variance is returned
#########################################################################
  imeth <- charmatch(method, c("exp", "gaus"), nomatch = NA)
  if(is.na(imeth)) stop("method should be \"exp\", or \"gaus\"") 
  ch.fn0 <-function(vmat, tol = 1e-007){
    vs <- svd(vmat)
    vsqrt <- t(vs$v %*% (t(vs$u) * sqrt(vs$d)))
    vsqrt
  }
  ch.fn <- function(n, mu, vmat, tol = 1e-007, vsqrt=NULL,
                    useMASS=FALSE){
    if (!useMASS) {
      p <- ncol(vmat)
      if (is.null(vsqrt))
        vsqrt <- ch.fn0(vmat,tol)
      ans <- matrix(rnorm(n * p), nrow = n) %*% vsqrt
      ans <- sweep(ans, 2, mu, "+")
      dimnames(ans) <- list(NULL, dimnames(vmat)[[2]])
      ans
    } else {
      mvrnorm(n,mu,Sigma=vmat,empirical=empirical)
    }
  }
  n <- length(x)
  z <- matrix(0, ncol = 1, nrow = n)
  tmpmat <- cbind(x + 0, y + 0)
  xy <- tmpmat[, 1:2] + 0
  dmat <- sqrt(outer(x,x, "-")^2+outer(y,y,"-")^2)
  if(imeth == 1) {
    covmat <- exp( - dmat/p)
  }
  else {
    if(imeth == 2) {
      covmat <- exp( - (dmat/p)^2)
    }
  }
  if (ret.vsqrt)
    return(ch.fn0(vmat=covmat))
  else {
    z <- ch.fn(1, mu = rep(0, n), vmat = covmat, vsqrt=vsqrt,
               useMASS=useMASS)
    if (!useMASS) t(z)[, 1]
    else z
  }
}

