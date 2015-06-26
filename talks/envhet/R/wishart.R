
## utility function for adding up the outer products of a set of vectors
sumouter <- function(alpha) {
  matrix(apply(apply(alpha,1,function(x)outer(x,x)),1,sum),nrow=ncol(alpha))
}
## would cov(x)*n be more efficient?

## recipe from Gelman et al. p. 481
## rwishart <- function(n=1,nu,Sigma) {
##  require(MASS)  ## for mvrnorm()
##  if (nu<dim(Sigma)[1]) stop("can't simulate for nu<k")
##  ## check for square, pos. def. Sigma?
##  alpha <- mvrnorm(nu,mu=rep(0,nrow(Sigma)),Sigma)
##  return(sumouter(alpha))
## }

## from Bill Venables, via S-PLUS help list, via google:
## "The following code should generate a single p x p Wishart matrix
## realisation with df degrees of freedom. If SqrtSigma is not
## specified the parent variance matrix is assumed to be the
## identity matrix; if it is specified it must be some matrix
## "square root" of the variance matrix, Sigma, e.g. chol(Sigma).
## In other words crossprod(SqrtSigma) must be equal to Sigma. (As
## this is a big-ticket computation I figured it was better to do it
## once outside the function than inside every time.)
## Limited testing suggests it is not completely stupid, at any
## rate, but that's all. Caveat emptor."

rwishart <- function(df, p = nrow(SqrtSigma), SqrtSigma = diag(p)) {
  if ((Ident <- missing(SqrtSigma)) && missing(p))
    stop("either p or SqrtSigma must be specified")
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, df:(df-p+1)))
  if(p > 1) {
    pseq <- 1:(p-1)
    Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
  }
  if(Ident)
    crossprod(Z)
  else
    crossprod(Z %*% SqrtSigma)
}

## according to Gelman et al. p. 81:
## to resample from the variance-covariance matrix
## S = \sum_{i=1}^n (y_i-\bar y) (y_i-\bar y)^T
## Samples from the joint posterior distribution of (\mu,\Sigma)
## are easily obtained using the following procedure:
## first, draw \Sigma|y \~ Inv-Wishart_{\nu_n}(\Lambda_n^{-1}),
## then draw \mu|\Sigma,y \~ N(\mu_n, \Sigma_\kappa_n)

## with multiv. Jeffreys prior density, p(\mu,\Sigma) \propto |\Sigma|^{-(d+1)/2}
## the corresponding posterior distribution can be written as
## \Sigma|y \~ Inv-Wishart^{n-1}(S)

batch <- TRUE
if (!batch) {
  Sigma <- diag(1:3)
  cholzinv <- chol(solve(Sigma))  ##
  nt <- 3
  rw <- try(solve(rwishart(nt-1,SqrtSigma=cholzinv),tol=1e-10))
}

## parametric bootstrap of a set of space-time data:
## *MUST ASSUME* that covariance matrix is stationary,
## generates data with mean=0
paramcovboot <- function(x,n,scale=FALSE,randnorm=TRUE,
                         use.cov="pairwise.complete.obs",tol=1e-6) {
  require(MASS)  ## for mvrnorm()
  ## x: a space (columns) x time (rows) data matrix
  ## n: number of bootstrap reps to return
  ## n[1] is number of new VC samples, n[2] is number
  ## of (multivariate) sample per VC sample
  ## if length(n)==1  n-> c(n,1)
  ## scale: standardize each row (separately) to mean=0, sd=1?
  ## randomize: randomize at level of points/cov matrices?
  npts <- ncol(x)
  nt <- nrow(x)
  #### restriction only needed for wishart
  ##  if (nt <= npts+1)
  ##    stop("not enough time points")
  #### next 4 lines replaced by cov(x)
  ##  if (scale)  ## standardize to N(0,1)
  ##    scale(x,TRUE,TRUE)
  ##  xsum <- sumouter(x)
  ##  rw <- xsum/nt  ## empirical estimate of var-cov matrix
  if (scale) vc <- cor(x) else vc <- cov(x,use=use.cov)
  #### needed for rwishart
  ##  cholzinv <- chol(solve(xsum))  ##
  zsqrt <- matsqrt(vc)
  if (length(n)==1) n <- c(n,nt)
  if (n[2]<nt && !randnorm)
    stop("can't generate mvrnorm with exact (empirical) covariances for n[2]<nt")
  ret <- array(dim=c(n[2],npts,n[1]))
  for (i in 1:n[1]) {
### most commented stuff is crap from an attempt to use
### nt==npts+1 with rwishart, which doesn't work
      ##      ok <- FALSE; tryctr <- 0
      ##      while (!ok && tryctr<100) {
      ##        tryctr <- tryctr+1
### 
      ## rw <- try(solve(rwishart(nt-1,SqrtSigma=cholzinv),tol=1e-10))
      rw <- newvc(n=nt,sqrtSigma=zsqrt,tol=tol)
      ##        ok <- TRUE
      ##        cr <- class(rw)
      ##        if (!is.null(cr) && cr=="try-error") {
      ##          ok <- FALSE
      ##          browser()
      ##        }
      ##        }
    ## can't (doh!) use empirical=TRUE for n=1 in mvrnorm()
    ## leads to more levels-of-randomization questions
    ## choices: pick new v-c matrix or not (in this case,
    ##            it would be silly not to)
    ##   pick new v-c matrix for each rep? (yes)
    ##   pick mv normal with exact v-c matrix?
    ##   (can only do this if we're picking more than one
    ##   time step at a time)
      ret[,,i] <- mvrnorm(n[2],
                         mu=rep(0,npts),Sigma=rw,
                         empirical=!randnorm)
    }
  ret <- drop(ret)
  ret
}

paramcovboot.old <- function(x,n,scale=FALSE,randnorm=TRUE) {
  require(MASS)  ## for mvrnorm()
  ## x: a space (columns) x time (rows) data matrix
  ## n: number of bootstrap reps to return
  ## n[1] is number of new VC samples, n[2] is number
  ## of (multivariate) sample per VC sample
  ## if length(n)==1  n-> c(n,1)
  ## scale: standardize each row (separately) to mean=0, sd=1?
  ## randomize: randomize at level of points/cov matrices?
  npts <- ncol(x)
  nt <- nrow(x)
  #### restriction only needed for wishart
  if (nt <= npts+1)
    stop("not enough time points")
  if (scale) vc <- cor(x) else vc <- cov(x)
  cholzinv <- chol(solve(vc))  ##
  if (length(n)==1) n <- c(n,nt)
  if (n[2]<nt && !randnorm)
    stop("can't generate mvrnorm with exact (empirical) covariances for n[2]<nt")
  ret <- array(dim=c(n[2],npts,n[1]))
  for (i in 1:n[1]) {
    rw <- try(solve(rwishart(nt-1,SqrtSigma=cholzinv),tol=1e-10))
      ret[,,i] <- mvrnorm(n[2],
                         mu=rep(0,npts),Sigma=rw,
                         empirical=!randnorm)
  }
  ret <- drop(ret)
  ret
}

## copied from mvrnorm() in MASS library
matsqrt <- function(Sigma,tol=1e-6) {
  eS <- eigen(Sigma, sym = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1]))) 
    stop("Sigma is not positive definite")
  return(t(eS$vectors %*% diag(sqrt(pmax(ev, 0)))))
}

batch <- TRUE
if (!batch) {
  Sigma <- cov(X$z[,1:4])
  Sigma <- cov(X$z)
  m <- matsqrt(Sigma)
  m0 <- chol(Sigma)
  S2 <- t(m) %*% m
  max(abs(S2-Sigma))
}


newvc <- function(Sigma=NULL,n,sqrtSigma=NULL,tol=1e-6) {
  ## algorithm from J. Gentle (Springer, 1998), _Random Number
  ## Generation and Monte Carlo Methods_ p. 107, Algorithm 3.7:
  ## Hartley and Harris (1963) and Odell and Feiveson (1966)
  ## also cited
  if (is.null(sqrtSigma)) {
    sqrtSigma <- matsqrt(Sigma,tol=tol)
  }
  d <- nrow(sqrtSigma)
  nmat <- matrix(0,nrow=d,ncol=d)
  nmat[upper.tri(nmat,diag=TRUE)] <- c(0,rnorm((d*(d+1)/2-1)))
  cvec <- rchisq(d,df=n-1)
  B <- matrix(nrow=d,ncol=d)
  B[1,1] <- cvec[1]
  tmpmat <- matrix(0,nrow=d,ncol=d)
  tmpmat[upper.tri(tmpmat)] <- 1
  nmat2 <- nmat
  diag(nmat2) <- 0
  ## y_j + sum_{i=1}^{j-1} z_{ij}^2 for j=2,3,...,d
  diag(B)[2:d] <- cvec[2:d]+apply((nmat2^2)[,-1,drop=FALSE],2,sum)
  B[1,-1] <- nmat[1,-1]*sqrt(cvec[1])
  ## b_{ij} = z_{ij}*sqrt{y_i} + sum_{k=1}^{i-1} z_{ki} z_{kj}
  ## for i<j = 2,3,...,d
  cp <- crossprod(nmat2)
  if (d>2) {
    for (i in 2:(d-1))
      B[i,(i+1):d] <- cp[(i+1):d,i]+nmat[i,(i+1):d]*sqrt(cvec[i])
  }
  ## symmetrize
  for (i in 1:(d-1))
    B[(i+1):d,i] <- B[i,(i+1):d]
  1/n*(t(sqrtSigma) %*% B %*% sqrtSigma)
}

if (!batch) {
  Sigma <- matrix(c(2,1,0.5,1.5),nrow=2)
  library(MASS)
  x <- mvrnorm(50,mu=c(0,0),Sigma)
  plot(x)
  y <- paramcovboot(x,50)
}
