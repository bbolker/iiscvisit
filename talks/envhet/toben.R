##############################################################################################
mSncf<-function(x, y, z, df = 25, type = "boot", resamp = 0,
                npoints = 1024, save = FALSE, 
                max.it=25, xmax=FALSE, na.omit = FALSE,
                noisy = FALSE){
##############################################################################################
#mSncf is the abridged function to estimate the nonparametric covariance function
#(using a #smoothing spline as an equivalent kernel) 
#
#The function requires multiple observations at each location (use spline.correlog
#otherwise). 
#
#In R use library(modreg)
#
#REQUIRED ARGUMENTS
#x         vector of length n representing the x coordinates 
#y         vector of length n representing the y coordinates 
#z         matrix of dimension n x p representing p observation at each location
#
#df        degrees of freedom for the spline
#type      takes the value "boot" to generate a bootstrap distribution or "null" to generate a 
#             null distribution for the estimator under randomization
#resamp    is the number of resamples for the bootstrap or the null distribution
#npoints   is the number of points at which to save the value for the spline function (and
#             confidence envelope / null distribution)
#save      if True, the whole matrix of output from the resampling is saved (a resamp x npoints
#             dimensional matrix)
#xmax	   if FALSE the max observed in the data is used. Otherwise all distances greater
#	      than xmax is omitted
#na.omit   if TRUE, missing values is accommodated through a pairwise deletion.
#
#VALUE
#an object of class Sncf is returned consisted of the following components:
#real      $predicted$x is the x coordinates for the fitted covariance function
#          $predicted$y is the y values for the covariance function
#          $cbar is the regional average synchrony
#boot      gives the analogous output from the bootstrap or randomization resampling
#boot$summary  gives the full vector of output for cbar
#              and a quantile summary for the resampling distribution
#boot$boot     if save=T, the full raw matrices from the resampling is saved
############################################################################################
require(stats)  ## loads modreg (for smooth.spline) if necessary
#the following sets up the output:
  real <- list(cbar = NA, predicted = list(x = matrix(NA, nrow = 1, ncol = npoints), 
                            y = matrix(NA, nrow = 1, ncol = npoints)))

  inner.na <- function(x, y){
    ##this function is used instead of corr to do pairwise deletion
    ## in the presence of missing values. 
    ok <- (is.finite(x) & is.finite(y))
    x <- x[ok]
    y <- y[ok]
    n2 <- length(x)
    if(n2 <= 2)
      return(NA)
    coef <- cor(x, y)
    if(is.na(coef))
      return(NA)
    z <- coef
    return(z)
  }
  ##end of inner.na
  NAO <- F
  ##check for missing values
  if(any(!is.finite(z))) {
    if(na.omit){
      warning("Missing values exist; Pairwise deletion will be used")
      NAO <- T
    }
    else {
      stop("Missing values exist; use na.omit = T for pairwise deletion")
    }
  }
  
  ##This generates the moran distances
  ##the odd adding of zero is just to ensure that all vectors 
  ##are treated as numeric
  n <- dim(z)[1]
  p <- dim(z)[2]
  z <- as.matrix(z)+0
  if(!NAO){
    ##this is the correlation in the absence of NA's
    moran <- cor(t(z))
  }
  else{
    ##this part (in the presence of missing values) is SLOW (it's a double loop)
    ##it is probably possible to recode as some sort of outer product; but haven't
    ##had time to do it!
    moran <- matrix(1, nrow = n, ncol = n)
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        tmp <- inner.na(as.vector(z[i,]), as.vector(z[j,]))
        moran[j, i] <- tmp
        moran[i, j] <- tmp
      }
    }
  }

##then generating geographic distances
  xdist <- sqrt(outer(x,x, "-")^2+outer(y,y,"-")^2)
  maxdist <- ifelse(!xmax, max(na.omit(xdist)), xmax)

##The spline function
  triang <- lower.tri(xdist)

  #####
  ## break at this point, return xdist[triang], moran[triang]?
	
  u <- xdist[triang]
  v <- moran[triang]
  sel <- is.finite(v)
  u <- u[sel]
  v <- v[sel]
  v <- v[u <= maxdist]
  u <- u[u <= maxdist]
  sel <- is.finite(u)
  u <- u[sel]
  v <- v[sel]

  real$cbar <- mean(v)
  sobj <- smooth.spline(u, v, df = df)

  xpoints <- seq(0, maxdist, length = npoints)

  lx <- predict(sobj, x = xpoints)
  ly <- 1:length(lx$y)
  ## linear interpolation of location of zero (?)
  choise <- ly[lx$y < 0][1]
  pos <- lx$x[choise - 1]
  neg <- lx$x[choise]
  pos <- pos + (neg - pos)/2
##   tmp <- smooth.spline(lx)  ## unused?

  real$predicted <- list(x = xpoints, y = lx$y)

##End of spline fit

  boot <- list(NULL)
  boot$boot.summary <- list(NULL)
  if(resamp != 0) {
    ##here is the bootstrapping/randomization
    boot$boot.summary$cbar <- matrix(NA, nrow = resamp, ncol = 1)
    predicted <- list(x = matrix(NA, nrow = 1, ncol = npoints), y = matrix(NA, nrow = resamp, ncol = npoints))
    predicted$x[1,] <- xpoints
    type <- charmatch(type, c("boot", "perm"), 
                      nomatch = NA)
    if(is.na(type))
      stop("method should be \"boot\", or \"perm\"")
    
    for(i in 1:resamp) {
      if (noisy) cat(i, " of ", resamp, "\n")
      if(type == 1) {
        trekkx <- sample(1:n, replace = T)
        trekky <- trekkx
      }
      
      if(type == 2) {
        trekky <- sample(1:n, replace = F)
        trekkx <- 1:n
      }
      
      xdistb <- xdist[trekkx, trekkx]
      
      triang <- lower.tri(xdist)
      
      xdistb <- xdistb[triang]
      moranb <- moran[trekky, trekky][triang]	

      if(type == 1) {
        moranb <- moranb[!(xdistb == 0)]
        xdistb <- xdistb[!(xdistb == 0)]
      }
      
      u <- xdistb
      v <- moranb
      sel <- is.finite(v)
      u <- u[sel]
      v <- v[sel]
      v <- v[u <= maxdist]
      u <- u[u <= maxdist]
      sel <- is.finite(u)
      u <- u[sel]
      v <- v[sel]

      boot$boot.summary$cbar[i,1] <- mean(v)
      sobj <- smooth.spline(u, v, df = df)
      lx <- predict(sobj, x = xpoints)

      ly <- 1:length(lx$y)
      choise <- ly[lx$y < 0][1]
      pos <- lx$x[choise - 1]
      neg <- lx$x[choise]
      pos <- pos + (neg - pos)/2
##      tmp <- smooth.spline(lx) ## unused?
      predicted$y[i,] <- lx$y
    }
##end of bootstrap loop!

    if(save == T) {
      boot$boot <- list(predicted = predicted)
			}
    else {
      boot$boot <- NULL
    }

    ty <- apply(predicted$y, 2, quantile, probs = c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 
                                            0.75, 0.9, 0.95, 0.975, 1), na.rm = T)
    dimnames(ty) <- list(c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1), NULL)
    tx <- predicted$x
    boot$boot.summary$predicted <- list(x = tx,y = ty)
  }

##The following else is if resamp=0

  else {
          boot <- NULL
          boot.summary <- NULL
        }

  res <- list(real = real, boot = boot, max.distance = maxdist)
  class(res) <- "mSncf"
  res
}
##############################################################################
plot.mSncf <- function(obj, xmax = 0, text = TRUE, add=FALSE,
                       ylim = c(-1,1), xlab="Distance",
                       ylab="Correlation", ...) {
#############################################################################
#this is the generic plot function for Sncf objects
#############################################################################
  xmax <- ifelse(xmax == 0, obj$max.distance, xmax)
  cbar <- round(obj$real$cbar, 2)
  if(!is.null(obj$boot$boot.summary)){
    rul <- round(quantile(obj$boot$boot.summary$cbar[,1],
                          probs = c(0.025, 0.975), na.rm=T), 2)
  }
  if(!add){
    plot(obj$real$predict$x, obj$real$predict$y, xlim = c(0, xmax),
         ylim = ylim, type = "l", xlab = xlab, ylab = ylab ,...)
    if(text && !is.null(obj$boot$boot.summary)){
      ttext <- paste("Regional synch:", cbar,
                     " {", rul[1], ", ", rul[2], "}", sep = "")
      title(main=ttext)
    }
  }
  lines(obj$real$predict$x, obj$real$predict$y,...)
#  lines(c(0, max(obj$real$predict$x)), c(0, 0),...)
#  lines(c(0, max(obj$real$predict$x)), c(cbar, cbar),...)
  abline(h=c(0,cbar),lty=2,...)
  if(!is.null(obj$boot$boot.summary)){
    lines(obj$boot$boot.summary$predict$x,
          obj$boot$boot.summary$predict$y["0.025", ],...)
    lines(obj$boot$boot.summary$predict$x,
          obj$boot$boot.summary$predict$y["0.975", ],...)}
  
}



