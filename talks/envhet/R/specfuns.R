## TO DO: merge/harmonize spec.hankel and spec.hankel2

## returns (slow) Hankel transform
## [faster coded in C? or not, because already pretty vectorized?]
hankel <- function(x,rvec,rhovec,debug=FALSE,
                   ang.freq=FALSE,
                   norm="none",
                   weight=FALSE,
                   spec=FALSE) {
  if (ang.freq) freqfac <- 1  else freqfac <- 2*pi
  inorm <- charmatch(norm,c("none","norm","integrate","intnorm"),nomatch=0)
  if (inorm==0) {
    warning("unknown normalization type, using 'none'")
    inorm <- 1
  }
  x <- x[order(rvec)]
  if (any(is.na(rvec))) stop("can't handle NAs in rvec")
  rvec <- sort(rvec)
  ruvec <- unique(rvec)
  wtvec <- rep(1,length(x))
  if (any(duplicated(rvec))) {
    ## account for duplicated values
    ruvec.w <- apply(outer(rvec,ruvec,"=="),2,sum) ## number of duplicates
    wtvec <- 1/ruvec.w
    wtvec <- wtvec[match(rvec,ruvec)]  ## expand to match
  }
  if (weight) {
    ## figure out the "domain" of each point
    ## (bar rule, could use trapezoid/Simpson's rule?)
    rbreak <- c(0,(ruvec[-1]+ruvec[-length(ruvec)])/2,max(rvec))
    wtvec2 <- diff(rbreak)
    wtvec2 <- wtvec2[match(rvec,ruvec)] ## expand
    wtvec  <- wtvec*wtvec2
  }
  wtvec <- wtvec/sum(wtvec)*length(rvec)  ## normalize weight vector
  bmat <- besselJ(outer(freqfac*rhovec,rvec),0)
  n.out <- length(rhovec)
  ret <- as.vector(bmat %*% (rvec*x*wtvec))
  dr <- diff(range(rvec))/length(rvec)
  drho <- diff(range(rhovec))/(n.out)
  if (spec) ret <- ret^2
  if (debug) cat("dr: ",dr,"\n")
  if (debug) cat("drho: ",drho,"\n")
  normval <- switch(inorm,
                    1,                ## no normalization
                    1/sum(ret),       ## sum(x)==1
                    dr,               ## integrate: multiply by dr
                    1/(sum(ret)*drho))  ## integral of transform=sum(x*d(rho))=1
  if (debug) cat("normval: ",normval,"\n")
  ret <- ret*normval
  ret
}

## wrapper function for hankel: figures out r vector,
## optionally scales, squares the transform (to get power
## spectrum), smooths
spec.hankel2 <- function(x,y,z,freq="auto",
                         fmax=NULL,length.out=NULL,
                         smooth=FALSE,smooth.df=NULL,
                         ang.freq=FALSE,demean=FALSE,
                         norm="none",
                         spec=FALSE,weight=FALSE,
                         debug=FALSE) {
  x <- as.vector(x); y <- as.vector(y); z <- as.vector(z)
  if (demean) z <- z-mean(z)
  rvec <- sqrt(x^2+y^2)
  if (freq=="auto") {
    f <- calc.freq(x,length.out=length.out,fmax=fmax,ang.freq=ang.freq)
    ## should this be calc.freq(rvec) ??
  } else{
    f <- freq
  }
  f <- f[f>0] ## eliminate zero ordinate
  h <- hankel(z,rvec,f,ang.freq=ang.freq,norm=norm,weight=weight,
              spec=spec,
              debug=debug)
  ##  if (spec) h <- h^2
  if (smooth) {
    h <- smooth.spec(f,h,smooth.df)
  }
  ret <- list(freq=f,spec=h)
  class(ret) <- "spec.hankel"
  ret
}

calc.freq <- function(x,length.out=NULL,fmax=NULL,ang.freq=FALSE,
                      debug=FALSE) {
  dr <- mean(diff(unique(sort(x))))
  r.rng <- diff(range(x))
  if (is.null(fmax)) fmax <- 1/(2*dr)
  if (is.null(length.out)) length.out <- round(r.rng/dr)
  if (debug) {
    cat("average dr=",dr,"\nr range=",r.rng,
        "\nlength.out=",length.out,"\nfmax=",fmax,"\n",sep="")
  }
  f <- seq(0,fmax,length=length.out+1)  ## will drop 0 frequency shortly
  if (ang.freq) f <- f*(2*pi)
  f
}

smooth.spec <- function(f,h,smooth.df=NULL,smooth.n=NULL,type="bs") {
  if (type=="bs") {
    require(splines)
    if (is.null(smooth.df))
      smooth.df <- round(sqrt(length(h)))
    tmp <- data.frame(h=h,f=f)
    tmps <- try(glm(h~bs(f,df=smooth.df),
                    family=quasi("log","mu^2"),data=tmp,maxit=100))
    ca <- class(tmps)
    if (!is.null(ca) && ca=="try-error") {
      h <- rep(NA,length(tmp$f))
    } else {
      h <- predict(tmps,data.frame(f=tmp$f),type="response")
    }
    ret <- h
  } else if (type=="smooth.spline") {
      ## require(modreg) ## smooth.spline has moved to stats
    sobj <- smooth.spline(f, h, df = smooth.df)
    xpoints <- seq(0, max(f), length = smooth.n)
    p <- predict(sobj, x = xpoints)
    h <- p$y
    f <- p$x
    ret <- list(f=f,h=h)
  }
  ret
}

## hankel transform on a 1-D vector
spec.hankel <- function(x,rvec,rhovec=NULL,
                        fmax=NULL,length.out=NULL,weight=FALSE,
                        smooth=FALSE,smooth.n=NULL,smooth.df=25,
                        spec=TRUE,norm="none",
                        debug=FALSE,
                        ang.freq=FALSE) {
  if (is.null(rhovec))
    rhovec <- calc.freq(rvec,fmax=fmax,length.out=length.out,debug=debug)
  if (debug) {
    cat("spec.hankel: ",length(rvec)," input, ",
        length(rhovec)," output, ",sep="")
    if (!smooth)
      cat("no smooth\n")
    else
      cat(smooth.n,"/",smooth.df," smooth pts/df\n",sep="")
  }
  h <- hankel(x,rvec,rhovec,ang.freq=ang.freq,norm=norm,weight=weight,
              spec=spec)
  ##  if (spec) h <- h^2
  ## smoothing copied from Bjornstad mSncf
  if (smooth) {
    sm <- smooth.spec(rhovec,h,smooth.df=smooth.df,smooth.n=smooth.n,
                      type="smooth.spline")
    h <- sm$h
    rhovec <- sm$f
  }
  ret <- list(freq=rhovec,spec=h)
  class(ret) <- "spec.hankel"
  ret
}

## "envelope" function: x is a matrix with columns=replicates
env <- function(x,ctrfun="median",q=c(0,1),MARGIN=1,na.rm=TRUE) {
  if (ctrfun=="none") {
    r <- t(apply(x,MARGIN,quantile,q,na.rm=na.rm))
  } else if (ctrfun=="median") {
    r <- t(apply(x,MARGIN,quantile,c(q[1],0.5,q[2]),na.rm=na.rm))
 } else if (ctrfun=="mean") {
    tmp <- t(apply(x,MARGIN,quantile,q,na.rm=na.rm))
    r <- cbind(tmp[,1],apply(x,MARGIN,mean),tmp[,2])
  }
  r
}


## bivariate exponential (missing 1/(2*pi))
biexp <- function(r,rate=1) {
  rate^2*exp(-rate*r)
}

## bivariate Gaussian (missing 1/(2*pi))
bigaus <- function(r,sigma=1) {
  1/(sigma^2)*exp(-(r/sigma)^2/2)
}



## Hankel (2D) transform of exp(-r/scale)/scale
exp.hankel <- function(rho,rate=1,ang.freq=FALSE) {
  if (!ang.freq)
    rho <- 2*pi*rho
  (1+(rho/rate)^2)^(-3/2)
}

gaus.hankel <- function(rho,sigma=1,ang.freq=FALSE) {
  if (!ang.freq)
    rho <- 2*pi*rho
  exp(-((rho*sigma)^2)/2)
}

mspec.hankel.ratio <- function(popspec,envspec,smooth.ratio=TRUE,
                               smooth.df=NULL,spec=FALSE) {
  if (!all(popspec$freq==envspec$freq))
    warning("unequal freqs in mspec.hankel.ratio")
  ratio.0 <- envspec$spec/popspec$spec
  ratio.boot <- NULL
  if (is.null(popspec$allboot))
    warning("allboot not available, can't compute boot ratio")
  if (!is.null(popspec$allboot)) ratio.allboot <- envspec$allboot/popspec$allboot
  if (smooth.ratio) {
    ratio.0 <- smooth.spec(popspec$freq,ratio.0)
    if (!is.null(popspec$allboot))
      ratio.allboot <- t(apply(ratio.allboot,1,
                            function(x)smooth.spec(popspec$freq,
                                                   x,smooth.df)))
    ratio.boot <- env(t(ratio.allboot),ctrfun="none",q=c(0.025,0.975))
  } else {
    ratio.allboot <- NULL
    ratio.boot <- NULL
  }
  list(freq=popspec$freq,ratio=ratio.0,boot=ratio.boot,allboot=ratio.allboot)
}

mspec.hankel <- function(x,y,z,
                         boot=TRUE,nboot=20,
                         seed=1001,freq="auto",
                         fmax=NULL,length.out=NULL,ang.freq=FALSE,
                         norm="none",weight=TRUE,
                         spec=FALSE,
                         smooth=FALSE,
                         all.spec=FALSE,
                         all.boot=FALSE,
                         old.pboot=FALSE,
                         boot.scale=FALSE,
                         boot.tol=1e-6,
                         debug=FALSE) {
  len <- length(x)
  nt <- nrow(z)
  ## standardize z by spatial mean and std. dev
  z <- scale(z,scale=FALSE)
  ## do spectra by direct Hankel transform
  if (freq=="auto") {
    f <- calc.freq(x,length.out=length.out,fmax=fmax,ang.freq=ang.freq)
  } else{
    f <- freq
  }
  f <- f[f>0] ## eliminate zero ordinate
  flen <- length(f)
  smooth.df <- round(sqrt(flen))
  allspec <- apply(z,1,function(Z)spec.hankel2(x,y,Z,
                                               demean=FALSE,
                                               smooth=smooth,smooth.df=smooth.df,
                                               weight=weight,norm=norm,
                                               freq=f,spec=spec,
                                               debug=debug)$spec)
  spec0 <- apply(allspec,1,mean,na.rm=TRUE)
  mlims <- spec0+t(matrix(rep(apply(allspec,1,sd,na.rm=TRUE),2),nrow=2,byrow=TRUE)*qt(c(0.025,0.975),nt-1))
  if (boot) {
    set.seed(seed)
    boot.res <- matrix(nrow=nboot,ncol=flen)
    boot.nok <- numeric(nboot)
    if (old.pboot)
      boot.z <- paramcovboot.old(z,nboot)
    else
      boot.z <- paramcovboot(z,nboot,tol=boot.tol)
    ## take direct Hankel periodograms of all bootstraps
    for (b in 1:nboot) {
      if (debug) cat(b,"\n")
      allspec <- apply(boot.z[,,b],1,
                       function(Z){ if (boot.scale) Z <- scale(Z,scale=FALSE);
                                    spec.hankel2(x,y,Z,
                                                 demean=FALSE,
                                                 smooth=smooth,smooth.df=smooth.df,
                                                 weight=weight,norm=norm,
                                                 freq=f,spec=spec)$spec})
      boot.res[b,] <- apply(allspec,1,mean,na.rm=TRUE)
      boot.nok[b] <- sum(!is.na(allspec[1,]))
      boot.env <-  env(t(boot.res),ctrfun="none",q=c(0.025,0.975))
    }
  } else { ## if not boot
    boot.res <- NA
    boot.nok <- NA
    boot.env <- NULL
  }
  if (!debug) {
    ret <- list(freq=f,spec=spec0,
                speclim=mlims,boot=boot.env,boot.nok=boot.nok)
    if (all.spec) ret$allspec <- t(allspec)
    if (all.boot) ret$allboot <- boot.res
  }
  else {
    ret <- list(freq=f,spec=spec0,bootvals=boot.z)
  }
  class(ret) <- "spec.hankel"
  ret
}

plot.spec.hankel <- function(x,bootq=c(0.025,0.975),
                             lty=c(1,2,2),col=1,
                             plot.boot=TRUE,
                             plot.mlim=TRUE,
                             add=FALSE,...) {
  if (plot.boot && is.null(x$boot))
    warning("bootstrap information missing, not plotted")
  plot.boot <- (plot.boot && !is.null(x$boot))
  if (plot.boot) Z <- cbind(x$spec,x$boot)
  else if (plot.mlim) Z <- cbind(x$spec,x$speclim)
  if (!add) {
    if (plot.boot || plot.mlim) {
      matplot(x$freq,Z,lty=lty,col=col,
              xlab="Frequency",ylab="Power",
              type="l",...)
    } else
    plot(x$freq,x$spec,
         lty=lty[1],col=col[1],
         xlab="Frequency",ylab="Power",
         type="l",...)
  } else {  ## (add==TRUE)
    if (plot.boot || plot.mlim) {
      matlines(x$freq,cbind(x$spec,x$boot),
               lty=lty,col=col,...)
    } else
    lines(x$freq,x$spec,
          lty=lty[1],col=col[1],...)
  }
}

bimspec.hankel <- function(pop.x,pop.y,pop.z,
                           env.x,env.y,env.z,
                           fmax=NULL,length.out=NULL,
                           ...) {
  if (is.null(fmax) && is.null(length.out))
    

  smooth.df <- round(sqrt(length(spec.0$freq)))
  allspec <- apply(z,1,function(Z)spec.hankel2(x,y,Z,
                                               demean=FALSE,
                                               smooth=smooth,smooth.df=smooth.df,
                                               weight=weight,norm=norm,
                                               fmax=fmax,length.out=length.out,
                                               spec=spec)$spec)
  if (collapse) {
    spec.0$spec <- apply(allspec,1,mean,na.rm=TRUE)
  } else {
    spec.0$spec <- t(allspec)
  }
  flen <- length(spec.0$freq)
  if (boot) {
    set.seed(seed)
    boot.res <- matrix(nrow=nboot,ncol=flen)
    boot.nok <- numeric(nboot)
    if (old.pboot)
      boot.z <- paramcovboot.old(z,nboot)
    else
      boot.z <- paramcovboot(z,nboot,tol=boot.tol)
    ## take direct Hankel periodograms of all bootstraps
    for (b in 1:nboot) {
      if (debug) cat(b,"\n")
      allspec <- apply(boot.z[,,b],1,
                       function(Z){ if (boot.scale) Z <- scale(Z,scale=FALSE);
                                    spec.hankel2(x,y,Z,
                                                 demean=FALSE,
                                                 smooth=smooth,smooth.df=smooth.df,
                                                 weight=weight,norm=norm,
                                                 fmax=fmax,length.out=length.out,
                                                 spec=spec)$spec})
      boot.res[b,] <- apply(allspec,1,mean,na.rm=TRUE)
      boot.nok[b] <- sum(!is.na(allspec[1,]))
    }
  } else { ## if not boot
    boot.res <- NA
    boot.nok <- NA
  }
  fvec <- spec.0$freq
  if (!debug) {
    ret <- list(freq=fvec,spec=spec.0$spec,boot=boot.res,boot.nok=boot.nok)
  }
  else {
    ret <- list(freq=fvec,spec=spec.0$spec,bootvals=boot.z)
  }
  class(ret) <- "spec.hankel"
  ret
}

## like mSncf but with no smoothing
empcov2 <- function(x,y,z,mean=FALSE,bias=TRUE,
                    scale=FALSE,cor=FALSE) {
  npts <- length(x)
  if (bias) {
    bmean <- mean.default
  } else {
    bmean <- function(x) {sum(x)/(length(x)-1)}
  }
  if (scale) z <- scale(z,scale=FALSE)
  rmat <- matrix(nrow=npts,ncol=npts)
  for (i in 1:(npts-1)) {
    for (j in (i+1):npts) {
      rmat[j,i] <- sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)
    }
  }
  distvec <- rmat[lower.tri(rmat)]
  if (!cor)
    xcov <- cov(z,use="pairwise.complete.obs")
  else
    xcov <- cor(z,use="pairwise.complete.obs")
  xcov2 <- xcov[lower.tri(xcov)]
  if (!mean) {
    ret <- list(r=sort(distvec),acf=xcov2[order(distvec)])
  } else {
    acov <- sapply(split(xcov2,factor(distvec)),bmean)
    ret <- list(r=sort(unique(distvec)),acf=acov)
  }
  ret
}
