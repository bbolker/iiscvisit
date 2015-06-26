##library(landsc)  ## for evenize(), acf2d()
## library(bbmisc)  ## for chop(), clean.args()
## generic chop() function for making nearly-complex values real
## consider chopping small Real values too?  (e.g., 1e-15+0.3*I -> 0.3*I) -- more cosmetic
chop <- function(x,fuzz=1e-10) {
  if (is.numeric(x)) {  ## pass through unchanged
    x
  }
  else if (is.complex(x)) {
    ## Arg should work, but runs into troubles if ~ pi
    ifelse(abs(Im(x))<abs(Re(x))*fuzz,Re(x),x)
  }
}  

## repeated from bbfuns/misc.R
clean.args <- function(argstr,fn,exclude.repeats=FALSE,
                       exclude.other=NULL,dots.ok=TRUE) {
  fnargs <- names(formals(fn))
  if (length(argstr)>0 && !("..." %in% fnargs && dots.ok))  {
    badargs <- names(argstr)[!sapply(names(argstr),"%in%",c(fnargs,""))]
    for (i in badargs)
      argstr[[i]] <- NULL
  }
  if (exclude.repeats) {
    ntab <- table(names(argstr))
    badargs <- names(ntab)[ntab>1 & names(ntab)!=""]
    for (i in badargs)
      argstr[[i]] <- NULL
  }
  for (i in exclude.other)  ## additional arguments to exclude.other
    argstr[[i]] <- NULL
  argstr
}

library(RandomFields) ## for simulating
library(splines)  ## for smoothed spectra

source("toben.r")
source("rickercml.r")
#source("/home/ben/work/moran/moran-funs.R")
source("rndfield.R")
#source("/home/ben/work/moran/lomb.R")

spfun6.ratio <- function(popspec,envspec,smooth.ratio=TRUE,
                         smooth.df=NULL,spec=FALSE) {
  if (!all(popspec$freq==envspec$freq))
    warning("unequal freqs in spfun6.ratio")
  ratio.0 <- popspec$spec/envspec$spec
  ratio.boot <- NULL
  if (!is.na(popspec$boot)) ratio.boot <- popspec$boot/envspec$boot
  if (smooth.ratio) {
    ratio.0 <- smooth.spec(popspec$freq,ratio.0)
    if (!is.na(popspec$boot))
      ratio.boot <- t(apply(ratio.boot,1,
                            function(x)smooth.spec(popspec$freq,
                                                   x,smooth.df)))
  }
  list(freq=popspec$freq,ratio=ratio.0,boot=ratio.boot)
}

sum.run <- function(sim,dim=1,xmax=40,resamp=20,
                   noisy=FALSE,
                   calc.cov=TRUE,smooth.spec=FALSE,
                   smooth.ratio=TRUE,norm="integrate",
                   specfun=mspec.hankel,boot=TRUE,
                   ratiofun=mspec.hankel.ratio,
                   sqrtspec=FALSE) {
    ## calc cov of pop and env
  if (calc.cov) {
    if (noisy) cat("pop cov\n")
  tstr.cov <- mSncf(sim$x,sim$y,t(sim$pop),xmax=xmax,resamp=resamp,
                    save=TRUE)
  if (noisy) cat("env cov\n")
  tstr.ecov <- mSncf(sim$x,sim$y,t(sim$env),xmax=xmax,resamp=resamp,
                     save=TRUE)
  } else {
    tstr.cov <- NULL
    tstr.ecov <- NULL
  }
  if (noisy) cat("spectra\n")
  ## calculate spectra
    if (noisy) cat("... pop\n")
  tstr.spec <- specfun(x=sim$x,y=sim$y,z=sim$pop,
                       nboot=resamp,
                       spec=TRUE,norm=norm,
                       boot=boot,allboot=TRUE,
                       smooth=smooth.spec)
  if (noisy) cat("... env\n")
  tstr.espec <- specfun(x=sim$x,y=sim$y,z=sim$env,
                        nboot=resamp,
                        spec=TRUE,norm=norm,
                        boot=boot,allboot=TRUE,
                        smooth=smooth.spec)
  ## make sure to use raw spectra for ratio
  if (smooth.ratio && smooth.spec) {
      if (noisy) cat("raw spectra\n")
        if (noisy) cat("... pop\n")
      tstr.rspec <- specfun(x=sim$x,y=sim$y,z=sim$pop,
                            nboot=resamp,
                            spec=TRUE,norm=norm,boot=boot,
                            allboot=TRUE,
                            smooth=FALSE)
  if (noisy) cat("... env\n")
      tstr.respec <- specfun(x=sim$x,y=sim$y,z=sim$env,
                             nboot=resamp,
                             spec=TRUE,norm=norm,boot=boot,
                             allboot=TRUE,
                             smooth=FALSE)
    } else {
      tstr.rspec <- tstr.spec
      tstr.respec <- tstr.espec
    }
  if (sqrtspec) {
    tstr.spec <- sqrt(tstr.spec)
    tstr.espec <- sqrt(tstr.espec)
    tstr.rspec <- sqrt(tstr.rspec)
    tstr.respec <- sqrt(tstr.respec)
  }
  tstr.ratio <- ratiofun(tstr.rspec,tstr.respec,smooth.ratio=smooth.ratio)
  ret <- list(sim=sim,cov=list(pop=tstr.cov,env=tstr.ecov),
       spec=list(pop=tstr.spec,env=tstr.espec),
       specratio=tstr.ratio)
  ret
}

do.run <- function(...,dim=1,xmax=40,resamp=20,
                   plot.it=TRUE,
                   plot.theor=TRUE,noisy=FALSE,
                   calc.cov=TRUE,smooth.spec=FALSE,
                   smooth.ratio=TRUE,
                   specfun=mspec.hankel,
                   norm="integrate",boot=TRUE,weight=TRUE,
                   ratiofun=mspec.hankel.ratio,
                   sqrtspec=FALSE) {
  require(bbmisc) ## for clean.args()
  if (noisy) cat("run sim\n")
  tstr <- runsim(...,dim=dim,save.env=TRUE)
  ## calc cov of pop and env
  if (noisy) cat("pop cov\n")
  tstr.cov <- mSncf(tstr$x,tstr$y,t(tstr$pop),xmax=xmax,resamp=resamp,
                    save=TRUE)
  if (noisy) cat("env cov\n")
  tstr.ecov <- mSncf(tstr$x,tstr$y,t(tstr$env),xmax=xmax,resamp=resamp,
                     save=TRUE)
  if (noisy) cat("spectra\n")
  ## calculate spectra
  if (noisy) cat("... pop\n")
  tstr.spec <- specfun(x=tstr$x,y=tstr$y,z=tstr$pop,
                       nboot=resamp,all.boot=TRUE,
                       spec=TRUE,norm=norm,boot=boot,weight=weight,
                       smooth=smooth.spec)
  if (noisy) cat("... env\n")
  tstr.espec <- specfun(x=tstr$x,y=tstr$y,z=tstr$env,
                        nboot=resamp,all.boot=TRUE,
                        spec=TRUE,norm=norm,boot=boot,weight=weight,
                        smooth=smooth.spec)
  ## make sure to use raw spectra for ratio
  if (smooth.ratio && smooth.spec) {
      if (noisy) cat("raw spectra\n")
      if (noisy) cat("... pop\n")
      tstr.rspec <- specfun(x=tstr$x,y=tstr$y,z=tstr$pop,
                            nboot=resamp,all.boot=TRUE,
                            boot=boot,weight=weight,norm=norm,
                            spec=TRUE,
                            smooth=FALSE)
      if (noisy) cat("... env\n")
      tstr.respec <- specfun(x=tstr$x,y=tstr$y,z=tstr$env,
                             nboot=resamp,all.boot=TRUE,
                             boot=boot,weight=weight,norm=norm,
                             spec=TRUE,
                             smooth=FALSE)
    } else {
      tstr.rspec <- tstr.spec
      tstr.respec <- tstr.espec
    }
  transspec <- function(x) {
    list(freq=x$freq,spec=sqrt(x$spec),boot=sqrt(x$boot),boot.nok=boot.nok)
  }
  if (sqrtspec) {  ## transform/take square roots of spectra
    tstr.spec <- transspec(tstr.spec)
    tstr.espec <- transspec(tstr.espec)
    tstr.rspec <- transspec(tstr.rspec)
    tstr.respec <- transspec(tstr.respec)
  }
  tstr.ratio <- ratiofun(tstr.rspec,tstr.respec,smooth.ratio=smooth.ratio)
  ## calculate theoretical values
  arglist <- list(...)
  ## delete sim-only values: KLUGE (now done by clean.args())
#   arglist$"nloc" <- NULL
#   arglist$"nt" <- NULL
#   arglist$"xlen" <- NULL
#   arglist$"ylen" <- NULL
#   arglist$"wrap" <- NULL
#   arglist$"equispace" <- NULL
#   arglist$"save.env" <- NULL
#   arglist$"init.sd" <- NULL
#   arglist$"init.mean" <- NULL
#  arglist$"trans" <- NULL
  ## need to protect vectors with list()
  theor <- do.call("theor.vals",c(clean.args(arglist,theor.vals),
                                  lag=list(tstr.cov$real$predicted$x),
                                  freq=list(tstr.spec$freq),
                                  norm=norm))
  ##
  ret <- list(sim=tstr,cov=list(pop=tstr.cov,env=tstr.ecov),
       spec=list(pop=tstr.spec,env=tstr.espec),
       specratio=tstr.ratio,theor=theor)
  if (plot.it) {
    plot.run(ret,plot.theor=plot.theor,dim=dim)
  }
  invisible(ret)
}

plot.run <- function(x,fmax=NULL,plot.theor=TRUE,
                     psfile=NULL,legend=FALSE,
                     imgcol=gray((0:100/100)),
                     title=FALSE,mar=c(5,4,2,2)+0.1,
                     dim=1,lwd=1,...) {
  if (is.null(fmax)) fmax <- max(x$spec$pop$freq)
  par(bty="l",cex=1.5)
  ## pop densities
  if (!is.null(psfile))
    postscript(file=paste(psfile,"a.ps",sep=""))
  par(mar=mar)
  if (dim==1) {
      image(t(x$sim$pop),xlab="Space",ylab="Time",main="Pop. densities",
            col=imgcol)
    } else {
      ## how can I plot 2D sim results???
      plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab="",ylab="")
      box()
      text(0.5,0.5,c("No simulation plot"),cex=2)
    }
  if (!is.null(psfile)) dev.off()
    ## covariances
  if (!is.null(x$cov$pop) || !is.null(x$cov$env)) {
    if (!is.null(psfile)) postscript(file=paste(psfile,"b.ps",sep=""))
    par(mar=mar)
    plot(x$cov$pop,text=FALSE,lwd=lwd,...)
    if (title==TRUE) title(main="Covariances")
    plot(x$cov$env,add=TRUE,col="red",text=FALSE,...)
    if (plot.theor) {
      lines(x$theor$env$cov$lag,x$theor$env$cov$acf,col="red",lty=3,lwd=lwd)
      lines(x$theor$pop$cov$lag,x$theor$pop$cov$acf,col="black",lty=3,lwd=lwd)
    }
    if (!is.null(psfile)) dev.off()
  }
  ## spectra
  if (!is.null(psfile)) postscript(file=paste(psfile,"c.ps",sep=""))
  par(mar=mar)
  if (!is.na(x$spec$pop$boot)) {
    matplot(x$spec$pop$freq,sqrt(x$spec$pop$boot),type="l",
            lty=2,col=1,log="y",xlim=c(0,fmax),
            main="",xlab="Frequency",ylab="Transform",lwd=lwd,...)
    lines(x$spec$pop$freq,sqrt(x$spec$pop$spec),col="black",lwd=lwd)
  }
  else plot(x$spec$pop$freq,x$spec$pop$spec,col="black",type="l",log="y",
            xlab="Frequency",ylab="Power",lwd=lwd,...)
  if (!is.na(x$spec$env$boot)) {
    matlines(x$spec$env$freq,sqrt(x$spec$env$boot),
             col="red",lty=2,lwd=lwd)
    lines(x$spec$env$freq,sqrt(x$spec$env$spec),col="red",lwd=lwd)
  }
  else lines(x$spec$env$freq,x$spec$env$spec,col="red")
  if (plot.theor) {
    lines(x$theor$pop$trans$freq,x$theor$pop$trans$trans*sqrt(x$spec$pop$spec[1]),lty=3,col="black")
    lines(x$theor$env$trans$freq,x$theor$env$trans$trans*sqrt(x$spec$env$spec[1]),lty=3,col="red")
  }
  if (!is.null(psfile)) dev.off()
  ## ratios  
  if (!is.null(psfile))
    postscript(file=paste(psfile,"d.ps",sep=""))
  par(mar=mar)
  plot.ratio(x,plot.theor=plot.theor,fmax=fmax,lwd=lwd,neg=TRUE,...)
  if (!is.null(psfile)) dev.off()
}

plot.ratio <- function(x,plot.theor=FALSE,
                       wavelen=FALSE,fmax=NULL,title=FALSE,
                       neg=FALSE,xlab="Frequency",ylab="Ratio",...) {
  if (is.null(fmax)) fmax <- max(x$specratio$freq)
  if (wavelen) {
    xval <- 1/x$specratio$freq
    xlim <- NULL
  } else {
    xval <- x$specratio$freq
    xlim <- c(0,fmax)
  }
  if (title) main <-  "Spectral ratio & bootstrap envelope" else main <- ""
  if (!is.null(x$specratio$boot)) {
    ratio.env <- sqrt(x$specratio$boot)
    ratio.val <- sqrt(x$specratio$ratio)
    if (neg) {
      ratio.env <- -ratio.env
      ratio.val <- -ratio.val
    }
    matplot(xval,ratio.env,
            type="l",lty=2,col=1,ylim=range(ratio.val),
            xlim=xlim,
            xlab=xlab,ylab=ylab,
            main=main,
            axes=FALSE,
            ...)
    axis(side=1)
    if (neg) axis(side=2,at=pretty(ratio.val),labels=as.character(-pretty(ratio.val)))
    else axis(side=2)
    polygon(c(xval,rev(xval)),
            c(ratio.env[,1],rev(ratio.env[,2])),col="gray",border=0)
    matlines(xval,ratio.env,
             lty=2,col=1)
    lines(xval,ratio.val,col="red")
  } else plot(xval,ratio.val,col="red",type="l",
              xlab="Frequency",ylab="Ratio")
  abline(h=mean(ratio.val),col="blue",lty=2)
  if (plot.theor && !is.null(x$theor)) {
    thval <- sqrt(x$theor$ratio$spec.ratio)
    if (neg) thval <- -thval
    if (wavelen) {
      lines(1/x$theor$ratio$freq,thval,lty=1,col="green")
    } else {
      lines(x$theor$ratio$freq,thval,lty=1,col="green")
    }
  }
}

plot.theor <- function(x,wavelen=FALSE) {
  matplot(x$theor$env$cov$lag,
          cbind(x$theor$env$cov$acf,x$theor$pop$cov$acf),
          col=c("red","black"),lty=1,type="l",
          xlab="lag",ylab="ACF")
  matplot(x$theor$env$spec$freq,
          cbind(x$theor$env$spec$spec,x$theor$pop$spec$spec),
          col=c("red","black"),lty=1,type="l",
          xlab="Frequency",ylab="Power",log="y")
  if (wavelen) {
      plot(1/x$theor$ratio$freq,x$theor$ratio$spec.ratio,lty=1,col="green",
           xlab="Frequency",ylab="Ratio")
    } else {
      plot(x$theor$ratio$freq,x$theor$ratio$spec.ratio,lty=1,col="green",
           xlab="Frequency",ylab="Ratio",type="l")
    }
}

stoch.ricker2 <- function(x, r, noise, K=1, add.noise=FALSE) {
  if (add.noise) {
    r <-   x * exp(r * (K - x)) + noise
  } else {
    r <- x * exp(r * (K - x) + noise)
  }
  r
}

stoch.clogist <- function(x, r, noise, K=1, add.noise=FALSE) {
  if (add.noise) {
    dx <-   r*x * (1-x/K) + noise
  } else {
    dx <- r*x*(1 - x/(K*exp(noise)))
  }
  dx
}

nexp <- function(r,scale=1) {
  1/scale*exp(-r/scale)
}

runsim <- function(nloc=100,nt=30,
                   trans=0,
                   xlen=nloc,
                   ylen=nloc,
                   r=2,
                   K=1,
                   init.sd=0.1,
                   init.mean=K,
                   move=0.4,
                   disp.fun=nexp,
                   disp.scale=1,
                   wrap=TRUE,
                   env.scale=1,
                   env.cov="exp",
                   env.sd=1,
                   env.noise=1,
                   equispace=TRUE,
                   dim=1,
                   seed=1001,
                   save.env=FALSE,
                   disc=TRUE,
                   dt=0.05) {
  set.seed(seed)
  if (equispace) xloc <- (1:nloc)*xlen/nloc
  else xloc <- runif(nloc)*xlen
  if (dim==1)
    yloc <- rep(0,nloc)
  else if (equispace)
    yloc <- (1:nloc)*ylen/nloc
  else
    yloc <- runif(nloc)*ylen
  ## distance matrix
  dist <- sqrt(outer(xloc,xloc,"-")^2+outer(yloc,yloc,"-")^2)
  movemat <- apply(dist,c(1,2),disp.fun,disp.scale)
  ## normalize: relative probs of moving to each location, given move
  diag(movemat) <- 0
  movemat <- sweep(movemat,1,apply(movemat,1,sum),"/")*move
  if (!disc) movemat <- movemat*dt
  if (disc) {
    diag(movemat) <- 1-move
  } else {
    diag(movemat) <- 1-(move*dt)
  }
  X <- matrix(ncol=nloc,nrow=nt)
  envmat <- NULL
  if (save.env)
    envmat <- matrix(ncol=nloc,nrow=nt-1)
  ## SVD matrices for environmental conditions
  EV <- rmvn.spa(xloc,yloc,p=env.scale,ret.vsqrt=TRUE)
  ## initial conditions: random
  X[1,] <- rnorm(nloc,mean=init.mean,sd=init.sd)
  if (trans>0) {
    Xtrans <- X[1,]
    for (time in 0:trans) {
      env <- rmvn.spa(xloc,yloc,p=env.scale,vsqrt=EV)
      if (disc) { 
        Xtrans <- stoch.ricker2(Xtrans,r,env.sd*env,K)
        ## movement
        Xtrans <- movemat %*% Xtrans
      } else {
        for (j in 1:round(1/dt)) {
          dx <- stoch.clogist(Xtrans,r,env.sd*env,K)*dt
          + movemat %*% Xtrans
          Xtrans <- pmax(0,Xtrans+dx)
        }
      }
    }
    X[1,] <- Xtrans
  }
  if (nt>=2)
    for (time in 2:nt) {
      env <- rmvn.spa(xloc,yloc,p=env.scale,vsqrt=EV)
      if (save.env)
        envmat[time-1,] <- env
      if (disc) {
        ## dynamics+noise
        X[time,] <- stoch.ricker2(X[time-1,],r,env.sd*env,K)
        ## movement
        X[time,] <- movemat %*% X[time,]
      } else {
        X[time,] <- X[time-1,]
        for (j in 1:round(1/dt)) {
          env <- rmvn.spa(xloc,yloc,p=env.scale,vsqrt=EV)
          dx <- stoch.clogist(X[time,],r,env.sd*env,K)*dt
          + movemat %*% X[time,]
          X[time,] <- pmax(0,X[time,]+dx)
        }
      }
    }
  list(x=xloc,y=yloc,pop=X,env=envmat)
}

## calculate theoretical covariances, transforms, and spectra
theor.vals <- function(lag=NULL,dlag=NULL,maxlag=NULL,
                       freq=NULL,dfreq=NULL,maxfreq=NULL,
                       r=2,move=0.4,
                       disp.fun=nexp,disp.scale=1,
                       env.scale=1,env.cov="exp",
                       env.sd=1,env.noise=1,
                       dim=1,
                       norm="integrate",
                       ang.freq=FALSE) {
  ## figure out max dist, lag
  if (is.null(lag))  lag <- seq(0,maxlag,by=dlag)
  if (is.null(freq)) freq <- seq(0,maxfreq,by=dfreq)
  if (dim==1) {
    specfun <- exptr.1d
  } else {
    specfun <- exp.hankel
  }
  if (env.cov=="exp") {
    ## set up environmental covariance and spectrum
    cov.env <- exp(-lag/env.scale)
    spec.env <- env.sd^2*specfun(freq,1/env.scale,ang.freq=ang.freq)
    trans.env <- sqrt(spec.env)
  }
  ## assume exp. dispersal for now
  trans.disp <- specfun(freq,1/disp.scale,ang.freq=ang.freq)
  spec.disp <- trans.disp^2
  trans.pop <- trans.env/(2*((r+move)-move*trans.disp))
  spec.pop <- trans.pop^2
  ## numerical inversion of pop spectrum
  if (dim==1) {
    cov.pop <-  (Re(fft(spec.pop,inverse=TRUE))*dfreq)[1:length(lag)]
  } else {
    ## invert
    cov.pop <-  spec.hankel(spec.pop,freq,lag)$spec ## normalization etc.??
    cov.pop <- cov.pop/cov.pop[1]
  }
  trans.ratio <- trans.env/trans.pop
  spec.ratio <- spec.env/spec.pop
  list(env=list(cov=list(lag=lag,acf=cov.env),
         spec=list(freq=freq,spec=spec.env),
         trans=list(freq=freq,trans=trans.env)),
       pop=list(cov=list(lag=lag,acf=cov.pop),
         spec=list(freq=freq,spec=spec.pop),
         trans=list(freq=freq,trans=trans.pop)),
       disp=list(spec=list(freq=freq,spec=spec.disp),
         trans=list(freq=freq,trans=trans.disp)),
       ratio=list(freq=freq,trans.ratio=trans.ratio,
         spec.ratio=spec.ratio))
}

## set up a 1D landscape, evenly spaced (exp. covariance)
xcov.1d <- function(dx=0.1,sclen=20,scale=1) {
  y <- seq(0,sclen,by=dx)
  x <- rep(0,length(y))
  z <- rmvn.spa(x,y,p=scale)
  list(x=x,y=y,z=z)
}

## 1D landscape, unevenly spaced (uniform distribution along transect),
## exponential covariance
xcov2.1d <- function(dx=0.1,sclen=20,scale=1) {
  n <- round(sclen/dx)+1
  y <- runif(n)*sclen
  x <- rep(0,n)
  z <- rmvn.spa(x,y,p=scale)
  list(x=x,y=y,z=z)
}

## 2D landscape, unevenly spaced (uniform distribution along transect),
## exponential covariance
xcov2 <- function(dx=0.1,sclen=20,scale=1,dim=1) {
  n <- round(sclen/dx^dim)+1
  y <- runif(n)*sclen
  if (dim==1) {
      x <- rep(0,n)
  } else {
    x <- runif(n)*sclen
  }
  z <- rmvn.spa(x,y,p=scale)
  list(x=x,y=y,z=z)
}


## 1D transform of exponential covariance (Cauchy)
exptr.1d <- function(x,scale=1,ang.freq=FALSE,
                     dx=0.1,norm=FALSE) {
  f <- 1
  if (!ang.freq) f <- 2*pi
  r <- scale/(scale^2+(f*x)^2)+dx/2
  if (norm) r <- r*4/f
  r
}

matrix.to.array <- function(m,nblocks,debug=FALSE) {
  ## break a matrix m into an array in n blocks determined
  ## by rows
  nr <- nrow(m)
  if (debug) cat("matrix.to.array:",dim(m),ncol(m),nr/nblocks,nblocks,"\n")
  aperm(array(t(m),c(ncol(m),nr/nblocks,nblocks)),
        c(2,1,3))
}

## spectral function number 6: applies to multivariate
## spatial data (x,y are spatial positions, z is a matrix
## with (length(x)) columns and nt rows
## stripped out a lot of the extra crap:
## only does two dimensions
## only does one kind of bootstrap (parametric/Wishart)
## OBSOLETE
spfun6 <- function(x,y,z,
                   boot=TRUE,nboot=20,
                   seed=1001,
                   debug=FALSE,
                   norm=FALSE,
                   grid=FALSE, dim=2,
                   spec=FALSE,
                   smooth=FALSE,
                   weight=TRUE,
                   smooth.kluge=0) ## dummy parameters
{
  if (grid==TRUE) stop("can't do spfun6 on grid")
  if (dim==1) stop("can't do spfun6 for dim=1")
  if (norm) {
    normspec <- function(x)x/sum(x)
  }
  else {
    normspec <- function(x)x
  }
  len <- length(x)
  nt <- nrow(z)
  ## standardize z by spatial mean and std. dev
  z <- scale(z)
  ## do spectra by direct Hankel transform
  ## do spectrum of first row just to figure out appropriate frequencies
  spec.0 <- spec.hankel2(x,y,z[1,],evenize=FALSE,demean=FALSE)
  smooth.df <- round(sqrt(length(spec.0$freq)))+smooth.kluge
  allspec <- apply(z,1,function(Z)spec.hankel2(x,y,Z,evenize=FALSE,
                                               demean=FALSE,
                                               smooth=smooth,
                                               weight=weight,
                                               smooth.df=smooth.df,
                                               spec=spec)$spec)
  spec.0$spec <- normspec(apply(allspec,1,mean,na.rm=TRUE))
  flen <- length(spec.0$spec)
  if (boot) {
    set.seed(seed)
    boot.res <- matrix(nrow=nboot,ncol=flen)
    boot.nok <- numeric(nboot)
    boot.z <- paramcovboot(z,nboot,randnorm=TRUE)
    ## paramcovboot now returns results as an array -- no conversion
    ## necessary?
    ## matrix.to.array(paramcovboot(z,nboot*nt),nboot)
    ## take direct Hankel periodograms of all bootstraps
    for (b in 1:nboot) {
      if (debug) cat(b,"\n")
#       allspec <- try(apply(boot.z[,,b],1,
#                        function(Z)spec.hankel2(x,y,Z,evenize=FALSE,
#                                                demean=FALSE,
#                                                smooth=smooth,
#                                                smooth.df=smooth.df,
#                                                spec=spec)$spec))
#       ca <- class(allspec)
#       if (!is.null(ca) && ca=="try-error") {
#         boot.res[b,] <- rep(NA,flen)
#       } else {
#         boot.res[b,] <- normspec(apply(allspec,1,mean))
#      }
      allspec <- apply(boot.z[,,b],1,
                       function(Z)spec.hankel2(x,y,Z,
                                               demean=FALSE,
                                               smooth=smooth,
                                               weight=weight,
                                               smooth.df=smooth.df,
                                               spec=spec)$spec)
      boot.res[b,] <- normspec(apply(allspec,1,mean,na.rm=TRUE))
      boot.nok[b] <- sum(!is.na(allspec[1,]))
    }
  } else { ## if not boot
    boot.res <- NA
    boot.nok <- NA
  }
  fvec <- spec.0$freq
  if (!debug)
    list(freq=fvec,spec=spec.0$spec,boot=boot.res,boot.nok=boot.nok)
  else
    list(freq=fvec,spec=spec.0$spec,bootvals=boot.z)
}
