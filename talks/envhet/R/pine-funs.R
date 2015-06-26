## Would be better to get nlme to compute
## the correlation function values for us,
## but this is deeply nested and obscure.
##  actually, it's not as bad as I thought:
##   "Variogram" will construct/plot the
##    (semi?)variogram, and this could be
##    modified/rewritten to produce a correlogram ...
##

require(plyr)  ## for rename
require(reshape) ## for melt
require(ggplot2)
suppressPackageStartupMessages(library(RandomFields))

glsmultfit <- function(...,cormodel,
                       distrange,nugget,form,
                       type="random",ntry=10,
                       value,which=c("range","nugget"),
                       verbose=FALSE) {
  if (missing(value)) {
    value <- c(range=0,nugget=0)
  }
  ranges <- nuggets <- NULL
  if (random) {
    if ("range" %in% which) ranges <- exp(runif(ntry,log(distrange[1]),log(distrange[2])))
    if ("nugget" %in% which) nuggets <- runif(ntry)
    values <- cbind(ranges,nuggets)
  } else {
    ## problem here -- do we need to expand.grid?  ugh ...
    if ("range" %in% which) ranges <- exp(seq(log(distrange[1]),log(distrange[2]),length=ntry))
    if ("nugget" %in% which) nuggets <- seq(0.001,0.999,length=ntry)
  }
  ## UNFINISHED
}

calc_cor <- function(model,params,x,nugget=FALSE) {
  ## scale x by r? would simplify formulas ...
  ## cat(model,"\n")
  if (length(unique(model))>1) stop("something wrong with model")
  corfun <- switch(model[1],
                   corExp=function(x,r) { exp(-x/r) },
                   corRatio=function(x,r) { 1/(1+(x/r)^2) },
                   corSpher = function(x,r) { ifelse(x>r,0,1-1.5*(x/r)+0.5*(x/r)^3) },
                   corLin = function(x,r) { ifelse(x>r,0,1-(x/r))},
                   corGaus = function(x,r) { exp(-(x/r)^2) },
                   stop("unknown model type",model[1])
                   )
  corval <- corfun(x,params[1])
  if (nugget) corval <- corval*(1-params[2])
  corval
}

sim_cor <- function(cs,apVar,dist,nsim=1000,envelope=FALSE) {
  require(MASS)
  nc <- length(coef(cs))
  if (inherits(apVar,"character")) stop("non-pos def var-cov")
  vals <- mvrnorm(nsim,mu=coef(cs),Sigma=apVar[1:nc,1:nc])
  ## back-transform from 'unconstrained' scale
  vals[,1] <- exp(vals[,1])
  vals[,2] <- plogis(vals[,2])
  r <- apply(vals,1,calc_cor,model=class(cs)[1],x=dist,nugget=attr(cs,"nugget"))
  if (envelope) {
    t(apply(r,1,quantile,c(0.025,0.975)))
  } else r
}

plot_pred_cor <- function(X,add=FALSE, 
                          obsticks=FALSE, ## mark locations of observed data?
                          use.nlme=FALSE, ## use nlme functions to extract cov?
                          nsim=0,         ## number of sims
                          envelope=TRUE,  ## sim envelope only?
                          col=1,
                          col.sim=col,
                          lty=1,
                          lty.sim=2,
                          npts=256,...) { ## # of points between min and max
  cs <- X$modelStruct$corStruct
  ## form <- attr(cs,"formula")
  ## extract data variables from formula (ugly)
  ## x <- get(as.character(form[[2]])[2],env=environment(form))
  ##  y <- get(as.character(form[[2]])[3],env=environment(form))
  ## vars <- cbind(x,y)
  ## d1 <- as.numeric(dist(vars))
  d1 <- unlist(getCovariate(cs)) ## better way to get distances!
  dord <- order(d1)
  d1 <- d10 <- sort(d1)
  if (use.nlme) {
    m <- corMatrix(cs)
    corval <- (m[lower.tri(m)])[dord]
  } else {
    d1 <- seq(0,max(d1),length=npts)
    corval <- calc_cor(model=class(cs)[1],x=d1,
                       params=coef(cs,unconstrained=FALSE),
                       nugget=attr(cs,"nugget"))
    if (nsim>0) {
      s <- sim_cor(cs,X$apVar,d1,nsim=nsim,envelope=envelope)
    }
  }
  if (!add) {
    plot(d1,corval,type="l",col=col,lty=lty,...)
    if (obsticks) {
       n <- length(d10)
       u <- par("usr")
       segments(d10,rep(u[3],n),d10,rep(u[3]+0.02,n))
     }
  } else {
    lines(d1,corval,col=col,lty=lty,...)
  }
  if (nsim>0 && !use.nlme) {
    matlines(d1,s,lty=lty.sim,col=col.sim)
  }
}


## Function to run a GLS model for a specified
## model (character as above), nugget (TRUE/FALSE),
## starting value:

gls2 <- function(model,data,cor_modstr,
                 cor_value,cor_nugget,cor_formula,
                 weights=NULL,
                 fun=gls,verbose=FALSE,...) {
  ## browser()
  f <- try(fun(model,data,correlation=get(cor_modstr)(
                            value=cor_value,
                            nugget=cor_nugget,
                            form=cor_formula),
               weights=weights,...))
  if (verbose) cat("xx\n")
  if (inherits(f,"try-error")) return(NA)
  f$call[-1] <- lapply(f$call[-1],eval.parent)
  f$call[[1]] <- substitute(fun)
  f
}

scoef <- function(x,unconstrained=FALSE,...) {
  is.model <- function(x) {
    inherits(x,"gls") || inherits(x,"lme")
  }
  if (is.model(x)) {
    cc <- coef(x$modelStruct$corStruct,unconstrained=unconstrained,...)
    if (length(cc)==1) cc <- c(cc,nugget=NA)  ## add nugget
    cc
  } else if (length(x)>1 && is.model(x[[1]])) {
    t(sapply(x,scoef))
  } else stop("unknown type")
}

head.ICtab <- function(x,delta=10,n,weight,...) {
  if (!missing(weight) & !is.null(x$weight)) {
    keep <- x$weight>weight
  }
  x2 <- lapply(x,"[",keep)
  attr(x2,"row.names") <- attr(x,"row.names")[keep]
  class(x2) <- "ICtab"
  x2
}

glsList <- function(data,f,model,...) {
  z <- lapply(split(data,f),
              gls,model=model,
              ...)
  class(z) <- "glsList"
  z
}

logLik.glsList <- function(object,...) {
  sum(sapply(object,"[[","logLik"))
}

AIC.glsList <- function(object,...) {
  sum(sapply(object,AIC))
}

## would have tried beanplot, but gets hung up trying to compute densities for sparse data
mdiag <- function(model,pscale=c(1,6),alpha=0.7,eps=1e-6) {
  require(ggplot2)
  require(gridExtra)
  data <- model$call$data
  plotord <- order(with(data,tapply(count,list(plot),mean)))
  data <- within(data, {
    resid <- residuals(model,type="pearson")
    sresid <- sqrt(abs(resid))
    fitted <- fitted(model)
    fitted[fitted<eps] <- 0
  })
  g1 <- ggplot(data,aes(x=reorder(plot,count),y=count))+stat_sum(aes(size=..n..),alpha=alpha)+
    scale_size_continuous(to=pscale,legend=FALSE)+theme_bw()
  g2 <- ggplot(data,aes(x=reorder(plot,count),y=resid))+stat_sum(aes(size=..n..),alpha=alpha)+
    scale_size_continuous(to=pscale,legend=FALSE)+theme_bw()+
      geom_smooth(aes(group=1))
  ## browser()
  nlsfit <- try(nls(sresid~fitted^b,data=data,
                algorithm="plinear",start=list(b=0.25),trace=FALSE),silent=TRUE)
  g3 <- ggplot(data,aes(x=fitted,y=sresid))+stat_sum(aes(size=..n..),alpha=alpha)+
    scale_size_continuous(to=pscale,legend=FALSE)+theme_bw()+
      geom_smooth(aes(group=1))+theme_bw()+
        stat_function(fun=sqrt,colour="red")+xlim(0,max(data$fitted))

  if (inherits(nlsfit,"try-error")) warning("nls fit failed") else
  {
    predframe <- data.frame(fitted=seq(0,max(data$fitted),length=101))
    predframe <- cbind(predframe,sresid=predict(nlsfit,newdata=predframe))
    g3 <- g3 + geom_line(data=predframe,colour="purple")
  }
  invisible(nlsfit)
}


## calculate a likelihood profile on the range parameter
rangeprof <- function(model,logr=seq(-1,6.5,by=0.05)) {
  ## don't need this, AND it messes up when apVar is bad
  ## cc <- intervals(model,which="var-cov")$corStruct["range",]
  cc2 <- model$modelStruct$corStruct
  Lfun <- function(r) {
    attr(cc2,"fixed") <- c(TRUE,FALSE) ## allow nugget to adjust
    cc2[1] <- r
    L <- try(update(model,correlation=cc2))
    if (inherits(L,"try-error")) NA else logLik(L)
  }
  Lvals <- sapply(logr,Lfun)
  ## get relative value (subtract max log-likelihood)
  dvec <- -2*(Lvals-max(Lvals,na.rm=TRUE))
  data.frame(logr,logLik=Lvals,deviance=dvec)
}

rangevarprof <- function(model,logr=seq(0,6.5,by=0.05),varp=seq(0,1,by=0.05)) {
  cc <- intervals(model,which="var-cov")$corStruct["range",]
  cc2 <- model$modelStruct$corStruct
  vv2 <- 
  Lfun <- function(r,v) {
    ## may need to make explicit local copy?
    ## experimental:
    ## CC <- model$call
    cc3 <- cc2
    attr(cc3,"fixed") <- c(TRUE,FALSE) ## allow nugget to adjust
    cc3[1] <- r
    ## CC$correlation <- cc3
    ## CC$weights <- varPower(fixed=v)
    u1 <- try(update(model,correlation=cc3,weights=varPower(fixed=v)),silent=TRUE)
    ## eval(CC)
    if (inherits(u1,"try-error")) NA else logLik(u1)
  }
  if (require(tcltk)) {
    pb <- tkProgressBar("range/varp calc")
  } else pb <- NULL
  ntot <- length(logr)*length(varp)
  Lvals <- numeric(ntot)
  k <- 0
  for (j in (seq_along(varp))) {
    for (i in (seq_along(logr))) {  ## i varies faster in expand.grid
      k <- k+1
      Lvals[k] <- Lfun(logr[i],varp[j])
      setTkProgressBar(pb,k/ntot)
    }
  }
  dvec <- -2*(Lvals-max(Lvals,na.rm=TRUE))
  close(pb)
  data.frame(expand.grid(logr=logr,varpower=varp),logLik=Lvals,deviance=dvec)
}



## run likelihood profiles for each of a set of models
allrprof <- function(modellist) {
  rplist <- lapply(modellist,
         function(w) {
           rp <- try(rangeprof(w))
           if (inherits(rp,"try-error")) NULL else rp })
  ## drop bad models
  rplist[!sapply(rplist,is.null)]
}

## convert list of likelihood profiles to 'long form'
convrprof <- function(rproflist) {
  cbind(model=rep(names(rproflist),
          unlist(sapply(rproflist,nrow))),
        do.call(rbind,rproflist))
}

which_cat <- function(x,categories) {
  s <- !is.na(sapply(categories,match,x=x))
  which(s)
}

## extension of base::rapply
##  1. 'classes' may be non-atomic
##  2. 'f' and 'classes' may be specified as lists

## ?? allow generic tests?

Rapply <- function(object, f, classes= "ANY", deflt = NULL, how=c("replace","unlist"), debug=FALSE) {
  if (debug) cat("*\n")
  if (!is.list(classes)) {
    if (class(object) %in% classes || (is.atomic(object) && classes=="ANY")) {
      if (debug) cat(class(object),"found; returning f(object)","\n")
      return(f(object))
    }
  } else {
    if (!is.list(f) || length(f)!=length(classes)) stop("as many functions as sets of classes must be defined")
    which_class <- which_cat(class(object),classes)
    if (debug) cat(class(object),length(which_class),which_class,"\n")
    if (length(which_class)>1) stop("match to multiple classes")
    if (length(which_class)==1) {
      if (debug) cat(class(object)," found; returning f[[",which_class,"]](object)","\n",sep="")
      return(f[[which_class]](object))
    }
  }
  how <- match.arg(how)
  if (is.atomic(object)) {
    if (debug) cat("atomic object, classes!='ANY': returning default value\n")
    return(deflt)
  }
  z <- lapply(object,Rapply,f=f,classes=classes,deflt=deflt,debug=debug)
  if (is.function(how)) how(z) else {
    if (how=="unlist") unlist(z) else z
  }
}

## repeat a function n times: 'generalized power'
gpow <- function(object,f,n,i=0,...) {
  if (i == n) return(object)
  gpow(f(object,...),f,n,i+1,...)
}

## 'generalized unlist'
gunlist <- function(x,n=4) {
  gpow(x,unlist,n,recursive=FALSE)
}

collapse_proflist <- function(proflist,rframe0) {
  P1 <- proflist[[1]][[1]][[1]][[1]][[1]]
  logrvec <- P1$logr
  proflist2 <- Rapply(proflist,list(function(x) x,
                                    function(x) data.frame(logr=logrvec,
                                                           logLik=rep(NA,length(logrvec)),
                                                           deviance=rep(NA,length(logrvec)))),
                      list("data.frame","try-error"))
  proflist_err <- Rapply(proflist,list(function(x) "OK",
                                       as.character),
                         list("data.frame","try-error"),
                         how="unlist")
  ## OK or not?
  proflist_OK <- Rapply(proflist,list(function(x) TRUE,
                                      function(x) FALSE),
                        list("data.frame","try-error"),how="unlist")
  ## collapse down to flat list of data frames (and errors?)
  proflist2B <- gunlist(proflist2)
  ## null or not?
  proflist_notnull <- !sapply(proflist2B,is.null)
  proflist3 <- do.call("rbind",proflist2B)
  ## ID variables
  rframe2 <- as.data.frame(lapply(rframe0[proflist_notnull,],rep,
                                  each=length(logrvec)))
  proflist4 <- na.omit(cbind(rframe2,proflist3))
  proflist4
}

do_sim <- function(range,
                   nugget,
                   m,
                   v,
                   rfun=rpois,
                   invlink=exp,
                   gn=100,
                   samples=75,
                   L=50, ## length scale
                   spmodel="exponential",
                   retval=c("grid","samples","list"),
                   input=NULL,  ## value to modify (e.g. seed numbers)
                   inputpar="size") {
  require(RandomFields)
  retval <- match.arg(retval)
  coords <- (1:gn)*L/gn
  i1 <- GaussRF(x=coords,y=coords,grid=TRUE,
                model=spmodel,
                param=c(mean=m,variance=v,nugget=nugget,scale=range))
  i2 <- invlink(i1)
  if (is.null(input)) {
    ## generate unconditional random output
    r <- rfun(length(i2),c(i2))
  } else {
    if (length(formals(rfun))<3) stop("Bad rfun")
    arglist <- list(length(i2),c(i2),input)
    names(arglist) <- c("n","",inputpar)
    r <- do.call(rfun,arglist)
  }
  if (retval=="grid") {
    matrix(r,nrow=nrow(i2))
  } else if (retval=="list") {
    list(output=matrix(r,nrow=nrow(i2)),env=i2)
  } else {
    s <- sample(length(r),size=samples)
    m <- matrix(r,nrow=nrow(i2))
    d <- data.frame(count=r[s],
                    env=i2[s],
                    x0=coords[row(m)[s]],y0=coords[col(m)[s]])
    ## hack: shouldn't bother having both x/x0 and y/y0
    d$x <- d$x0
    d$y <- d$y0
    d
  }
}

  

rbind.id <- function(x,labels=seq_along(x)) {
  data.frame(.id=rep(labels[seq_along(x)],
               sapply(x,nrow)),
             do.call(rbind.fill,x))
}


## collapse a list into a data frame, labeling as you go
labfun <- function(x,labels) {
  data.frame(do.call(rbind,x),plot=rep(labels[seq_along(x)],
                                sapply(x,nrow)))
}

xapply <- function(FUN,...,FLATTEN=TRUE,MoreArgs=NULL,pb=FALSE,
                   verbose=TRUE,
                   cores=2) {  ## hard-code 'cores' to 2 here
  if (verbose) pb <- FALSE
  hasMC <- require(doMC)
  if (hasMC) {
    registerDoMC(cores=cores)
    pb <- FALSE
  }
  L <- list(...)
  n <- sapply(L,length)
  inds <- do.call(expand.grid,lapply(n,seq))
  if (verbose) cat("xxaply:",nrow(inds),"total rows\n")
  N <- nrow(inds)
  if (pb) pb1 <- txtProgressBar(max=N)
  ## FIXME:: use plyr::mlply instead?  (builds in parallel backend)
  if (!hasMC) {
    retlist <- list()
    for (i in 1:nrow(inds)) {
      if (pb) setTxtProgressBar(pb1,i)
      if (verbose) {
        do.call(cat,c(list("xapply:",i),inds[i,],list("\n")))
      }
      arglist <- mapply(function(x,j) x[[j]],L,as.list(inds[i,]),SIMPLIFY=FALSE)
      if (FLATTEN) {
        retlist[[i]] <- do.call(FUN,c(arglist,MoreArgs))
      }
    }
  } else {
    retlist <- foreach(i=1:nrow(inds)) %do%
    {
      arglist <- mapply(function(x,j) x[[j]],L,as.list(inds[i,]),SIMPLIFY=FALSE)
      if (verbose) {
        do.call(cat,c(list("xapply:",i),inds[i,],list("\n")))
      }
      do.call(FUN,c(arglist,MoreArgs))
    }
  }
  retlist
}

## example
## xapply(paste,as.list(LETTERS[1:4]),as.list(letters[1:3]))

do_analysis <- function(dat,
                        modList, ## fixed-effect models
                        mod2,    ## correlation models
                        valList, ## starting values
                        vpList,  ## heteroscedasticity model (varPower)
                        rframe0,
                        clip=NA,profAIC=10,pb=TRUE,
                        save_full=FALSE,
                        save_best=FALSE) {
  ## must match up correctly with rframe0, based on expand.grid in sim_anal.R/(!!)
  dmodels_batch <- xapply(gls2,modList,mod2,valList,vpList,
                          MoreArgs=list(data=dat,
                            control=glsControl(),
                            cor_nugget=TRUE,cor_formula=~x0+y0|plot),pb=pb)
  ffun <- function(x,fun,n=1) {
    sapply(x,function(x) if ((is.logical(x) && is.na(x)) || is.character(x)) rep(NA,n) else c(fun(x)))
  }
  LLvec <- ffun(dmodels_batch,logLik)
  AICvec <- ffun(dmodels_batch,AIC)
  rangevec <- ffun(dmodels_batch,function(x) scoef(x)["range"])
  nuggetvec <- ffun(dmodels_batch,function(x) scoef(x)["nugget"])
  d_rframe <- cbind(rframe0,LL=LLvec,dAIC=AICvec-min(AICvec,na.rm=TRUE),
                    range=rangevec,nugget=nuggetvec)
  d_rframe$bothmod <- with(d_rframe,vpmod:fixmod)

  ## fancy wrapper to apply a function and return NA when one
  ##  of several bad things goes wrong ...
  ffun2 <- function(x,fun) {
    lapply(x,function(x) {
      if (is.na(x) || is.logical(x)) return(NA)
      z <- c(fun(x))
      if (!is.numeric(z)) return(NA)
      z
    })
  }

  ## apVar holds information on the approximate variance-covariance matrix of
  ##   the parameters; the attr(,"Pars") component is ??? I think it's the
  ##   estimates of the parameters, but what scale are they on??
  plist <- ffun2(dmodels_batch,function(x) c(x$apVar,attr(x$apVar,"Pars")))
  nmax <- max(sapply(plist,length))
  pmat <- t(sapply(plist,rep,length.out=nmax))

  ## N.B.!
  if (!is.na(clip)) {
    clipvals <- abs(log10(d_rframe$range))>clip
    d_rframe$range[clipvals] <- NA  ## get rid of extreme values
  }
  
  ## profiling ... takes a little while ...


  ## only profile 'adequate' fits and those that don't basically give the same
  ##   range/nugget/AIC values ...
  tmpframe0 <- as.data.frame(lapply(subset(d_rframe,select=c(dAIC,range,nugget)),round,3))
  whichprof <- which(tmpframe0$dAIC<profAIC & !duplicated(tmpframe0))
  
  logr <- seq(0,5,by=0.05)
  ## calculate range profiles for each model
  d_proflist <- ffun2(dmodels_batch[whichprof],
                      function(x) try(rangeprof(x,logr),silent=TRUE))
  ## construct a matrix of all of the deviances (as a function of range)
  ##  of all of the models
  m0 <- do.call(cbind,lapply(d_proflist,
                             function(x) {
                               if (is.logical(x) && is.na(x)) rep(NA,length(logr)) else x[["deviance"]]
                             }))
                 
  ## for (i in seq_along(d_proflist)) {
  ##   x <- d_proflist[[i]]
  ##   m0[[i]] <- if (is.logical(x) && is.na(x)) NULL else {
  ##     data.frame(i,do.call(cbind,x))
  ##   }
  ## }
  ## m0 <- do.call(rbind,m0)

  rownames(m0) <- logr
  d_profmat <- list(x=logr,y=m0)
  d_profframe <- merge(data.frame(n=seq(whichprof),rframe0[whichprof,]),rename(melt(m0),c(X1="logr",X2="n")))
                         
  L <- list(rframe=d_rframe,
       proflist=d_proflist,whichprof=whichprof,profframe=d_profframe,profmat=d_profmat,pmat=pmat)
  ## ?? could COLLAPSE all similar fits that differ only in starting points ...
  if (save_full) {
    L <- c(list(fullbatch=dmodels_batch),L)
  } else if (save_best) {
    best <- which.min(L$rframe$dAIC)
    L <- c(list(bestfit=dmodels_batch[[best]]),L)
  }
  L
}

get_cpars <- function(r,y,nugget=TRUE,plot.it=FALSE) {
  if (!nugget) {
    lm1 <- lm(log(y)~r-1)
    rng <- -1/coef(lm1)[1]
    nugget <- 0
  } else {
    lm1 <- lm(log(y)~r)
    rng <- -1/coef(lm1)[2]
    nugget <- 1-exp(coef(lm1)[1])
  }
  if(plot.it) {
    par(las=1,bty="l")
    plot(r,y,log="y")
    curve((1-nugget)*exp(-x/rng),add=TRUE)
  }
  c(rng,nugget)
}

plot_sdata <- function(d,reorder=TRUE,empty=TRUE) {
  if (reorder) d <- transform(d,plot=reorder(plot,count))
  if (empty) d$empty <- factor(d$count==0,levels=c(TRUE,FALSE))
  g1 <- ggplot(d,aes(x=x0,y=y0,size=count))+
    facet_wrap(~plot)+
      theme_bw()+theme(panel.margin=unit(0,"lines"))+
        xlab("")+ylab("")
  if (empty) {
    g1+ geom_point(aes(colour=empty),alpha=0.5)+
      scale_colour_manual(values=c("red","blue"))
  } else {
    g1+ geom_point(alpha=0.5)
  }
}

prodfun <- function(n,prob,size) {
    size*prob  ## pseudo-binomial but just take *expected* value (product)
}
ifun <- function(n,lambda) {
    lambda   ## pseudo-Poisson but just take expected value
}



fixpmat <- function(x) {
  nmax <- max(sapply(x,length))
  fixfun <- function(x) {
    if (!is.numeric(x) || length(x)<nmax) {
      rep(NA,nmax)
    }
  }
  t(sapply(x,fixfun))
}

### ANALYSIS FUNCTIONS
plot_maps <- function() {
## PLOT seed and seedling (simulated) data
  list(plot_sdata(seeds),plot_sdata(seedlings))
}

round_numeric <- function(L,digits=3) {
  as.data.frame(lapply(L,
         function(x) if(is.numeric(x)) round(x,digits) else x))
}

viewL <- function(r,dAIC_cut=10) {
  r <- subset(r,select=c(spmodel,startval,dAIC,range,nugget,vpmod),
         dAIC<dAIC_cut)
  r[order(r$dAIC),]
}

plot_AICvals <- function(X,AICcut=10) {
  x <- X$rframe
  x <- droplevels(subset(x,dAIC<AICcut))
  g1 <- ggplot(x,aes(x=spmodel,y=dAIC,colour=startval,pch=method))+
        geom_point()+
        facet_grid(bothmod~.)+coord_flip()+xlab("Correlation model")
  g2 <- ggplot(x,aes(x=spmodel,y=range,colour=startval,pch=method))+
        geom_point()+facet_grid(bothmod~.)+
        coord_flip()+xlab("Correlation model")
  list(g1,g2)
}

plot_profs <- function(x,AICcut=10,xlim=NULL) {
  ch <- with(x$rframe,dAIC<AICcut)
  logr <- x$proflist[[1]]$logr
  pmat <-with(x$rframe,
              sapply(x$proflist[],
                     "[[","logLik"))

  ## different starting values give DIFFERENT profiles -- a little alarming.
  minval <- min(-pmat)
  matplot(logr,-pmat,type="l",xlim=xlim,ylim=c(minval,minval+6))  ## two different sets:
  points(logr,-pmat[,1])
  abline(v=log(seedparms$r))
  abline(h=minval+1.92)
}

## L = n^2+n = n*(n+1)
## n^2+n-L=0
## (-1 + sqrt(1+4*L))/2
## L = 20 -> 8/2
##
## L = 12 -> 3


get_pmat <- function(x) {
  n <- round((sqrt(4*length(x)+1)-1)/2)
  z <- matrix(x[1:n^2],nrow=n)
  if (is.na(z[1,2]*z[2,1])) z[1,2] <- z[2,1] <- 0 ## KLUGE
  if (any(is.na(z))) return(NA)
  z[1:2,1:2]
}

get_best <- function(x,dAIC_cut=0.01) {
  OK_pmat <- apply(x$pmat,1,function(x) is.numeric(get_pmat(x)))
  w <- which(x$rframe$dAIC<dAIC_cut & OK_pmat)[1]
  if (is.na(w)) {
    mindAIC <- min(x$rframe$dAIC[OK_pmat])
    stop("can't find a suitable example: increase dAIC_cut? ",
         round(mindAIC,3))
  }
  w
}

fixcormod <- function(m) paste("cor",gsub("N$","",m),sep="")
  

## generate a single set of cor (or cov) ests
gen_cors <- function(ss_seed,ss_seedling,Nvar,Svar,Nbar,Ebar,rvec=rvec) {
  cc_seed_est     <- cf2(rvec,Nvar,unlist(ss_seed["nugget"]),unlist(ss_seed["range"]),
                         fixcormod(unlist(ss_seed["spmodel"])))
  cc_seedling_est <- cf2(rvec,Svar,unlist(ss_seedling["nugget"]),unlist(ss_seedling["range"]),
                         fixcormod(unlist(ss_seedling["spmodel"])))
  cc_env_est <- (cc_seedling_est-Ebar^2*cc_seed_est)/Nbar^2
  corfun <- function(x) x/x[1]
  data.frame(w="est",
             rbind(data.frame(type="seedling",value=corfun(cc_seedling_est)),
                   data.frame(type="seed",value=corfun(cc_seed_est)),
                   data.frame(type="env",value=corfun(cc_env_est))))
}


tfun <- function(m) {
  if (is.matrix(m)) 
    cbind(exp(m[,1]),plogis(m[,2]))
  else c(exp(m[1]),plogis(m[2]))
}

## generate ENVELOPES
rancors <- function(apvars_seed,apvars_seedling,model_seed,model_seedling,
                    ...,nsim=500,calc_range=FALSE) {
  n <- round((sqrt(4*length(apvars_seed)+1)-1)/2)
  v_seed <- matrix(apvars_seed[1:n^2],nrow=n)[1:2,1:2]
  v_seedling <- matrix(apvars_seedling[1:n^2],nrow=n)[1:2,1:2]
  p_seed <- apvars_seed[n^2+(1:2)]
  p_seedling <- apvars_seedling[n^2+(1:2)]
  mseed <- tfun(rbind(p_seed,MASS::mvrnorm(500,mu=p_seed,Sigma=v_seed)))
  mseedling <- tfun(rbind(p_seedling,MASS::mvrnorm(500,mu=p_seedling,Sigma=v_seedling)))
  ## cat(mseed[1,],model_seed,"\n")
  ## cat(mseedling[1,],model_seedling,"\n")
  mc <- cbind(mseed,mseedling)
  colnames(mc) <- NULL ## argh
  ff <- function(x,...) {
    gen_cors(data.frame(range=x[1],nugget=x[2],spmodel=model_seed),
             data.frame(range=x[3],nugget=x[4],spmodel=model_seedling),...)$value
  }
  z <- apply(mc,1,ff,...)
  ## browser()
  if (!calc_range) {
    zz <- t(apply(z,1,quantile,c(0.025,0.975),na.rm=TRUE))
    colnames(zz) <- c("lo","hi")
  } else {
    x <- mc[1,]
    g0 <- gen_cors(data.frame(range=x[1],nugget=x[2],spmodel=model_seed),
             data.frame(range=x[3],nugget=x[4],spmodel=model_seedling),...)
    ## z <- z[,!apply(is.na(z),2,all)]
    mm <- with(g0,data.frame(w,type,z))
    tmpf <- function(zz) {
      zz <- as.matrix(zz[,-(1:2)])
      zzz <- apply(zz,2,get_cpars)
    }
    ##ddply(mm,c("est","seedling"),
  }
  zz
}

rfglsmatch <- list(rf=c("exponential","rational","gaussian","spherical","linear"),
                   gls=c("Exp","Ratio","Gaus","Spher","Lin"))
  
## correspondence between 
cf <- function(rvec,var,nugget,range,m="exponential") {
  CovarianceFct(rvec,model=m,
                param=c(mean=0,var=var*(1-nugget),nugget=var*nugget,range=range))
}

cf2 <- function(rvec,var,nugget,range,m="corExp") {
  var*calc_cor(m,c(range,nugget),rvec,nugget=TRUE)
}

## FIXME: 'true' env might not be exponential ...
true_cors <- function(L,Nvar,Nbar,Ebar,Evar.true,
                      ss_env,rvec) {
  if (!do_poisson) {
    ## FIXME: make sure to choose correct model (not always exponential?)
    cc_seed_true <- with(L$seedparms,
                         cf(rvec,v,n,r)) ## /(v+n))
      ## divide by v [or cc_seed_true[1]] to get corr fn
      cc_env_true <- with(L$seedlingparms,
                          cf(rvec,v,n,r))
      cc_seedling_true <- L$seedparms$m^2*cc_env_true+L$seedlingparms$m^2*cc_seed_true
    } else {
      cc_seed_true <- with(L$seedparms,
                           cf(rvec,Nvar,n,r))
      cc_env_true <- cf(rvec,Evar.true,ss_env["nugget"],ss_env["range"])
      cc_seedling_true <- Nbar^2*cc_env_true+Ebar^2*cc_seed_true
      ## divide by v [or cc_seed_true[1]] to get corr fn
    }
  corfun <- function(x) x/x[1]
  data.frame(w="true",
             rbind(data.frame(type="seedling",value=corfun(cc_seedling_true)),
                   data.frame(type="seed",value=corfun(cc_seed_true)),
                   data.frame(type="env",value=corfun(cc_env_true))))
}

cv2 <- function(x) var(x)/mean(x)^2


get_moments <- function(L) {
  Nbar.est <- mean(L$seeds$count)  ## what about trap area??
  Nvar.est <- var(L$seeds$count)  ## what about trap area??
  Ebar.est <- mean(L$seedlings$count)/Nbar.est  ## what about trap area??
  Svar.est <- var(L$seedlings$count)
  Sbar.est <- mean(L$seedlings$count)
  Ebar.true <- if (is.null(L$seedlings$env)) NULL else mean(L$seedlings$env)

  list(Nbar=Nbar.est,Nvar=Nvar.est,Ebar=Ebar.est,Ebar.true=Ebar.true,Svar=Svar.est,Sbar=Sbar.est)
}


get_all_cors <- function(L,best_seedling=NULL,best_seed=NULL,
                         do_sim=TRUE,do_poisson=FALSE,
                         dAIC_cut=0.01,
                         rvec=seq(0,25,by=0.1)) {
  ## tmpf <- function(m) {
  ##   if (is.list(m)) {
  ##     m <- do.call(rbind,m)
  ##   } else if (is.matrix(m) && ncol(m)!=12) {
  ##     m <- matrix(m,ncol=12,byrow=TRUE)
  ##   }
  ##   m
  ## }

  ## L$seed_anal$pmat <- tmpf(L$seed_anal$pmat)
  ## L$seedling_anal$pmat <- tmpf(L$seedling_anal$pmat)
  ## ## DON'T FIX IF NOT BROKEN -- WILL BREAK IT
  ## ## L$seedling_anal$pmat <- matrix(L$seedling_anal$pmat,nrow=36,byrow=TRUE)

  if (is.null(best_seedling)) best_seedling <- get_best(L$seedling_anal,dAIC_cut=dAIC_cut)
  ## best_seedling <- 9 ## override
  if (is.null(best_seed)) best_seed <- get_best(L$seed_anal,dAIC_cut=dAIC_cut)
  ## ss_seed <- unlist(L$seed_anal$rframe[best_seed,c("range","nugget")])
  ss_seed <- L$seed_anal$rframe[best_seed,c("range","nugget","spmodel")]

  ss_seedling <- L$seedling_anal$rframe[best_seedling,c("range","nugget","spmodel")]
  apvars_seed <- L$seed_anal$pmat[best_seed,]
  apvars_seedling <- L$seedling_anal$pmat[best_seedling,]


  mm <- get_moments(L)

  if (do_sim) {
    if (do_poisson) {
      best_env <- get_best(L$env_anal)
      ss_env <- unlist(L$env_anal$rframe[best_env,c("range","nugget")])
      if (!is.null(L$seedlings$env)) {
        Evar.true <- var(L$seedlings$env)
        truecors <- with(mm,data.frame(rvec,true_cors(L,Nvar,Nbar,Ebar,
                                                      Evar.true,
                                                      ss_env,rvec=rvec),lo=NA,hi=NA))
      }
    } else {
      Evar.true <- L$seedlingparms$v
      truecors <- with(mm,data.frame(rvec,true_cors(L,Nvar,Nbar,Ebar,
                                                    Evar.true,
                                                      ss_env,rvec=rvec),lo=NA,hi=NA))
    }
  } else truecors <- NULL
  
  

  with(mm,
       rbind(truecors,
             data.frame(rvec,gen_cors(ss_seed,ss_seedling,
                                      Nvar=Nvar,Svar=Svar,Nbar=Nbar,Ebar=Ebar,
                                      rvec=rvec),
                        rancors(apvars_seed,apvars_seedling,
                                as.character(unlist(ss_seed["spmodel"])),
                                as.character(unlist(ss_seedling["spmodel"])),
                                Nvar=Nvar,Svar=Svar,
                                Nbar=Nbar,Ebar=Ebar,rvec=rvec))))

}

gplot <- function(m,ribbon=TRUE,grp=FALSE,alpha=0.5) {
  g0 <- ggplot(m,aes(x=rvec,y=value,
                     colour=type,
                     fill=type,
                     linetype=w,
                     ymin=lo,
                     ymax=hi))+
                       theme_bw(base_size=16)+
                         theme(panel.margin=unit(0,"lines"))+
                           labs(x="Distance (m)",y="Autocorrelation")+
                             scale_colour_brewer()+scale_fill_brewer(palette="Dark2")
  if (grp) {
    g0 <- g0 + geom_line(aes(group=grp))
    if (ribbon) g0 <- g0 +geom_ribbon(aes(group=grp),alpha=alpha,colour=NA)
  } else {
    g0 <- g0 + geom_line()
    if (ribbon) g0 <- g0 +geom_ribbon(alpha=alpha,colour=NA)
  }
  g0
}
