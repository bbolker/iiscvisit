require(deSolve)
## functions for transmission and deriv. of transmission
##   as a function of virulence (alpha)
beta = function(alpha,c0,gamma=2) { c0*alpha^(1/gamma) }
derivbeta = function(alpha,c0,gamma=2) { c0/gamma*alpha^(1/gamma-1) }

## uh-oh, this may be too old to reconstruct from scratch ...
if (file.exists("myxoz.RData")) load("myxoz.RData")

## override/hack old myxoz() function 
myxoz2 <- function (x, deriv = 0) 
{
    deriv <- as.integer(deriv)
    if (deriv < 0 || deriv > 3) 
        stop("'deriv' must be between 0 and 3")
    if (deriv > 0) {
        z0 <- double(z$n)
        z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
            z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
            c = z0), list(y = 6 * z$d, b = z0, c = z0))
        z[["d"]] <- z0
    }
    res <- stats:::.splinefun(x,z)
        if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
        res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
    res
    ## .C("spline_eval", z$method, as.integer(length(x)), x = as.double(x), 
    ## y = double(length(x)), z$n, z$x, z$y, z$b, z$c, z$d, 
    ## PACKAGE = "stats")$y
}
environment(myxoz2) <- environment(myxoz)
myxoz <- myxoz2

## derivatives of ($S$, $I$, $\bar\alpha$)

derivfun1 = function(t,y,parms) {
  derivs = with(c(as.list(y),as.list(parms)),
    c(m*(N-S)-beta(alpha,c0,gamma)*S*I,
      (beta(alpha,c0,gamma)*S-(m+alpha))*I,
      h*(S*derivbeta(alpha,c0,gamma)-1)))
  list(derivs,NULL)
}

## log-I
derivfun1L = function(t,y,parms) {
  derivs = with(c(as.list(y),as.list(parms)),
    c(m*(N-S)-beta(alpha,c0,gamma)*S*exp(logI),
      (beta(alpha,c0,gamma)*S-(m+alpha)),
      h*(S*derivbeta(alpha,c0,gamma)-1)))
  ## cat(t,y,derivs,"\n")
  ## if (any(!is.finite(derivs))) browser()
  list(derivs,NULL)
}

## no birth
derivfun2L = function(t,y,parms) {
  derivs = with(c(as.list(y),as.list(parms)),
    c(-m*S-beta(alpha,c0,gamma)*S*exp(logI),
      (beta(alpha,c0,gamma)*S-(m+alpha)),
      h*(S*derivbeta(alpha,c0,gamma)-1)))
  ## cat(t,y,derivs,"\n")
  ## if (any(!is.finite(derivs))) browser()
  list(derivs,NULL)
}

## log-I and log-alpha
derivfun1L2 = function(t,y,parms) {
  derivs = with(c(as.list(y),as.list(parms)),
    {alpha=exp(logalpha)
     c(m*(N-S)-beta(alpha,c0,gamma)*S*exp(logI),
      (beta(alpha,c0,gamma)*S-(m+alpha)),
      h/alpha*(S*derivbeta(alpha,c0,gamma)-1))})
  ## cat(t,y,derivs,"\n")
  ## if (any(!is.finite(derivs))) browser()
  list(derivs,NULL)
}

casemortfun <- function(T) {
  pmax(0,pmin(1,exp((334.7-T)/70.3)/100))
}

## where do these come from??
myxoparms1 <- c(r=0.008,m=0.006,c0=9.38*8e-4*250,
  gamma=1.22,h=10,N=1)
## high-alpha start
myxostart1 = c(S=0.99,I=0.01,R=0,alpha=0.15)
myxoderiv1 = function(t,y,parms) {
  derivs = with(c(as.list(y),as.list(parms)),
    c(r*(N-I)-m*S-beta(alpha,c0,gamma)*S*I,
      beta(alpha,c0,gamma)*S*I-(m+alpha)*I,
      alpha*I*(1-casemortfun(1/alpha)),
      h*(S*derivbeta(alpha,c0,gamma)-1)))
  list(derivs,NULL)
}

## use "realistic" average-beta curve
myxoderiv2 = function(t,y,parms) {
  derivs = with(c(as.list(y),as.list(parms)),
    c(r*(N-I)-m*S-z(alpha)*S*I,
      myxoz(alpha)*S*I-(m+alpha)*I,
      alpha*I*(1-casemortfun(1/alpha)),
      h*(S*myxoz(alpha,1)-1)))
  list(derivs,NULL)
}

## use "realistic" average-beta curve
myxoderiv3 = function(t,y,parms) {
  derivs = with(c(as.list(y),as.list(parms)),
    c(r*(N-I)-m*S-myxoz(alpha)*S*I,
      (myxoz(alpha)*S-(m+alpha))*I,
      ## alpha*I*(1-casemortfun(1/alpha)),
      h*(S*myxoz(alpha,1)-1)))
  list(derivs,NULL)
}

## default parameters, starting values:
params1 = c(h=1,c0=3,m=1.5,N=1,gamma=2)
startval1 = c(S=0.9,I=0.1,alpha=1.5)

## more extreme
startval2 = c(S=0.99,I=0.01)
params2 = c(h=5,c0=3,m=1.5,N=1,alpha=1.5,gamma=2)

gfun = function(g,m=1) {m^(1/g-1)*(g-1)^(1-1/g)/g }

eqvals = function(parms) {
  with(as.list(parms),
       { R0=c0*N*gfun(gamma,m)
         beta=c0*(m/(gamma-1))^(1/gamma);
         c(S=N/R0,I=(m/beta)*(R0-1),
           R0=R0,alpha=m/(gamma-1),beta=beta)
       })
}

## compute max virulence/time to max virulence/initial dV/dt
findvals = function(
  params, ## parameters
  start,  ## vector of starting S0,I0  
  I0,     ## alternative: specify I0 only (S0=1-I0)
  times=seq(0,20,by=0.025), ## initial set of times to run
  scale=TRUE,  ## scale peak vir by starting vir?
  scaletime=TRUE, ## scale peak time?
  logI=FALSE,     ## use deriv on log-I scale?
  logalpha=FALSE, ## use deriv on log-alpha scale?
  births=TRUE,    ## SIR model with or without births?
  debug=FALSE,    ## return just time series
  S0,             ## specify S0
  search=TRUE,maxit=20,
  plot.it=FALSE,col2=2,
  xlab="Time",
  ylab1="Infection prevalence (I)",
  ylab2=expression(paste("Virulence ",(alpha)))) {
  if (!logI && !births) stop("must use log-I with no-births option")
  eq = eqvals(params)
  if (missing(start)) {
    if (missing(I0)) stop("must specify either start or I0")
    ## start = c(eq["S"]*(1-I0),eq["I"]*I0,eq["alpha"])
    if (missing(S0)) S0 = 1-I0
    start = c(S=params["N"]*S0,I=params["N"]*I0,alpha=eq["alpha"])
    names(start)=c("S","I","alpha") ## argh
  } else I0=start["I"]/eq["I"]
  start0=start
  dfun = if(!births) derivfun2L else if (!logI) derivfun1 else if (!logalpha) derivfun1L else derivfun1L2
  if (logI) { start["I"] = log(start["I"])
              names(start)[names(start)=="I"] = "logI" }
  if (logalpha) { start["alpha"] = log(start["alpha"])
                  names(start)[names(start)=="alpha"] = "logalpha"
                }
  runsim = function(times) {
    OK = FALSE
    it = 1
    while (!OK && it<maxit) {
      L = lsoda(y=start,times=times,parms=params,func=dfun)
      wa = grep("alpha",colnames(L))
      wi = grep("I",colnames(L))
      ## KLUGE!! should figure out why these crash in this way ...
      if (any(is.na(L))) {
        last.OK = min(row(L)[is.na(L)])
        last.time = times[last.OK]
        times = seq(min(times),last.time*0.9,length=length(times))
        cat("NAs encountered, rerunning (max time",last.time,")\n")
      } else if (search) {
        OK = !(all(diff(L[,wa])>0) || all(diff(L[,wi])>0))
        ## might be better to deal with inf peak and alpha peak
        ##   separately, but too complicated ... worried a bit
        ## that running for long times to get inf peak will miss
        ## alpha peak
        ## could also try to get initial order-of-magnitude estimate
        ##  of when inf and alpha peak occurred ...
        if (!OK) {
          times = seq(min(times),2*max(times),length=length(times))
          cat("haven't reached vir/inf peak, doubling time to",max(times),"\n")
        }
      } else OK=TRUE
      if (logI) {
        L[,"logI"] = exp(L[,"logI"])
        colnames(L)[colnames(L)=="logI"]="I"
      }
      if (logalpha) {
        L[,"logalpha"] = exp(L[,"logalpha"])
        colnames(L)[colnames(L)=="logalpha"]="alpha"
      }
      it=it+1
    } ## while !OK
    list(times=times,L=L)
  }
  ret = runsim(times)
  L = ret$L
  times = ret$times
  if (debug) return(L)
##   if (search) {
##     ## re-run for more precise time
##     lasttime = L[which.max(L[,"alpha"])+1,"time"]
##     cat("re-running until",lasttime,"\n")
##     times = seq(min(times),lasttime,length=length(times))
##     ret = runsim(times)
##     L = ret$L
##     times = ret$times
##   }
  maxalpha = max(L[,"alpha"])
  if (scale) maxalpha <- maxalpha/eq["alpha"]
  maxtime = L[which.max(L[,"alpha"]),"time"]
  eqI = if (!births) eqI = 0 else eq["I"]
  if (scaletime) {
    maxtime = L[which.max(L[,"alpha"]),"I"]/max(c(eqI,L[,"I"]))
  }
  initderiv = derivfun1(0,start0,params)[[1]][3]
  if (plot.it) {
    op = par(mar=c(5,4,2,4)+0.1)
    with(as.data.frame(L),
         {
           plot(I~time,type="l",xlab=xlab,ylab=ylab1,
                ylim=c(min(I),max(max(I),eqI)))
           t0 = time[which.max(alpha)]
           i0 = I[which.max(alpha)]
           abline(v=t0,lty=3)
           x0 = par("usr")[1]
           cat(x0,i0,t0,"\n")
           segments(x0,i0,t0,i0,lty=3)
           par(new=TRUE)         
           plot(alpha~time,type="l",axes=FALSE,col=col2,
                ann=FALSE,ylim=c(min(alpha),max(alpha)*1.1))
           axis(side=4,col.axis=col2,col=col2)
           mtext(side=4,line=par("mgp")[1],ylab2,col=col2,
                 cex=par("cex"),las=0)
           abline(h=max(alpha),lty=3,col=col2)
         })
    par(op)
  }
  list(params=params,I0=I0,maxalpha=maxalpha,maxtime=maxtime,
       initderiv=initderiv,eq=eq)
}

nR0 = 30
nalpha = 30

R0trans = function(x,inv=FALSE) {
  if (inv) 10^x+1 else log10(x-1)
}
R0trvec = seq(R0trans(1.1),R0trans(50),length=nR0)
R0vec = R0trans(R0trvec,inv=TRUE)
alphatrvec = seq(0,3,length=nalpha)
alphavec = 10^alphatrvec
