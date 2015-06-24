require(deSolve)

## change values within a list or vector
transform.list <- function(`_data`,...) {
    e <- eval(substitute(list(...)), `_data`, parent.frame())
    tags <- names(e)
    inx <- match(tags, names(`_data`))
    matched <- !is.na(inx)
    if (any(matched)) {
        `_data`[inx[matched]] <- e[matched]
    }
    res <- if (!all(matched)) {
        c(`_data`, e[!matched])
    } else `_data`
    ## preserve original class
    class(res) <- class(`_data`)
    res
}

## substitute values into a list
Lsub <- function(L,subs) {
    n <- names(subs)
    for (n in names(subs))
        L[[n]] <- subs[[n]]
    L
}

transform.WNVparams <- transform.list

transform.numeric <- function(x, ...) {
    unlist(eval(transform(as.list(x),...),parent.frame()))
}

## single-strain Jacobian
jacobian1 <- function (params) {
    m <- with (as.list (params),{
        ## a 3 by 3 matrix:
        matrix (c(-mu-v,
                  betaM*vhratio,  ## change to betaM*vhratio
                  0,
                  0,-Dm-k, k,
                  betaB, 0,-Dm),
                ncol=3,
                dimnames=list(c("I","E","B"),c("I","E","B")))
    }
               )
    return(m)
}

## convert a single-strain list of parameters into
## an n-strain list of parameters, all equal
nstrainParams <- function(params,n) {
    np <- transform(params,
              betaM=rep(betaM,n),
              betaB=rep(betaB,n),
              k=rep(k,n),
              v=rep(v,n),
              mu=rep(mu,n),
              nstrain=n)
    class(np) <- "WNVparams"
    np
}


## extract parameters for a particular species
getParams <- function(params,n) {
    transform(params,
              betaM=betaM[n],
              betaB=betaB[n],
              mu=mu[n],
              v=v[n],
              k=k[n]
              )
}

## r value for single-strain parameters
rvalue_0 <- function (params) {
    e <- eigen(jacobian1(params))
    max(Re(e$values))
}

## r value for strain n from a multi-strain set of parameters
rvalue <- function(params,n=1) {
    rvalue_0(getParams(params,n))
}

## convert parameters to input parameters for functions that
## need {betaM, betaB, v, mu} rather than {a1,a2,b,c,casemort,infper}
param_conv2 <- function(L,base=params_robinX) {
    base2 <- do.call(transform.list,c(list(base),L))  ## substitute new values
    transform(base2,
              betaM=a1*a2*b,
              betaB=a1*a2*c,
              v=(1-casemort)/infper,
              mu=casemort/infper)
}

## note: code for producing r contour grid is at the end
## of virdat.rmd

value_grid <- function(clearancevec, ## vector of values of 1/(mu+v)
                       bvec, ## vector of values of betaB
                       baseParams=transform(params_robinX,vhratio=400),
                       which="r")
                        ## bump up R0 of base case to ~ 10
    
{
    results <- matrix(nrow=length(clearancevec),ncol=length(bvec))
    for (i in 1:length(clearancevec)) {
        for (j in 1:length(bvec)) {
            ## calculate appropriate r value here
            params <- param_conv2(list(infper=1/clearancevec[i],b=bvec[j]),
                                  base=baseParams)
            results[i,j] <- summary(params)[[which]]
        }
    }
    L <- list(clearance=clearancevec,beta=bvec,results)
    names(L)[3] <- which
    L
}


rvalue_grid <- function(...) value_grid(...,which="r")
r0value_grid <- function(...) value_grid(...,which="r0")

if (FALSE) {
    ## TESTING
    clearancevec <- seq(0.1,1,length=51) ## infectious period from 10 days to 1 day
    ## we want to specify b (bird -> mosq transm probability) and transform
    ##   that to betaB (a1*a2*b = bites/day*fraction of bites*transm prob)
    bvec <- seq(0,0.7,length=51)
    rgrid <- rvalue_grid(clearancevec,bvec)
    r0grid <- r0value_grid(clearancevec,bvec)
    par(cex=3)
    with(rgrid,contour(clearance,beta,r,labcex=2))
    with(r0grid,contour(clearance,beta,r0,labcex=2,levels=c(1,2,5,10),
                        col="red",lty=2,add=TRUE))
    with(r0grid,contour(clearance,beta,r0,labcex=2))
}
    
## r0 value for single-strain parameters: n.b. UNSCALED (i.e. using k)
r0value_0 <- function (params,type="classical") {
    r0 <- with(as.list (params),
     {
         (betaB*betaM*vhratio*k)/
              (Dm*(mu+v)*(k+Dm))
     })
    if (type=="classical") r0 else sqrt(r0)
}

## r0 value for strain n
r0value <- function(params,n=1,...) {
    r0value_0(getParams(params,n),...)
}

## baseline parameter table from Wonham et al.
wonhamPar1 <- matrix(c(0.09,0.03,0.16,     ## biting rate
                       0.88,0.8,1,         ## M->B transm
                       0.16,0.02,0.24,     ## B->M transm
                       0.029,0.016,0.070,  ## adult mort
                       ## AMK: increase to 0.15
                       ## 0.15,0.016,0.070,
                       0.106,0.087,0.125,  ## M exp->inf
                       0.143,0.125,0.2,    ## B mort
                       0,NA,NA),           ## B recovery
                     ncol=3,
                     dimnames=list(c("a","b","c",
                     "muA","k","muV","g"),
                     c("est","lwr","upr")),
                     byrow=TRUE)
descrip <- c("M per capita biting rate on crows",
             "WN transmission probability, M to B",
             "WN transmission probability, B to M",
             "M adult per capita mortality rate",
             "M per capita transition rate, exposed to infected",
             "B per capita mortality rate from WN",
             "B per capita recovery rate from WN")

## our baseline parameters: renaming and collapsing parameters
## we're not changing anything about Wonham's parameters here,
## just adapting it to our formulation (we add a tiny recovery rate)
params0L <-
    with(as.list(wonhamPar1[,"est"]),
         list(betaB=a*b,    ## M to B transmission
              Dm=muA,       ## adult M death rate
              v=1e-6,       ## B recovery rate (make >0 to avoid errors)
              mu=muV,       ## B disease-induced mortality
              betaM=a*c,    ## B to M transmission
              k=k,          ## M exposed -> infected rate
              vhratio=1,    ## vector/host ratio
              nstrain=1))   ## number of strains

class(params0L) <- "WNVparams"

latex_pars <- c(betaB="\\beta_B",Dm="D_m",
                v="\nu",mu="\\mu_V",betaM="\\beta_M",
                k="k",
                vhratio="V/H",
                nstrain="n")
latex_pars[] <- paste0("$",latex_pars,"$")

fparams <- function(x,latex=FALSE) {
    m <- do.call(rbind,x)
    colnames(m) <- paste0("sp",seq(ncol(m)))
    if (latex) {
        rownames(m) <- latex_pars[rownames(m)]
    }
    m
}

prtab <- function(p,unique=FALSE,...) {
    require("Hmisc")
    ff <- fparams(p,latex=TRUE)
    if (unique) ff <- ff[apply(ff,1,function(x) x[1]!=x[2]),]
    latex(ff,file="",digits=2,...)
}

print.WNVparams <- function(x,digits=3,horiz=FALSE,...) {
    if (!horiz) {
        print(fparams(x),digits=digits)
    } else {
        print(t(fparams(x)),digits=digits)
    }
}


## create a 'skeleton' -- a list with the correct shape
pskelFun <- function(nsp) {
    transform(params0L,
              betaB=numeric(nsp),
              betaM=numeric(nsp),
              k=numeric(nsp),
              mu=numeric(nsp),
              v=numeric(nsp))
}

## create a 'skeleton'
xskelFun <- function(nsp) {
    list(S=numeric(1),
         I=numeric(nsp),
         R=numeric(nsp),
         X=numeric(nsp),
         A=numeric(1),
         E=numeric(nsp),
         B=numeric(nsp))
}
         
paramList <- function(params,nstrain=NULL) {
    if (is.null(nstrain)) nstrain <- params[["nstrain"]]
    relist(params,pskelFun(nstrain))
}

summary.WNVparams <- function(object,full=TRUE,...) {
    n <- object[["nstrain"]]
    r <- sapply(seq(n),rvalue,params=unclass(object))
    r0 <- unname(sapply(seq(n),r0value,params=unclass(object)))
    ## FIXME: make a nice summary function that also prints
    ## out the case mortality and infectious period
    res <- list(r=r,r0=r0)
    ## if (full
    class(res) <- "summary.WNVparams"
    res
}

print.summary.WNVparams <- function(x,digits=3,...) {
    op <- options(digits=digits)
    on.exit(options(op))
    cat("r:",x$r,"\n",sep=" ")
    cat("R0:",x$r0,"\n",sep=" ")   
}

scale_R0 <- function(params,sc,n) {
    ## change R0 of species n by an amount sc
    ## n.b.: avoid collision between this name and the elements of params!
    v0 <- rep(1,params[["nstrain"]])
    ## OR: params[["betaM"]][n] <- params[["betaM"]][n]*sc
    v0[n] <- sc
    transform(params,betaM=betaM*v0)
}

## scale r of strain n by a specified amount
scale_r <- function(params,sc,n) {
    v0 <- rep(1,params[["nstrain"]])
    v0[n] <- sc
    transform(params,
              betaM=betaM*v0,
              mu=mu*v0,
              v=v*v0)
}

## show the difference between the scaled r and the target r
r_targetfun <- function(sc,params,n,target) {
    rvalue(scale_r(params,sc,n),n)-target
}

## find scaling factor to achieve target r for strain n;
##  return NA if unsuccessful
find_rsc <- function(params,target,n,interval=c(0,10)) {
    uu <- try(uniroot(r_targetfun,interval,params=params,n=n,target=target),
              silent=TRUE)
    if (inherits(uu,"try-error")) NA else uu$root
}

## take baseline parameters 'params'; set r value to 'r' for strain n
set_r <- function(params,r,n,...) {
    sc <- find_rsc(params,r,n,...)
    scale_r(params,sc,n)
}

w_grad <- function (t,x,params){
    n <- params[["nstrain"]]
    Params <- relist(params,pskelFun(n))
    X <- relist(x,xskelFun(n))
    ## single-strain model
    with(as.list (c(params,X)),
     {
         ## {S,I,R,X} = birds {susceptible, infected, recovered, dead}
         ## {A,E,B} = mosq {larvae, susc, exposed, infected}

         hostinf <- (betaB*B*S)/(S+I+R)
         ##birds
         dS <- -sum(hostinf)
         dI <-  hostinf - (mu+v)*I
         dR <-  v*I
         dX <-  mu*I

         minf <- (betaM*A*I)/(S+I+R)

         ## need birth rate=per capita *N  (!!)
         dA <- Dm*(A+sum(E)+sum(B))-sum(minf) - Dm*A  
         dE <- minf - Dm*E - k*E
         dB <- k*E - Dm*B

         ## order matters here!
         ## don't know if we can make it not matter by
         ##   naming this vector appropriately
         dx <- c(dS, dI, dR, dX, dA, dE, dB, hostinf)
         list(dx)}
         )
}



startfun <- function(params,nstrain=params[["nstrain"]],
                     startfrac=c(M=5e-5,B=0),N=1) {
    ## N: bird density

    mknvec <- function(x,base,svec=seq(nstrain)) {
        setNames(x,paste0(base,svec))
    }
    Aval <- with(as.list(params),vhratio*N)
    Bstart <- rep(startfrac["B"],length=nstrain)
    Mstart <- rep(startfrac["M"],length=nstrain)
    c(S=N*(1-sum(Bstart)),
      mknvec(N*Bstart,"I"),
      mknvec(rep(0,nstrain),"R"),
      mknvec(rep(0,nstrain),"X"),
      A=Aval*(1-sum(Mstart)),
      mknvec(rep(0,nstrain),"E"),
      mknvec(Aval*Mstart,"B"),
      mknvec(rep(0,nstrain),"hostinf"))
}


## robin parameters: sensible/adjusted for AMK advice, but V/H ratio
##  is still 1 ...
adult_mort <- 0.15
params_robin1 <- transform(params0L,
                           betaB=0.09, ## should be proportion of robins in pop;
                                     ## following W+, look up CBC data on 
                                     ## prop of AMRO in NY area, divide by
                                     ## 3 (1/3 bite per female per day)
                                     ## multiply by M->B transmission, 
                                     ## ~1 according to AMK
                           ## birth_rate,
                           Dm=adult_mort,  ## AMK says 0.1-0.2, 
                                     ## see jones_rainfall_2012
                           ## modify mu and v to keep infectious period
                           ## mu/(v+mu)=2/17=case mortality rate
                           ## mu+v (=1/inf period) stays fixed at old value
                           ## mu1/(v1+mu1)=2/17
                           ## mu0+v0=mu1+v1
                           ## v0 ~ 0 (we only made it >0 to avoid problems)
                           ## mu0 = mu1+v1
                           ## v1=mu0-mu1
                           ## mu1/mu0=2/17
                           ##
                           v=mu*(1-2/17),  ## to be safe, change v first
                           mu=mu*2/17,
                           betaM=betaM*1/0.88 )
                                      ## suppose B->M is identical for robins
                                     ## then betaM changes due to change in
                                     ## encounter rate = W+ 'a' parameter
                                     ## prop to (robin density/crow density)
                                     ## leave it alone for now

## transform parameters
## a1=prop of bites on focal species
## a2= bites per mosquito per day
## b=M to B transmission prob
## c=B to M transmission prob
params_robinX <- c(list(a1=0.5,a2=1/8,
                        b=1,c=1),
                   params_robin1)
class(params_robinX) <- "WNVparams"
## the same but with casemort and infper added
params_robinX <- transform(params_robinX,
                           casemort=mu/(mu+v),
                           infper=1/(mu+v))


runSim <- function(params,tvec=0:120,...) {
    res <- ode(startfun(params,...),tvec,w_grad,params)
    class(res) <- c("WNVsim",class(res))
    res
}
plot.WNVsim <- function(x,y,type="base",
                        vars=c("S","I1","B1"),
                        colvec=c(1,2,4,5,6,3),
                        logy=TRUE,mosqprev=TRUE,...) {
    ## need y for generic compatibility
    if (type=="base") {
        op <- par(las=1,bty="l")
        on.exit(par(op))
        if (mosqprev) {
            ## this will FAIL for a multistrain model!
            if ("E2" %in% names(x)) {
                stop("mosquito rescaling doesn't work for multi-strain models")
            }
            ## rescale mosquitoes by density
            x <- plyr::mutate(as.data.frame(x),
                              Nmosq=A+E1+B1,
                              A=A/Nmosq,
                              B1=B1/Nmosq,
                              E1=E1/Nmosq)
        }
        log <- if (logy) "y" else ""
        matplot(x[,"time"],x[,vars],log=log,
                type="l",lty=1,col=colvec,...)
    } else {
        stop("unfinished")
        require(ggplot2)
    }
}

calc_obj <- function(x) {
    prev <- x[,"I1"]
    peakval <- which.max(prev)
    meanval <- mean(prev)
    list(peakval=peakval,meanval=meanval)
}

objfun3 <- function(pars,targets=c(meanprev=0.22,peakday=90,peakscale=100),
                    debug=TRUE) {
    pp <- pconv2(as.list(pars))
    tmprun <- ode(startfun(pp),
                  seq(0,150),w_grad,pp)
    obj <- with(c(as.list(targets),calc_obj(tmprun)),
                ((peakval-peakday)/peakscale)^2+(meanval-meanprev)^2)
    debug <- as.numeric(debug)
    if (debug>0) {
        if (debug==1) cat(pars,obj,"\n")
        else {
            dstr <- paste(names(pars),round(pars,3),sep="=",collapse=" ")
            objvec <- round(c(prev=meanval,peak=peakval,obj=obj),3)
            dstr <- paste(dstr,
                          paste(names(objvec),objvec,sep="=",collapse=" "))
            cat(dstr,"\n")
        }
    }
    obj
}

## these are what we think the RANGES of the parameters could/should be for robins
newrobinpars <- list(a1=c(0.5,0.3,0.8), ## prop of bites on focal species
                a2=c(1/8,1/10,1/5),     ## bites per mosq per day
                b=c(0.99,0,1),          ## mosquito to bird
                c=c(0.99,0,1),          ## bird to mosquito
                Dm=c(0.15,0.1,0.2),     ## adult mosquito death
                casemort=c(2/17,1/17,2/21),  ## bird case mort
                infper=c(4,3,5),             ## infectious period
                vhratio=c(100,1,3000))
newrobinval <- sapply(newrobinpars,"[[",1)
newrobinlwr <- sapply(newrobinpars,"[[",2)
newrobinupr <- sapply(newrobinpars,"[[",3)


## simple SIR stuff
convparm <- function(par) {
   ## R0=beta/gamma; r=beta-gamma
    ## beta = r + gamma
    ## R0 = (r+gamma)/gamma = r/gamma+1
    ## R0-1 = r/gamma
    ## gamma = r/(R0-1)
    ## beta = r+r/(R0-1)
     with(as.list(par),
         c(gamma=r/(R0-1),
           beta=r+r/(R0-1)))
}
    
SIRgrad <- function(t,y,params) {
    g <- with(as.list(c(y,convparm(params))),c(S=-beta*S*I,
                                               I=beta*S*I-gamma*I))
    list(g,NULL)
}

objfun4 <- function(pars,targets=c(meanprev=0.22,peakday=90,peakscale=100),debug=TRUE) {
    tmprun <- ode(c(S=0.95,I=0.05),
                  seq(0,150),SIRgrad,pars)
    prev <- tmprun[,"I"]
    peakval <- which.max(prev)
    meanval <- mean(prev)
    obj <- with(as.list(targets),
                ((peakval-peakday)/peakscale)^2+(meanval-meanprev)^2)
    obj
}


## generic plotting functions for ODE runs
plotfun <- function(run,type=c("base","ggplot")) {
    type <- match.arg(type)
    rr <- as.data.frame(run)
    if (type=="base") {
        par(mfrow=c(1,3),las=1,bty="l")
        colvec <- c("gray","red","blue")
        ltyvec <- c(3,1,2)
        ## susc & inf densities
        with(rr,matplot(time,cbind(S,I1,I2),
                        type="l",col=colvec,
                        ylab="Density",
                        lty=ltyvec,
                        main="S and I prevalence"))
        legend("topright",col=colvec,lty=ltyvec,
               c("S","I1","I2"),bty="n")
        ## str 1 vs 2
        with(rr,plot(time,I1/(I1+I2),type="l",ylim=c(0,1),
                     main="rel. prevalence of strain 1"))
        abline(h=0.5,col=adjustcolor("gray",alpha=0.5))
        with(rr,matplot(time, cbind(X1,X2), log="y", type="l", col=colvec[-1], 
                        lty=ltyvec[-1],ylab="Density",
                        main="cumulative infection"))
        legend("topright",col=colvec[-1],lty=ltyvec[-1],
               c("X1","X2"),bty="n")
    } else {
        require(ggplot2)
        require(reshape2)
        rr <- transform(rr,prev.frac=I1/(I1+I2))
        mm <- melt(subset(rr,select=c(time,I1,I2)),id.var="time")
        prevplot <- ggplot(mm,aes(time,value,colour=variable))+
            geom_line()+scale_colour_brewer(palette="Set1")+
                ylab("density")
        pfracplot <- ggplot(rr,aes(time,prev.frac))+
            geom_line()+ ylab("fraction of strain 1")+
                geom_hline(yintercept=0.5,lty=2)
        return(list(prev=prevplot,pfrac=pfracplot))
    }
}

plotfun2 <- function(run, mu=c(0,0), v=c(0,0)) {
    par(mfrow=c(1,2),las=1,bty="l")
    colvec <- c("red","blue")
    ltyvec <- c(3,1,2)
    rr <- as.data.frame(run)
    ## susc & inf densities
    with(rr,matplot(time,cbind(I1,I2),
                    log="y",type="l",col=colvec,
                    ylab="Density",
                    lty=ltyvec,
                    main="I prevalence"))
    legend("topright",col=colvec,lty=ltyvec,
           c("I1","I2"),bty="n")
    ## str 1 vs 2
    with(rr,plot(time,(I1*mu[1])/(I1+I2)+(I2*mu[2])/(I1+I2),type="l",ylim=c(0,1),
         main="rel. prevalence of strain 1"))
    abline(h=0.5,col=adjustcolor("gray",alpha=0.5))
    with(rr,matplot(time, cbind(X1,X2), log="y", type="l", col=colvec[-1], 
                    lty=ltyvec[-1],ylab="Density",
                    main="cumulative infection"))
    legend("topright",col=colvec[-1],lty=ltyvec[-1],
           c("X1","X2"),bty="n")
}
