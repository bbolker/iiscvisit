## requires pine-funs.R
run_simanal <- function(do_sims=TRUE,
                        do_seeds=TRUE,
                        do_seedlings=TRUE,
                        do_env=TRUE,
                        seedlingplotmin=50,
                        seed.p=10, ## number of plots
                        seedsamp=100,    ## quadrats per plot
                        seedlingsamp=100,
                        save_full=FALSE,
                        do_poisson=FALSE,
                        transform.fun=NULL,
                        seed=1001,
                        ## PARAMETERS (for Gaussian runs)
                        seedparms=list(m=13,v=95,n=0,r=6.5,s=seedsamp,p=10),
                        seedlingparms=list(m=1/13,
                          v=0.013,n=0,r=2,s=seedlingsamp,p=10),
                        ## gives mean approx. 1, variance approx. 4 (I think) for seedlings
                        ## (for Poisson/binomial:)
                        ## seedparms=list(m=2.36,v=0.4,n=0,r=6.5,s=20,p=10)
                        ## seedlingparms=list(m=-3.5,v=1.5^2,n=0,r=2,s=75,p=10)
                        ## i.e., these should give mean/variance equivalent to the previous
                        ##    gaussian example
                        ## mod <- c("Exp","Exp","Ratio","Spher","Lin","Gaus")
                        mod=c("Exp","Ratio"),
                        nugget=c(TRUE,TRUE),
                        valList=expand.grid(c(2,8,32,64),c(0.001,0.1)),
                        autostart=TRUE,  ## add autostart to list of start vals?
                        vpList=list(pow=varPower()),
                        ## pow1=varPower(value=1), NULL
                        ## 
                        modList=list(plot=count~plot,plotxy=count~plot*(x0+y0)),
                        glscList=list(nlminb=glsControl())
                        ) {

  retlist <- list()
  ## source("pine-funs.R")
  if (!do_sims) {
    load("pine2.RData")
    seedparms <- seedlingparms <- NULL
    retlist <- c(retlist,list(seeds=seeds,seedlings=seedlings))
  } else {
    set.seed(seed)
    if (do_poisson) {
      seed_invlink <- exp
      seed_rfun <- rpois
      seedling_invlink <- plogis
      seedling_rfun <- rbinom
    } else {
      seed_invlink <- identity
      seed_rfun <- ifun
      seedling_invlink <- identity
      seedling_rfun <- prodfun
    }
    ##  seedCV2 <- with(seedparms,v/m^2)
    ## estabCV2 <- with(seedlingparms,v/m^2)
    ## cat(seedCV2,estabCV2,"\n")
    
    ## convenience function for simulating with a specified list of
    ##   parameters (use do.call instead?)
    sfun <- function(retval="sample",parms=seedparms) {
      with(parms,
           do_sim(m=m,v=v,nugget=n,range=r,
                  rfun=seed_rfun,invlink=seed_invlink,
                  retval=retval,sample=s))
    }
    
    ## create N plots
    seedsims_list <- replicate(seed.p,sfun(),simplify=FALSE)
    ## stick seedling info together with plot info, label by plots
    seeds <- labfun(seedsims_list,LETTERS)
    
    ## SIMULATE SEEDLINGS (independent seed inputs)
    seedlingsims_list <- replicate(seedlingparms$p,
                                   with(seedlingparms,
                                        do_sim(m=m,v=v,nugget=n,
                                               range=r,
                                               rfun=seedling_rfun,invlink=seedling_invlink,
                                               retval="sample",sample=s,
                                               input=sfun("grid"))),
                                   simplify=FALSE)
    seedlings <- labfun(seedlingsims_list,LETTERS)
    ## save("seedsims",file=paste(fname,"seedsims.RData",sep="_"))
    ## save("seedlingsims",file=paste(fname,"seedlingsims.RData",sep="_"))
    retlist <- c(retlist,list(seeds=seeds,seedlings=seedlings))
  }
  ## SIMS GENERATED OR DATA LOADED, NOW ANALYZE
  library(nlme)
  ## define grid of models etc. to try
  ## 1. SPATIAL MODELS: corr model +/- nugget
  mod2 <- paste("cor",mod,sep="")
  ## nugget <- c(TRUE,FALSE,rep(TRUE,4))
  model.names <- paste(mod,ifelse(nugget,"N",""),sep="")
  ## 2. STARTING VALUES
  valList <- split(as.matrix(valList),row(valList))
  if (autostart) valList <- c(valList,list(numeric(0))) ## add default/automatic start
  ## 3. use varPower?
  ## vpList <- list(varPower(),NULL)
  ## 4. LINEAR TREND IN FIXED EFFECTS?
  ## 5. OPTIMIZERS?
  ## glscList <- list(glsControl(),glsControl(opt="optim"))
  svnames <- sapply(valList,paste,collapse="/")
  ## mnames <- c("nlminb","optim")
  ## FIXME: these need to match arguments
  mnames <- names(glscList)
  fixmodnames <- names(modList)
  ## vpmodnames <- c("pow","const")
  vpmodnames <- names(vpList)
  ## rframe0 <- expand.grid(spmodel=model.names,startval=svnames,
  ##                        method=mnames,fixmod=fixmodnames,vpmod=vpmodnames)
  rframe0 <- expand.grid(fixmod=fixmodnames,
                         spmodel=model.names,
                         startval=svnames,
                         vpmod=vpmodnames,
                         method=mnames)

  if (!is.null(transform.fun)) {
    seeds$count <- transform.fun(seeds$count)
    seedlings$count <- transform.fun(seedlings$count)
  }
  
  if (do_seeds) {
    cat("SEEDS\n")
      seedstuff <- do_analysis(seeds,modList,mod2,valList,vpList,rframe0,save_full=save_full)
      retlist <- c(retlist,list(seedparms=seedparms,seed_anal=seedstuff))
    }
  
  if (do_seedlings) {
    cat("SEEDLINGS\n")
    if (!do_sims) {
      ## drop plots without enough counts
      totseedlings <- with(seedlings,tapply(count,plot,sum))
      OKplots <- with(seedlings,plot %in% levels(plot)[totseedlings>seedlingplotmin])
      seedlings2 <- gdata::drop.levels(subset(seedlings,OKplots),reorder=FALSE)
      seedlingstuff <- do_analysis(seedlings2,modList,mod2,valList,vpList,rframe0,save_full=save_full)
      retlist$seedlings <- seedlings2
    } else {
      seedlingstuff <- do_analysis(seedlings,modList,mod2,valList,vpList,rframe0)
    }

    if (do_sims && do_env && do_poisson) {
      cat("ENV\n")
      envmodList <- lapply(modList,
                           function(x) { x[[2]] <- quote(env); x })
      envstuff <- do_analysis(seedlings,envmodList,mod2,valList,vpList,rframe0)
    } else envstuff <- NULL
    
    retlist <-c(retlist,list(seedling_anal=seedlingstuff,
                             env_anal=envstuff,modList=modList,mod2=mod2,valList=valList,
                             vpList=vpList,
                             seedlingparms=seedlingparms,do_poisson=do_poisson))
  }
  retlist
}
