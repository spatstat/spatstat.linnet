#'
#'  lurklppm.R
#'
#'  Method for 'lurking' for lppm
#'
#'  $Revision: 1.7 $ $Date: 2025/12/20 04:26:26 $

lurking.lpp <- lurking.lppm <- function(object, covariate,
                         type="raw",
                         cumulative=TRUE,
                         ..., 
                         plot.it=TRUE,
                         plot.sd=is.poisson(object),
                         clipwindow=NULL,
                         rv=NULL,
                         envelope=FALSE, nsim=39, nrank=1,
                         typename,
                         covname, 
                         oldstyle=FALSE, check=TRUE,
                         verbose=TRUE,
                         nx=128,
                         splineargs=list(spar=0.5),
                         internal=NULL) {
  cl <- match.call()
  clenv <- parent.frame()

  if(is.lpp(object)) {
    X <- object
    object <- lppm(X ~ 1)
    dont.complain.about(X)
  } else verifyclass(object, "lppm")
  
  object2D <- as.ppm(object)
  
  ## match type argument
  type <- pickoption("type", type,
                     c(eem="eem",
                       raw="raw",
                       inverse="inverse",
                       pearson="pearson",
                       Pearson="pearson"))
  if(missing(typename))
    typename <- switch(type,
                       eem="exponential energy weights",
                       raw="raw residuals",
                       inverse="inverse-lambda residuals",
                       pearson="Pearson residuals")

  ## default name for covariate
  if(missing(covname) || is.null(covname)) {
    co <- cl$covariate
    covname <- if(is.name(co)) as.character(co) else
               if(is.expression(co)) format(co[[1]]) else 
               if(is.character(co) && 
                  length(co) == 1 &&
                  co %in% c("x", "y")) paste(co, "coordinate") else NULL
  }

  #' spatial covariates
  CovariatesDF <- getppmOriginalCovariates(object2D)

  ## may need to refit the model
  if(is.expression(covariate)) {
    ## expression could involve variables that are not stored in object
    neednames <- all.vars(covariate)
    if(!all(neednames %in% colnames(CovariatesDF))) {
      object <- update(object, save.all.vars=TRUE)
      object2D <- as.ppm(object)
      CovariatesDF <- getppmOriginalCovariates(object2D)
    }
  }

  ## simulation involved?
  Xsim <- NULL
  if(!isFALSE(envelope)) {
    ## compute simulation envelope
    Xsim <- NULL
    if(!isTRUE(envelope)) {
      ## some kind of object
      Y <- envelope
      if(is.list(Y) && all(sapply(Y, is.lpp))) {
        #' Argument 'envelope' is a list of point patterns
        Xsim <- Y
        envelope <- TRUE
      } else if(inherits(Y, "envelope")) {
        #' Argument 'envelope' is an envelope object
        Xsim <- attr(Y, "simpatterns")
        if(is.null(Xsim))
          stop("envelope object does not contain simulated point patterns")
        envelope <- TRUE
      } else stop("Unrecognised format of argument: envelope")
      nXsim <- length(Xsim)
      if(missing(nsim) && (nXsim < nsim)) {
        warning(paste("Only", nXsim, "simulated patterns available"))
        nsim <- nXsim
      }
    }
  }

  
  #################################################################
  ## extract data from fitted model
  ## original data pattern
  X <- response(object)
  ## spatial locations and weights used in fit
  Q <- quad.ppm(object2D)
  quadpoints <- union.quad(Q)
  Z <- is.data(Q)
  wts <- w.quad(Q)
  ## subset of quadrature points used to fit model
  subQset <- getglmsubset(object2D)
  if(is.null(subQset)) subQset <- rep.int(TRUE, n.quad(Q))
  
  #################################################################
  #' trap case where covariate = "<name>"
  if(is.character(covariate) && (length(covariate) == 1)) {
    #' covariate is a single string
    is.cartesian <- covariate %in% c("x", "y")
    if(!is.cartesian) {
      #' not a reserved name; convert to an expression and evaluate later
      covariate <- str2expression(covariate)
    }
  } else {
    is.cartesian <- FALSE
  }

  #################################################################
  ## compute the covariate
  covunits <- NULL
  if(is.cartesian) {
    #' covariate is name of cartesian coordinate
    switch(covariate,
           x = {
             covvalues <- quadpoints$x
             covrange <- Frame(quadpoints)$xrange
           },
           y = {
             covvalues <- quadpoints$y
             covrange <- Frame(quadpoints)$yrange
           })
    covunits <- unitname(quadpoints)
  } else if(is.im(covariate)) {
    covvalues <- covariate[quadpoints, drop=FALSE]
    covrange <- internal$covrange %orifnull% range(covariate, finite=TRUE)
  } else if(is.function(covariate)) {
    covvalues <- pointweights(quadpoints,
                              weights=covariate,
                              weightsname="covariate")
    covrange <- internal$covrange %orifnull% range(covvalues, finite=TRUE)
    if(inherits(covariate, c("distfun", "distfunlpp")))
      covunits <- unitname(quadpoints)
  } else if(is.vector(covariate) && is.numeric(covariate)) {
    covvalues <- covariate
    covrange <- internal$covrange %orifnull% range(covvalues, finite=TRUE)
    if(length(covvalues) != npoints(quadpoints))
      stop("Length of covariate vector,", length(covvalues), "!=",
           npoints(quadpoints), ", number of quadrature points")
  } else if(is.expression(covariate)) {
    ## Expression involving covariates in the fitted object
    if(!is.null(CovariatesDF)) {
      ## Expression may involve an external variable
      neednames <- all.vars(covariate)
      missingnames <- setdiff(neednames, colnames(CovariatesDF))
      if(length(missingnames)) {
        ## missing variables should be 'external'
        foundvars <- mget(missingnames,
                          parent.frame(),
                          ifnotfound=rep(list(NULL), length(missingnames)))
        bad <- sapply(foundvars, is.null)
        if(any(bad)) {
          nbad <- sum(bad)
          stop(paste(ngettext(nbad, "Variable", "Variables"),
                     commasep(sQuote(missingnames[bad])),
                     "not found"), call.=FALSE)
        }
        founddata <- mpl.get.covariates(foundvars, quadpoints)
        CovariatesDF <- cbind(CovariatesDF, founddata)
      }
    }
    ## Evaluate expression
    sp <- parent.frame()
    covvalues <- eval(covariate, envir=CovariatesDF, enclos=sp)
    covrange <- internal$covrange %orifnull% range(covvalues, finite=TRUE)
    if(!is.numeric(covvalues))
      stop("The evaluated covariate is not numeric")
  } else 
    stop(paste("The", sQuote("covariate"), "should be either",
               "a pixel image, an expression or a numeric vector"))

  #################################################################
  ## Secret exit
  if(identical(internal$getrange, TRUE))
    return(covrange)
    
  ################################################################
  ## Residuals/marks attached to appropriate locations.
  ## Stoyan-Grabarnik weights are attached to the data points only.
  ## Others (residuals) are attached to all quadrature points.
  resvalues <- 
    if(!is.null(rv)) rv
    else if(type=="eem") eem(object, check=check)
    else residuals(object, type=type, check=check)

  if(inherits(resvalues, "msr")) {
    ## signed or vector-valued measure; extract increment masses
    resvalues <- resvalues$val
    if(ncol(as.matrix(resvalues)) > 1)
      stop("Not implemented for vector measures; use [.msr to split into separate components")
  }

  if(inherits(resvalues, "imlist")) {
    if(length(resvalues) > 1) 
      stop("Not implemented for vector-valued residuals")
    resvalues <- resvalues[[1]]
  }

  if(is.im(resvalues))
     resvalues <- resvalues[quadpoints]

  ## NAMES OF THINGS
  ## name of the covariate
  if(is.null(covname)) 
    covname <- if(is.expression(covariate)) covariate else "covariate"
  ## type of residual/mark
  if(missing(typename)) 
    typename <- if(!is.null(rv)) "rv" else ""

  clip <- !is.null(clipwindow)

  ## CALCULATE
  stuff <- LurkEngine(object=object,
                      type=type, cumulative=cumulative, plot.sd=plot.sd,
                      quadpoints=quadpoints,
                      wts=wts,
                      Z=Z,
                      subQset=subQset,
                      covvalues=covvalues,
                      resvalues=resvalues,
                      clip=clip,
                      clipwindow=clipwindow,
                      cov.is.im=is.im(covariate),
                      covrange=covrange,
                      typename=typename,
                      covname=covname,
                      covunits=covunits,
                      cl=cl, clenv=clenv,
                      oldstyle=oldstyle, check=check, verbose=verbose,
                      nx=nx, splineargs=splineargs,
                      envelope=envelope, nsim=nsim, nrank=nrank, Xsim=Xsim,
                      internal=internal)
    
  ## ---------------  PLOT ----------------------------------
  if(plot.it && inherits(stuff, "lurk")) {
    plot(stuff, ...)
    return(invisible(stuff))
  } else {
    return(stuff)
  }
}
