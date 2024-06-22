#'
#'   evidencelppm.R
#'
#'   spatialCovariateEvidence method for class lppm
#'
#'   $Revision: 1.14 $ $Date: 2024/06/22 10:14:04 $


spatialCovariateEvidence.lppm <- local({

  spatialCovariateEvidence.lppm <- function(model, covariate, ...,
                             lambdatype=c("cif", "trend", "intensity"),
                             eps=NULL, dimyx=NULL, xy=NULL,
                             rule.eps=c("adjust.eps",
                                        "grow.frame", "shrink.frame"),
                             delta=NULL, nd=NULL,
                             interpolate=TRUE,
                             jitter=TRUE, jitterfactor=1, 
                             modelname=NULL, covname=NULL,
                             dataname=NULL, subset=NULL, clip.predict=TRUE) {
    lambdatype <- match.arg(lambdatype)
    rule.eps <- match.arg(rule.eps)
    #' evaluate covariate values at data points and at pixels
    ispois <- is.poisson(model)
    csr <- ispois && is.stationary(model)
    #' arguments controlling resolution
    pixels.given <- !(is.null(eps) && is.null(dimyx) && is.null(xy))
    sampling.given <- !(is.null(delta) && is.null(nd))
    resolution.given <- pixels.given || sampling.given

    #' determine names
    if(is.null(modelname))
      modelname <- if(csr) "CSR" else short.deparse(substitute(model))
    if(is.null(covname)) {
      covname <- singlestring(short.deparse(substitute(covariate)))
      if(is.character(covariate)) covname <- covariate
    }
    if(is.null(dataname))
      dataname <- model$Xname
    info <-  list(modelname=modelname, covname=covname,
                  dataname=dataname, csr=csr, ispois=ispois, 
                  spacename="linear network")

    #' convert character covariate to function
    if(is.character(covariate)) {
      #' One of the characters 'x' or 'y'
      #' Turn it into a function.
      ns <- length(covariate)
      if(ns == 0) stop("covariate is empty")
      if(ns > 1) stop("more than one covariate specified")
      covname <- covariate
      covariate <- switch(covariate,
                          x=xcoordfun,
                          y=ycoordfun,
                          stop(paste("Unrecognised covariate",
                                     dQuote(covariate))))
    }
  
    #' extract model components
    X <- model$X
    fit <- model$fit
    #'
    L <- as.linnet(X)
    Q <- quad.ppm(fit)
    #' restrict to subset if required
    if(!is.null(subset)) {
      X <- X[subset]
      Q <- Q[subset]
    }
    isdat <- is.data(Q)
    U <- union.quad(Q)
    wt <- w.quad(Q)
  
    #' evaluate covariate
    if(!is.marked(model)) {
      #' ...................  unmarked .......................
      if(is.im(covariate)) {
        if(is.linim(covariate)) {
          type <- "linim"
          Zimage <- covariate
          if(resolution.given) Zimage <- as.linim(Zimage,
                                                  eps=eps, dimyx=dimyx, xy=xy,
                                                  rule.eps=rule.eps,
                                                  delta=delta, nd=nd)
        } else {
          type <- "im"
          Zimage <- as.linim(covariate, L,
                             eps=eps, dimyx=dimyx, xy=xy,
                             rule.eps=rule.eps,
                             delta=delta, nd=nd)
        }
        if(!interpolate) {
          #' look up covariate values at quadrature points
          Zvalues <- safelookup(covariate, U)
        } else {
          #' evaluate at quadrature points by interpolation
          Zvalues <- interp.im(covariate, U$x, U$y)
          #' fix boundary glitches
          if(any(uhoh <- is.na(Zvalues)))
            Zvalues[uhoh] <- safelookup(covariate, U[uhoh])
        }
        #' extract data values
        ZX <- Zvalues[isdat]
      } else if(is.function(covariate)) {
        type <- "function"
        Zimage <- as.linim(covariate, L,
                           eps=eps, dimyx=dimyx, xy=xy,
                           rule.eps=rule.eps,
                           delta=delta, nd=nd)
        #' evaluate exactly at quadrature points
        Zvalues <- covariate(U$x, U$y)
        if(!all(is.finite(Zvalues)))
          warning("covariate function returned NA or Inf values")
        #' extract data values
        ZX <- Zvalues[isdat]
        #' collapse function body to single string
        covname <- singlestring(covname)
      } else if(is.null(covariate)) {
        stop("The covariate is NULL", call.=FALSE)
      } else stop(paste("The covariate should be",
                        "an image, a function(x,y)",
                        "or one of the characters",
                        sQuote("x"), "or", sQuote("y")),
                  call.=FALSE)
      #' corresponding fitted [conditional] intensity values
      lambda <- as.vector(predict(model, locations=U, type=lambdatype))
    } else {
      #' ...................  marked .......................
      if(!is.multitype(model))
        stop("Only implemented for multitype models (factor marks)")
      marx <- marks(U, dfok=FALSE)
      possmarks <- levels(marx)
      #' single image: replicate 
      if(is.im(covariate)) {
        covariate <- rep(list(covariate), length(possmarks))
        names(covariate) <- possmarks
      }
      #'
      if(is.list(covariate) && all(sapply(covariate, is.im))) {
        #' list of images
        if(length(covariate) != length(possmarks))
          stop("Number of images does not match number of possible marks")
        #' determine type of data
        islinim <- sapply(covariate, is.linim)
        type <- if(all(islinim)) "linim" else "im"
        Zimage <- as.solist(covariate)
        #' convert 2D pixel images to 'linim'
        Zimage[!islinim] <- lapply(Zimage[!islinim], as.linim, L=L,
                                   eps=eps, dimyx=dimyx, xy=xy,
                                   rule.eps=rule.eps,
                                   delta=delta, nd=nd)
        if(resolution.given) 
          Zimage[islinim] <- lapply(Zimage[!islinim], as.linim,
                                   eps=eps, dimyx=dimyx, xy=xy,
                                   rule.eps=rule.eps,
                                   delta=delta, nd=nd)
        #' evaluate covariate at each data point by interpolation
        Zvalues <- numeric(npoints(U))
        for(k in seq_along(possmarks)) {
          ii <- (marx == possmarks[k])
          covariate.k <- covariate[[k]]
          if(!interpolate) {
            #' direct lookup
            values <- safelookup(covariate.k, U[ii])
          } else {
            #' interpolation
            values <- interp.im(covariate.k, x=U$x[ii], y=U$y[ii])
            #' fix boundary glitches
            if(any(uhoh <- is.na(values)))
              values[uhoh] <- safelookup(covariate.k, U[ii][uhoh])
          }
          Zvalues[ii] <- values
        }
        #' extract data values
        ZX <- Zvalues[isdat]
        #' corresponding fitted [conditional] intensity values
        lambda <- predict(model, locations=U, type=lambdatype)
        if(length(lambda) != length(Zvalues))
          stop("Internal error: length(lambda) != length(Zvalues)")
      } else if(is.function(covariate)) {
        type <- "function"
        #' evaluate exactly at quadrature points
        Zvalues <- functioncaller(x=U$x, y=U$y, m=marx, f=covariate, ...)
        #' functioncaller: function(x,y,m,f,...) { f(x,y,m,...) }
        #' extract data values
        ZX <- Zvalues[isdat]
        #' corresponding fitted [conditional] intensity values
        lambda <- predict(model, locations=U, type=lambdatype)
        if(length(lambda) != length(Zvalues))
          stop("Internal error: length(lambda) != length(Zvalues)")
        #' images
        Zimage <- list()
        for(k in seq_along(possmarks))
          Zimage[[k]] <- as.linim(functioncaller, L=L, m=possmarks[k],
                                  f=covariate,
                                  eps=eps, dimyx=dimyx, xy=xy,
                                  rule.eps=rule.eps,
                                  delta=delta, nd=nd)
        #' collapse function body to single string
        covname <- singlestring(covname)
      } else if(is.null(covariate)) {
        stop("The covariate is NULL", call.=FALSE)
      } else stop(paste("For a multitype point process model,",
                        "the covariate should be an image, a list of images,",
                        "a function(x,y,m)", 
                        "or one of the characters",
                        sQuote("x"), "or", sQuote("y")),
                  call.=FALSE)
    }    
    #' ..........................................................

    #' apply jittering to avoid ties
    if(jitter) {
      ZX <- jitter(ZX, factor=jitterfactor)
      Zvalues <- jitter(Zvalues, factor=jitterfactor)
    }

    lambdaname <- if(is.poisson(model)) "intensity" else lambdatype
    lambdaname <- paste("the fitted", lambdaname)
    check.finite(lambda, xname=lambdaname, usergiven=FALSE)
    check.finite(Zvalues, xname="the covariate", usergiven=TRUE)

    #' lambda values at data points
    lambdaX <- predict(model, locations=X, type=lambdatype)

    #' lambda image(s)
    lambdaimage <- predict(model, type=lambdatype)
    
    #' restrict image to subset 
    if(clip.predict && !is.null(subset)) {
      Zimage      <- applySubset(Zimage, subset)
      lambdaimage <- applySubset(lambdaimage, subset)
    }

    #' wrap up 
    values <- list(Zimage      = Zimage,
                   lambdaimage = lambdaimage,
                   Zvalues     = Zvalues,
                   lambda      = lambda,
                   lambdaX     = lambdaX,
                   weights     = wt,
                   ZX          = ZX,
                   type        = type)
    return(list(values=values, info=info, X=X))
  }

  xcoordfun <- function(x,y,m){x}
  ycoordfun <- function(x,y,m){y}

  functioncaller <- function(x,y,m,f,...) {
    nf <- length(names(formals(f)))
    if(nf < 2) stop("Covariate function must have at least 2 arguments")
    value <- if(nf == 2) f(x,y) else if(nf == 3) f(x,y,m) else f(x,y,m,...)
    return(value)
  }

  applySubset <- function(X, subset) {
    if(is.im(X)) return(X[subset, drop=FALSE])
    if(is.imlist(X)) return(solapply(X, "[", i=subset, drop=FALSE))
    return(NULL)
  }
  spatialCovariateEvidence.lppm
})

spatialCovariateEvidence.exactlppm <- local({

  spatialCovariateEvidence.exactlppm <- function(model, covariate, ...,
                                 lambdatype=c("cif", "trend", "intensity"),
                                 interpolate=TRUE,
                                 jitter=TRUE, jitterfactor=1,
                                 modelname=NULL, covname=NULL,
                                 dataname=NULL, subset=NULL, clip.predict=TRUE) {
    lambdatype <- match.arg(lambdatype)
    dont.complain.about(lambdatype)
    if(is.null(covariate)) stop("The covariate is NULL")
    #' evaluate covariate values at data points and at pixels
    ispois <- TRUE
    csr <- is.null(model$baseline)
    #' determine names
    if(is.null(modelname))
      modelname <- if(csr) "CSR" else "model with baseline"
    if(is.null(covname)) {
      if(is.character(covariate)) covname <- covariate else
      covname <- singlestring(short.deparse(substitute(covariate)))
    }
    if(is.null(dataname))
      dataname <- "X"
    
    info <-  list(modelname=modelname, covname=covname,
                  dataname=dataname, csr=csr, ispois=ispois,
                  spacename="linear network")
  
    X <- model$X
    L <- domain(X)

    Lfull <- Zfull <- NULL
    
    if(!is.null(subset)) {
      #' restrict to subset
      if(!clip.predict) {
        ## use original window for prediction
        Lfull <- L
      }
      X <- X[subset]
      L <- L[subset, drop=FALSE]
    } 
    
    #' evaluate covariate 
    if(is.character(covariate)) {
      #' One of the characters 'x' or 'y'
      #' Turn it into a function.
      ns <- length(covariate)
      if(ns == 0) stop("covariate is empty")
      if(ns > 1) stop("more than one covariate specified")
      covname <- covariate
      covariate <- switch(covariate,
                          x=xcoordfun,
                          y=ycoordfun,
                          stop(paste("Unrecognised covariate",
                                     dQuote(covariate))))
    } 
  
    if(!is.marked(X)) {
      #' ...................  unmarked .......................
      if(is.im(covariate)) {
        type <- "im"
      } else if(is.function(covariate)) {
        type <- "function"
      } else stop(paste("The covariate should be",
                        "an image, a function(x,y)",
                        "or one of the characters",
                        sQuote("x"), "or", sQuote("y")),
                  call.=FALSE)
      #' evaluate covariate at data points
      ZX <- evaluateNetCovariate(covariate, X, ...)
      #' pixel image of covariate values on network 
      Z <- evaluateNetCovariate(covariate, L, ...)
      if(!is.null(Lfull))
        Zfull <- evaluateNetCovariate(covariate, Lfull, ...)
      #' extract covariate values at equally-spaced sample points
      df <- attr(Z, "df")
      Zvalues <- df$values
      #' corresponding fitted intensity values
      lambdaimage <- predict(model, locations=L)
      ldf <- attr(lambdaimage, "df")
      lambdavalues <- ldf$values
      if(length(Zvalues) != length(lambdavalues))
        stop("Internal error: length mismatch in spatial covariate evidence")
      #' average weight
      sampleweight <- volume(L)/nrow(df)
    } else {
      #' ...................  marked .......................
      if(!is.multitype(X))
        stop("Only implemented for multitype models (factor marks)")
      marx <- marks(X, dfok=FALSE)
      possmarks <- levels(marx)
      npts <- npoints(X)
      if(is.im(covariate)) {
        #' single image: replicate 
        covariate <- rep(list(covariate), times=length(possmarks))
        names(covariate) <- as.character(possmarks)
      }
      #'
      if(is.list(covariate) && all(sapply(covariate, is.im))) {
        #' list of images
        type <- "im"
        if(length(covariate) != length(possmarks))
          stop("Number of images does not match number of possible marks")
        #' evaluate covariate at each data point 
        ZX <- numeric(npts)
        for(k in seq_along(possmarks)) {
          ii <- (marx == possmarks[k])
          covariate.k <- covariate[[k]]
          values <- evaluateNetCovariate(covariate.k, X[ii], ...)
          ZX[ii] <- values
        }
        #' evaluate covariate images throughout network
        Z <- solapply(covariate, "[", i=L, drop=FALSE)
        #' and evaluate in full domain, for prediction
        if(!is.null(Lfull))
          Zfull <- solapply(covariate, "[", i=Lfull, drop=FALSE)
      } else if(is.function(covariate)) {
        type <- "function"
        #' collapse function body to single string
        covname <- singlestring(covname)
        #' evaluate exactly at data points
        coX <- coords(X)
        ZX <- functioncaller(x=coX$x, y=coX$y, m=marx, f=covariate, ...)
        #' functioncaller: function(x,y,m,f,...) { f(x,y,m,...) }
        #' covariate in window
        Z <- list()
        for(k in seq_along(possmarks))
          Z[[k]] <- as.linim(functioncaller, L=L,
                             m=possmarks[k], f=covariate, ...)
        if(!is.null(Lfull)) {
          #' covariate in original window, for prediction
          Zfull <- list()
          for(k in seq_along(possmarks))
            Zfull[[k]] <- as.linim(functioncaller, L=Lfull,
                                m=possmarks[k], f=covariate, ...)
        }
      } else stop(paste("For a multitype point process model,",
                        "the covariate should be an image, a list of images,",
                        "a function(x,y,m)", 
                        "or one of the characters",
                        sQuote("x"), "or", sQuote("y")),
                  call.=FALSE)

      #' extract covariate values at all locations on network
      Zframes <- lapply(Z, attr, which="df")
      Zvalues <- unlist(lapply(Zframes, getElement, name="value"))
      #' compute fitted intensity
      lambdaimages <- predict(model, locations=L, type=lambdatype)
      #' extract intensity values at all locations on network
      lambdaframes <- lapply(lambdaimages, attr, which="df")
      lambdavalues <- unlist(lapply(lambdaframes, getElement, name="value"))
      if(length(lambdavalues) != length(Zvalues))
        stop("Internal error: length mismatch in spatial covariate evidence")
      #' average weight per sample point
      sampleweight <- volume(L)/nrow(df)        
    }    
    #' ..........................................................
    #' apply jittering to avoid ties
    if(jitter) {
      ZX <- jitter(ZX, factor=jitterfactor)
      Zvalues <- jitter(Zvalues, factor=jitterfactor)
    }

    check.finite(lambdavalues, xname="the fitted intensity", usergiven=FALSE)
    check.finite(Zvalues, xname="the covariate", usergiven=TRUE)

    #' lambda values at data points
    lambdaX <- predict(model, locations=X)

    #' lambda image(s)
    lambdaimage <- predict(model, locations=Lfull %orifnull% L)
    
    #' wrap up 
    values <- list(Zimage      = Zfull %orifnull% Z,
                   lambdaimage = lambdaimage,
                   Zvalues     = Zvalues,
                   lambda      = lambdavalues,
                   lambdaX     = lambdaX,
                   weights     = sampleweight,
                   ZX          = ZX,
                   type        = type)
    return(list(values=values, info=info, X=X)) # X is possibly a subset of original
  }

  xcoordfun <- function(x,y,m){x}
  ycoordfun <- function(x,y,m){y}

  ## Function caller used for marked locations (x,y,m) only.
  functioncaller <- function(x,y,m,f,...) {
    nf <- length(names(formals(f)))
    if(nf < 2) stop("Covariate function must have at least 2 arguments")
    if(nf == 2) return(f(x,y))
    if(nf == 3) return(f(x,y,m))
    argh <- list(...)
    extra <- intersect(names(argh),
                       names(formals(f))[-(1:3)])
    value <- do.call(f, append(list(x,y,m), argh[extra]))
    return(value)
  }
            
  spatialCovariateEvidence.exactlppm
})

