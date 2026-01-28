#
#  lppm.R
#
#  Point process models on a linear network
#
#  $Revision: 1.70 $   $Date: 2026/01/21 06:26:39 $
#

lppm <- function(X, ...) {
  UseMethod("lppm")
}


lppm.formula <- function(X, interaction=NULL, ..., data=NULL) {
  ## remember call
  callstring <- paste(short.deparse(sys.call()), collapse = "")
  cl <- match.call()

  ########### INTERPRET FORMULA ##############################
  
  if(!inherits(X, "formula"))
    stop(paste("Argument 'X' should be a formula"),
         call.=FALSE)
  formula <- X
  
  if(spatstat.options("expand.polynom"))
    formula <- expand.polynom(formula)

  ## check formula has LHS and RHS. Extract them
  if(length(formula) < 3)
    stop(paste("Formula must have a left hand side"),
         call.=FALSE)
  Yexpr <- formula[[2L]]
  trend <- formula[c(1L,3L)]
  
  ## FIT #######################################
  thecall <- call("lppm", X=Yexpr, trend=trend,
                  data=data, interaction=interaction)
  ncall <- length(thecall)
  argh <- list(...)
  nargh <- length(argh)
  if(nargh > 0) {
    thecall[ncall + 1:nargh] <- argh
    names(thecall)[ncall + 1:nargh] <- names(argh)
  }
  callenv <- list2env(as.list(data), parent=parent.frame())
  result <- eval(thecall, envir=callenv)

  result$Xname <- short.deparse(Yexpr)
  
  result$call <- cl
  result$callstring <- callstring

  return(result)
}

lppm.lpp <- function(X, ..., eps=NULL, nd=1000, random=FALSE) {
  Xname <- short.deparse(substitute(X))
  callstring <- paste(short.deparse(sys.call()), collapse = "")
  cl <- match.call()
  nama <- names(list(...))
  resv <- c("method", "forcefit")
  if(any(clash <- resv %in% nama))
    warning(paste(ngettext(sum(clash), "Argument", "Arguments"),
                  commasep(sQuote(resv[clash])),
                  "must not be used"),
            call.=FALSE)
  stopifnot(inherits(X, "lpp"))
  quadarg <- substitute(linequad(X, eps=eps, nd=nd, random=random),
                        list(X=substitute(X),
                             eps=substitute(eps),
                             nd=substitute(nd),
                             random=substitute(random)))
  fit <- do.call(ppm,
                 list(quadarg, ..., method="mpl", forcefit=TRUE),
                 envir=parent.frame())
  if(!is.poisson.ppm(fit))
    warning("Non-Poisson models currently use Euclidean distance",
            call.=FALSE)
  out <- list(X=X, fit=fit, Xname=Xname, call=cl, callstring=callstring)
  class(out) <- "lppm"
  return(out)
}

is.lppm <- function(x) { inherits(x, "lppm") }

# undocumented
as.ppm.lppm <- function(object) { object$fit }

##

getglmdata.lppm <- function(object, ...) { getglmdata(as.ppm(object), ...) }

getglmfit.lppm <- function(object, ...) { getglmfit(as.ppm(object), ...) }

getglmsubset.lppm <- function(object, ...) { getglmsubset(as.ppm(object), ...) }

hasglmfit.lppm <- function(object) { hasglmfit(as.ppm(object)) }

##

fitted.lppm <- function(object, ..., dataonly=FALSE, new.coef=NULL,
                        leaveoneout=FALSE) {
  pfit <- object$fit
  v <- fitted(pfit, dataonly=dataonly, new.coef=new.coef,
              leaveoneout=leaveoneout)
  return(v)
}
  
predict.lppm <- local({
  
  predict.lppm <- function(object, ..., 
                           type="trend", locations=NULL, covariates=NULL,
                           se=FALSE,
                           new.coef=NULL) {

    type <- pickoption("type", type,
                       c(trend="trend", cif="cif", lambda="cif",
                         intensity="intensity"))

    if(type == "intensity" && !is.poisson(object)) 
      stop(paste("Intensity calculation is not yet implemented",
                 "for non-Poisson models on a network"),
           call.=FALSE)
      
    X <- object$X
    fit <- object$fit
    L <- as.linnet(X)

    ## functions to evaluate the local covariates
    LocalCoords <- list(seg = linfun(function(x,y,seg,tp) { seg }, L),
                        tp  = linfun(function(x,y,seg,tp) { seg }, L))

    if(!is.null(locations)) {
      ## determine whether 'locations' includes local coordinates
      if(is.data.frame(locations)) {
        ## data frame of spatial locations
        gotlocal <- all(c("seg", "tp") %in% names(locations))
      } else if(is.lpp(locations)) {
        ## point pattern on network
        gotlocal <- TRUE
        loci <- locations
        locations <- as.ppp(locations)
        attr(locations, "situ") <- loci
      } else {
        ## other spatial data
        gotlocal <- FALSE
      }
      values <- predict(fit, locations=locations, covariates=covariates,
                        type=type, se=se, new.coef=new.coef,
                        extracovariates=if(!gotlocal) LocalCoords else NULL)
      return(values)
    }
  
    #' locations not given; want a pixel image
    #' pixellate the lines
    Llines <- as.psp(L)
    linemask <- psp2mask(Llines, ...)
    lineimage <- as.im(linemask)
    #' extract pixel centres
    xx <- rasterx.mask(linemask)
    yy <- rastery.mask(linemask)
    mm <- linemask$m
    xx <- as.vector(xx[mm])
    yy <- as.vector(yy[mm])
    pixelcentres <- ppp(xx, yy, window=as.rectangle(linemask), check=FALSE)
    pixdf <- data.frame(xc=xx, yc=yy)
    #' project pixel centres onto lines
    p2s <- project2segment(pixelcentres, Llines)
    projloc <- as.data.frame(p2s$Xproj)
    projmap <- as.data.frame(p2s[c("mapXY", "tp")])
    projdata <- cbind(pixdf, projloc, projmap)
    #' predict at the projected points
    if(!is.multitype(fit)) {
      #' unmarked
      values <- predict(fit, locations=projloc, covariates=covariates,
                        type=type, se=se, new.coef=new.coef,
                        extracovariates=LocalCoords)
      if(!se) {
        out <- putvalues(values, lineimage, pixelcentres, projdata, L)
      } else {
        outfit <- putvalues(values[[1L]], lineimage, pixelcentres, projdata, L)
        outse  <- putvalues(values[[2L]], lineimage, pixelcentres, projdata, L)
        out <- solist(outfit, outse)
        names(out) <- c(type, "se")
      }
    } else {
      #' multitype 
      #' predict for each type
      lev <- levels(marks(data.ppm(fit)))
      out <- outfit <- outse <- list()
      for(k in seq(length(lev))) {
        markk <- factor(lev[k], levels=lev)
        locnk <- cbind(projloc, data.frame(marks=markk))
        values <- predict(fit, locations=locnk, covariates=covariates,
                          type=type, se=se, new.coef=new.coef,
                          extracovariates=LocalCoords)
        if(!se) {
          out[[k]] <- putvalues(values, lineimage, pixelcentres, projdata, L)
        } else {
          outfit[[k]] <-
            putvalues(values[[1L]], lineimage, pixelcentres, projdata, L)
          outse[[k]]  <-
            putvalues(values[[2L]], lineimage, pixelcentres, projdata, L)
        }
      }
      if(!se) {
        out <- as.solist(out)
        names(out) <- as.character(lev)
      } else {
        outfit <- as.solist(outfit)
        outse <- as.solist(outse)
        names(outfit) <- names(outse) <- as.character(lev)
        out <- list(outfit, outse)
        names(out) <- c(type, "se")
      }
    }
    return(out)
  }

  putvalues <- function(values, lineimage, pixelcentres, projdata, L) {
    ## insert numerical 'values' in the right places in a linim
    ## (1) image component: map to nearest pixels
    Z <- lineimage
    Z[pixelcentres] <- values
    ## (2) data frame: attach exact line position data
    df <- cbind(projdata, values)
    out <- linim(L, Z, df=df, restrict=FALSE)
    return(out)
  }

  predict.lppm
})


coef.lppm <- function(object, ...) {
  coef(object$fit)
}

print.lppm <- function(x, ...) {
  terselevel <- spatstat.options('terse')
  splat("Point process model on linear network")
  if(waxlyrical('extras', terselevel)) {
    splat("\tFitted to point pattern dataset", sQuote(x$Xname))
    parbreak(terselevel)
  }
  print(x$fit, showname=FALSE)
  parbreak(terselevel)
  if(waxlyrical('gory', terselevel)) {
    cat("Domain: ")
    print(as.linnet(x))
  }
  return(invisible(NULL))
}

summary.lppm <- function(object, ...) {
  splat("Point process model on linear network")
  splat("Fitted to point pattern dataset", sQuote(object$Xname))
  cat("Internal fit: ")
  print(summary(object$fit))
  terselevel <- spatstat.options('terse')
  if(waxlyrical('extras', terselevel)) {
    parbreak(terselevel)
    splat("Original data:", object$Xname)
    if(waxlyrical('gory', terselevel)) {
      cat("Domain: ")
      print(summary(as.linnet(object)))
    }
  }
  return(invisible(NULL))
}

plot.lppm <- function(x, ..., type="trend") {
  xname <- short.deparse(substitute(x))
  y <- predict(x, type=type)
  dont.complain.about(y)
  do.call(plot, resolve.defaults(list(quote(y)),
                                   list(...),
                                   list(main=xname)))
}
  
anova.lppm <- function(object, ..., test=NULL) {
  stuff <- list(object=object, ...)
  if(!is.na(hit <- match("override", names(stuff)))) {
    warning("Argument 'override' is outdated and was ignored", call.=FALSE)
    stuff <- stuff[-hit]
  }
  #' extract ppm objects where appropriate
  mod <- sapply(stuff, is.lppm)
  stuff[mod] <- lapply(stuff[mod], getElement, name="fit")
  #' analysis of deviance or adjusted composite deviance
  do.call(anova.ppm, append(stuff, list(test=test)))
}

update.lppm <- function(object, ...) {
  stopifnot(inherits(object, "lppm"))
  X <- object$X
  fit <- object$fit
  Xname <- object$Xname
  callframe <- environment(formula(fit))
  aargh <- list(...)
  islpp <- sapply(aargh, is.lpp)
  if(any(islpp)) {
    # trap point pattern argument & convert to quadscheme
    ii <- which(islpp)
    if((npp <- length(ii)) > 1)
      stop(paste("Arguments not understood:", npp, "lpp objects given"),
           call.=FALSE)
    X <- aargh[[ii]]
    aargh[[ii]] <- linequad(X)
    Xexpr <- sys.call()[[ii+2L]] %orifnull% expression(X)
    Xname <- short.deparse(Xexpr)
  }
  isfmla <- sapply(aargh, inherits, what="formula")
  if(any(isfmla)) {
    # trap formula pattern argument, update it, evaluate LHS if required
    jj <- which(isfmla)
    if((nf <- length(jj)) > 1)
      stop(paste("Arguments not understood:", nf, "formulae given"),
           call.=FALSE)
    fmla <- aargh[[jj]]
    fmla <- update(formula(object), fmla)
    if(!is.null(lhs <- lhs.of.formula(fmla))) {
      X <- eval(lhs, envir=list2env(list("."=X), parent=callframe))
      Qpos <- if(any(islpp)) ii else (length(aargh) + 1L)
      aargh[[Qpos]] <- linequad(X)
      Xname <- short.deparse(lhs)
    }
    aargh[[jj]] <- rhs.of.formula(fmla)
  }
  newfit <- do.call(update.ppm,
                    append(list(fit), aargh),
                    envir=callframe)
  if(!is.poisson.ppm(newfit))
    warning("Non-Poisson models currently use Euclidean distance",
            call.=FALSE)
  out <- list(X=X, fit=newfit, Xname=Xname)
  class(out) <- "lppm"
  return(out)
}

updateData.lppm <- function(model, X, ...) {
  update(model, X)
}

terms.lppm <- function(x, ...) {
  terms(x$fit, ...)
}

logLik.lppm <- function(object, ...) {
  logLik(object$fit, ...)
}

deviance.lppm <- function(object, ...) {
  deviance(object$fit, ...)
}

pseudoR2.lppm <- function(object, ..., keepoffset=TRUE) {
  dres <- deviance(object, ..., warn=FALSE)
  nullfmla <- . ~ 1
  if(keepoffset && has.offset.term(object)) {
    off <- attr(model.depends(object), "offset")
    offterms <- row.names(off)[apply(off, 1, any)]
    if(length(offterms)) {
      nullrhs <- paste(offterms, collapse=" + ") 
      nullfmla <- as.formula(paste(". ~ ", nullrhs))
    }
  } 
  nullmod <- update(object, nullfmla, forcefit=TRUE)
  dnul <- deviance(nullmod, warn=FALSE)
  return(1 - dres/dnul)
}

formula.lppm <- function(x, ...) {
  formula(x$fit, ...)
}

extractAIC.lppm <- function(fit, ...) {
  extractAIC(fit$fit, ...)
}

as.owin.lppm <- function(W, ..., fatal=TRUE) {
  stopifnot(inherits(W, "lppm"))
  as.owin(as.linnet(W), ..., fatal=fatal)
}

Window.lppm <- function(X, ...) { as.owin(X) }

data.lppm <- function(object) { object$X }

model.images.lppm <- local({

  model.images.lppm <- function(object, L=as.linnet(object), ...) {
    stopifnot(inherits(object, "lppm"))
    stopifnot(inherits(L, "linnet"))
    m <- model.images(object$fit, W=as.rectangle(L), ...)
    if(length(m)) {
      ## restrict images to L
      type <- if(is.hyperframe(m)) "hyperframe" else
              if(is.imlist(m)) "imlist" else
              if(is.list(m) && all(sapply(m, is.im))) "imlist" else
              stop("Internal error: model.images not understood", call.=FALSE)
      switch(type,
             imlist = {
               ## list of images: convert to list of linims
               ZL <- netmask(L, template=m[[1L]])
               m <- tolinims(m, L=L, imL=ZL)
             },
             hyperframe = {
               ## hyperframe, each column being a list of images
               ## extract columns
               rownam <- row.names(m)
               m <- as.list(m)
               ZL <- netmask(L, template=m[[1L]][[1L]])
               mm <- lapply(m, tolinims, L=L, imL=ZL)
               m <- do.call(hyperframe, mm)
               row.names(m) <- rownam
             })
    }
    return(m)
  }

  netmask <- function(L, template) {
    as.im(psp2mask(as.psp(L), xy=as.mask(template)))
  }
    
  tolinim <- function(x, L, imL) linim(L, eval.im(x * imL), restrict=FALSE)
  tolinims <- function(x, L, imL) solapply(x, tolinim, L=L, imL=imL)

  model.images.lppm
})

  
model.matrix.lppm <- function(object,
                              data=model.frame(object, na.action=NULL),
                             ..., keepNA=TRUE) {
  stopifnot(is.lppm(object))
  if(missing(data)) data <- NULL
  model.matrix(object$fit, data=data, ..., keepNA=keepNA)
}

model.frame.lppm <- function(formula, ...) {
  stopifnot(inherits(formula, "lppm"))
  model.frame(formula$fit, ...)
}

domain.lppm <- as.linnet.lppm <- function(X, ...) {
  if(is.NAobject(X)) NAobject("linnet") else as.linnet(X$X, ...)
}

nobs.lppm <- function(object, ...) {
  if(is.NAobject(object)) NA_integer_ else npoints(object$X)
}

is.poisson.lppm <- function(x) { is.poisson(x$fit) }

is.stationary.lppm <- function(x) { is.stationary(x$fit) }

is.multitype.lppm <- function(X, ...) { is.multitype(X$fit) }

is.marked.lppm <- function(X, ...) { is.marked(X$fit) }

vcov.lppm <- function(object, ...) {
  if(!is.poisson(object))
    stop("vcov.lppm is only implemented for Poisson models", call.=FALSE)
  vcov(object$fit, ...)
}

valid.lppm <- function(object, ...) {
  valid(object$fit, ...)
}

emend.lppm <- function(object, ...) {
  object$fit <- emend(object$fit, ...)
  return(object)
}

response.lppm <- function(object) { data.lppm(object) }

parres.lppm <- function(model, covariate, ...,
                        smooth.effect=FALSE, subregion=NULL,
                        bw="nrd0", adjust=1, from=NULL,to=NULL, n=512,
                        bw.input = c("points", "quad"),
                        bw.restrict = FALSE,
                        covname) {  
  callstring <- paste(deparse(sys.call()), collapse = "")
  modelname <- short.deparse(substitute(model))

  stopifnot(is.lppm(model))

  if(is.marked(model))
    stop("Sorry, this is not yet implemented for marked models")
  
  if(missing(covariate)) {
    mc <- model.covariates(model)
    if(length(mc) == 1) covariate <- mc else stop("covariate must be provided")
  }
  if(missing(covname)) 
    covname <- sensiblevarname(deparse(substitute(covariate)), "X")

  if(is.null(adjust)) adjust <- 1

  bw.input <- match.arg(bw.input)
  
  partialResidualEngine(model         = as.ppm(model),
                        covariate     = covariate,
                        ...,
                        smooth.effect = smooth.effect,
                        subregion     = subregion,
                        bw            = bw,
                        adjust        = adjust,
                        from          = from,
                        to            = to,
                        n             = n,
                        bw.input      = bw.input,
                        bw.restrict   = bw.restrict,
                        covname       = covname,
                        callstring    = callstring,
                        modelname     = modelname,
                        do.variance   = is.poisson(model))
}

eem.lppm <- function(fit, ...) {
  verifyclass(fit, "lppm")
  lambda <- fitted(fit, dataonly=TRUE, ...)
  eemarks <- 1/lambda
  attr(eemarks, "type") <- "eem"
  attr(eemarks, "typename") <- "exponential energy marks"
  return(eemarks)
}

residuals.lppm <- function(object, type="raw", ...) {
  res <- residuals(as.ppm(object), type=type, ...)
  attr(res, "situ") <- attr(quad.ppm(object), "situ")
  return(res)
}
