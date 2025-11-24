##
##                            diaglppm.R
##
## Makes diagnostic plots based on residuals or energy weights
## for a point process model on a network
##
## $Revision: 1.3 $ $Date: 2025/11/24 02:48:11 $
##

diagnose.lppm <- function(object, ..., type="raw", which="all", 
                         sigma=NULL, 
                         cumulative=TRUE,
                         plot.it = TRUE, rv = NULL, 
                         compute.sd=is.poisson(object), compute.cts=TRUE,
                         envelope=FALSE, nsim=39, nrank=1,
                         typename, oldstyle=FALSE, 
                         splineargs=list(spar=0.5))
{
  asked.newstyle <- !missing(oldstyle) && !oldstyle

  stopifnot(is.lppm(object))
  if(is.marked(object))
    stop("Sorry, this is not yet implemented for marked models")
  
  ## check whether model originally came from replicated data
  is.subfit <- (object$fit$method == "mppm")

  ## -------------  Interpret arguments --------------------------

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

  ## 'which' is multiple choice with exact matching 
  optionlist <- c("all", "marks", "smooth", "x", "y", "sum")

  if(!all(m <- which %in% optionlist))
    stop(paste("Unrecognised choice(s) of",
               paste(sQuote("which"), ":", sep=""),
               paste(which[!m], collapse=", ")))

  opt <- list()
  opt$all <- "all" %in% which
  opt$marks <-  ("marks" %in% which)   | opt$all
  opt$smooth <- ("smooth" %in% which)  | opt$all
  opt$xmargin <- (("x" %in% which)       | opt$all) && !cumulative
  opt$ymargin <- (("y" %in% which)       | opt$all) && !cumulative
  opt$xcumul <-  (("x" %in% which)       | opt$all) && cumulative
  opt$ycumul <-  (("y" %in% which)       | opt$all) && cumulative
  opt$sum <-     ("sum" %in% which)      | opt$all

  ## compute and plot estimated standard deviations?
  ## yes for Poisson, no for other models, unless overridden
  if(!missing(compute.sd))
    plot.sd <- compute.sd
  else
    plot.sd <- list(...)$plot.sd
  if(is.null(plot.sd))
    plot.sd <- is.poisson(object)
  if(missing(compute.sd))
    compute.sd <- plot.sd

  ## default for mppm objects is oldstyle=TRUE
  if(compute.sd && is.subfit) {
    if(!asked.newstyle) {
      ## silently change default
      oldstyle <- TRUE
    } else {
      stop(paste("Variance calculation for a subfit of an mppm object",
                 "is only implemented for oldstyle=TRUE"),
           call.=FALSE)
    }
  }
    
  ## interpolate the density of the residual measure?
  if(missing(compute.cts)) {
    plot.neg <- resolve.defaults(list(...),
                                 formals(plot.diagppm)["plot.neg"])$plot.neg
    ## only if it is needed for the mark plot
    compute.cts <- opt$marks && (plot.neg != "discrete")
  }

  ## -------  DO THE CALCULATIONS -----------------------------------
  RES <-  diagLppmEngine(object, type=type, typename=typename,
                         opt=opt, sigma=sigma, 
                         compute.sd=compute.sd,
                         compute.cts=compute.cts,
                         envelope=envelope, nsim=nsim, nrank=nrank,
                         rv=rv, oldstyle=oldstyle,
                         splineargs=splineargs,
                         ...)

  RES$typename <- typename
  RES$opt <- opt
  RES$compute.sd <- compute.sd
  RES$compute.cts <- compute.cts
  
  class(RES) <- "diaglppm"

  ## -------  PLOT --------------------------------------------------
  if(plot.it) 
    plot(RES, ...)

  return(RES)
}

## >>>>>>>>>>>>>>>>>>> main calculation code <<<<<<<<<<<<<<<<<<<<<<<<<

diagLppmEngine <- function(object, ..., type="eem", typename, opt,
                           sigma=NULL,
                           compute.sd=is.poisson(object),
                           compute.cts=TRUE,
                           envelope=FALSE, nsim=39, nrank=1,
                           rv=NULL, oldstyle=FALSE,
                           splineargs = list(spar=0.5),
                           verbose=TRUE)
{
  if(is.marked(object))
    stop("Sorry, this is not yet implemented for marked models")

  ## original point pattern data on network
  X <- data.lppm(object)
  L <- as.linnet(X)

  ## quadrature points (2D) and weights
  Q <- quad.ppm(object)
  U2D <- union.quad(Q)
  Qweights <- w.quad(Q)
  W <- Window(Q)

  ## quadrature points on the network
  U <- attr(Q, "plekken") %orifnull% lpp(U2D, L)

  ## 
  ppmfit <- as.ppm(object)

  ## -------------- Calculate residuals/weights -------------------

  ## Discretised residuals

  if(type == "eem") {
    residval <- if(!is.null(rv)) rv else eem(object, check=FALSE)
    residval <- as.numeric(residval)
    Y <- X %mark% residval
  } else {
    if(!is.null(rv) && !inherits(rv, "msr"))
      stop("rv should be a measure (object of class msr)")
    residobj <-
      if(!is.null(rv)) rv else residuals(object, type=type, check=FALSE)
    residval <- with(residobj, "increment")
    if(ncol(as.matrix(residval)) > 1L)
      stop("Not implemented for vector-valued residuals; use [.msr to split into separate components")
    Y <- U %mark% residval
  }

  ## Atoms and density of measure

  Ymass <- NULL
  Ycts  <- NULL
  Ydens <- NULL

  if(compute.cts) {
    if(type == "eem") {
      Ymass <- Y
      Ycts  <- U %mark% (-1)
      Ydens <- as.linim(-1, L)
    } else {
      atoms <- with(residobj, "is.atom")
      masses <- with(residobj, "discrete")
      cts    <- with(residobj, "density")
      if(!is.null(atoms) && !is.null(masses) && !is.null(cts)) {
        Ymass <- (U %mark% masses)[atoms]
        Ycts    <- U %mark% cts
        ## remove NAs (as opposed to zero cif points)
        if(!all(ok <- is.finite(cts))) {
          U <- U[ok]
          Ycts <- Ycts[ok]
          cts  <- cts[ok]
          Qweights <- Qweights[ok]
        }
        ## interpolate continuous part to yield an image for plotting
        if(type == "inverse" && all(cts != 0)) {
          Ydens <- as.linim(-1, L)
        } else if(is.stationary(object) && is.poisson(object)) {
          ## all values of `cts' will be equal
          Ydens <- as.linim(cts[1L], L)
        } else {
          smallsigma <- max(nndist(Ycts))
          Ujitter <- rjitter(U, smallsigma)
          Ydens <- Smooth(Ujitter %mark% marks(Ycts),
                          sigma=smallsigma,
                          weights=Qweights,
                          edge=TRUE, ...)
        }
      }
    }
  }
    
  ## ------------ start collecting results -------------------------
  
  result <- list(type=type,
                 Y=Y,
                 W=W,
                 Ymass=Ymass,
                 Ycts=Ycts,
                 Ydens=Ydens)

  ## ------------- smoothed field ------------------------------

  Z <- NULL
  if(opt$smooth || opt$xcumul || opt$ycumul || opt$xmargin || opt$ymargin) {
    if(is.null(sigma))
      sigma <- 0.07 * diameter(W)
    Ycm <- marks(Y)
    ok <- is.finite(Ycm)
    if(all(ok)) {
      Z <- density(Y, sigma, weights=Ycm, edge=TRUE, ...)
    } else {
      Z <- density(Y[ok], sigma, weights=Ycm[ok], edge=TRUE, ...)
    }
  }
  if(opt$smooth) {
    result$smooth <- list(Z = Z, sigma=sigma)
    if(type == "pearson")
      result$smooth$sdp <- 1/(2 * sigma * sqrt(pi))
  }

  ## -------------- marginals of smoothed field ------------------------

  if(opt$xmargin) {
    ma <- marginalIntegralOfLinim(Z, "x")
    xx <- ma$x
    Fx <- ma$integral
    if(type == "eem") {
      ## length of network to the left of x
      One <- eval.linim(1 + 0 * Z)
      Ema <- marginalIntegralOfLinim(One, "x")
      EFx <- Ema$integral
    } else {
      ## zero
      EFx <- numeric(length(xx))
    }
    result$xmargin <- list(x=xx, xZ=Fx, ExZ=EFx)
  }
  
  if(opt$ymargin) {
    ma <- marginalIntegralOfLinim(Z, "y")
    yy <- ma$y
    Fy <- ma$integral
    if(type == "eem") {
      ## length of network below height y
      One <- eval.linim(1 + 0 * Z)
      Ema <- marginalIntegralOfLinim(One, "y")
      EFy <- Ema$integral
    } else {
      ## zero
      EFy <- numeric(length(yy))
    }
    result$ymargin <- list(y=yy, yZ=Fy, EyZ=EFy)
  }
  

  
  ## -------------- cumulative (lurking variable) plots --------------

  ## precompute simulated patterns for envelopes
  if(identical(envelope, TRUE))
    envelope <- simulate(object, nsim=nsim, progress=verbose)

  if(opt$xcumul)
    result$xcumul <- 
    lurking(object, covariate=expression(x),
            type=type,
            rv=residval,
            plot.sd=compute.sd,
            envelope=envelope, nsim=nsim, nrank=nrank,
            plot.it=FALSE,
            typename=typename,
            covname="x coordinate",
            oldstyle=oldstyle,
            check=FALSE,
            splineargs=splineargs,
            ...)

  if(opt$ycumul)
    result$ycumul <- 
    lurking(object, covariate=expression(y),
            type=type,
            rv=residval,
            plot.sd=compute.sd,
            envelope=envelope, nsim=nsim, nrank=nrank,
            plot.it=FALSE,
            typename=typename,
            covname="y coordinate",
            oldstyle=oldstyle,
            check=FALSE,
            splineargs=splineargs,
            ...)

  ## -------------- summary numbers --------------
  
  if(opt$sum) 
    result$sum <- list(marksum=sum(marks(Y), na.rm=TRUE),
                       lengthL=volume(L),
                       sumquad=sum(Qweights),
                       range=if(!is.null(Z)) range(Z) else NULL)

  return(invisible(result))
}


################# support for class 'diagppm' ##########################

plot.diaglppm <-
  function(x, ..., which,
           plot.neg=c("image", "discrete"),
           plot.smooth=c("image", "persp"),
           plot.sd, spacing=0.1, outer=3, 
           srange=NULL, monochrome=FALSE, main=NULL)
{
  opt <- x$opt
  
  plot.neg <- match.arg(plot.neg)
  plot.smooth <- match.arg(plot.smooth)
  
  if(!missing(which)) {
    witches <- c("all", "marks", "smooth", "x", "y", "sum")
    unknown <- is.na(match(which, witches))
    if(any(unknown))
      warning(paste("Unrecognised",
                    ngettext(sum(unknown), "option", "options"),
                    "which =",
                    commasep(sQuote(which[unknown])),
                    ": valid options are",
                    commasep(sQuote(witches))), call.=FALSE)
    oldopt <- opt
    newopt <- list()
    newopt$all <- "all" %in% which
    newopt$marks <-  ("marks" %in% which)   | newopt$all
    newopt$smooth <- ("smooth" %in% which)  | newopt$all
    newopt$xmargin <- (("x" %in% which)       | newopt$all) && oldopt$xmargin
    newopt$ymargin <- (("y" %in% which)       | newopt$all) && oldopt$ymargin
    newopt$xcumul <-  (("x" %in% which)       | newopt$all) && oldopt$xcumul
    newopt$ycumul <-  (("y" %in% which)       | newopt$all)  && oldopt$ycumul
    newopt$sum <-     ("sum" %in% which)      | newopt$all

    illegal <- (unlist(newopt) > unlist(oldopt))
    if(any(illegal)) {
      offending <- paste(names(newopt)[illegal], collapse=", ")
      whinge <- paste("cannot display the following components;\n",
                      "they were not computed: - \n", offending, "\n")
      stop(whinge)
    }

    opt <- newopt
  }

  if(missing(plot.sd)) {
    plot.sd <- x$compute.sd
  } else if(plot.sd && !(x$compute.sd)) {
    warning("can't plot standard deviations; they were not computed")
    plot.sd <- FALSE
  }

  if(!(x$compute.cts) && (plot.neg != "discrete") && (opt$marks || opt$all)) {
    if(!missing(plot.neg))
      warning("can't plot continuous component of residuals; it was not computed")
    plot.neg <- "discrete"
  }
  
  if(opt$all) 
    resid4plot(RES=x,
               plot.neg=plot.neg, plot.smooth=plot.smooth,
               spacing=spacing, outer=outer,
               srange=srange, monochrome=monochrome, main=main, ...)
  else
    resid1plot(RES=x, opt=opt,
               plot.neg=plot.neg, plot.smooth=plot.smooth,
               srange=srange, monochrome=monochrome, main=main, ...)
}


print.diaglppm <- function(x, ...) {
  
  opt <- x$opt
  typename <- x$typename
  
  splat("Model diagnostics", paren(typename))
  splat("Diagnostics available:")
  optkey <- list(all="four-panel plot",
                 marks=paste("mark plot", if(!x$compute.cts)
                   "(discrete representation only)" else NULL),
                 smooth="smoothed residual field",
                 xmargin="x marginal density",
                 ymargin="y marginal density",
                 xcumul="x cumulative residuals",
                 ycumul="y cumulative residuals",
                 sum="sum of all residuals")
  avail <- unlist(optkey[names(opt)[unlist(opt)]])
  names(avail) <- NULL
  cat(paste("\t", paste(avail, collapse="\n\t"), "\n", sep=""))
  
  if(opt$sum) {
    xs <- x$sum
    netname <- "entire network"
    splat("sum of", typename, "in", netname, "=", signif(sum(xs$marksum),4))
    splat("length of", netname, "=", signif(xs$lengthL, 4))
    splat("quadrature total mass =", signif(xs$sumquad, 4))
  }
  if(opt$smooth) {
    splat("range of smoothed field = ", prange(signif(range(x$smooth$Z),4)))
    if(!is.null(sdp <- x$smooth$sdp))
      splat("Null standard deviation of smoothed Pearson residual field:",
            signif(sdp, 4))
  }
  return(invisible(NULL))
}
