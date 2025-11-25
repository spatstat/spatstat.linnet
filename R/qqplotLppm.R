#'
#'    QQ plot of smoothed residual field against model 
#'   For point process model on a linear network
#' 
#'  qqplot.lppm()       QQ plot (including simulation)
#'
#'  $Revision: 1.2 $   $Date: 2025/11/24 08:56:11 $
#'

qqplot.lppm <- local({

  ## How to refit the model
  refit <- function(fit, pattern) {
    force(pattern)
    updateData(fit, pattern)
  }
  
  ## how to compute the residual field
  residualfield <- function(fit, ..., addtype=FALSE) {
    d <- diagnose(fit, which="smooth",
                  plot.it=FALSE, compute.cts=FALSE, compute.sd=FALSE,
                  check=FALSE, ...)
    result <- d$smooth$Z # smoothed residual field, pixel image on network
    if(addtype) {
      attr(result, "type") <- d$type
      attr(result, "typename") <- d$typename
    }
    return(result)
  }

  ## how to extract pixel values from image
  pixelvalues <- function(x) { as.numeric(x$v) }
  
  qqplot.lppm <-
    function(fit, nsim=100, expr=NULL, ..., type="raw", style="mean",
             fast=TRUE, verbose=TRUE, plot.it=TRUE, probs=NULL,
             saveall=FALSE,
             monochrome=FALSE,
             limcol=if(monochrome) "black" else "red",
             maxerr=max(100, ceiling(nsim/10)),
             envir.expr) {

    verifyclass(fit, "lppm")

    if(fast) {
      ## use coarse pixel grid
      oldnpixel <- spatstat.options("npixel")
      newnpixel <- pmin(40, oldnpixel)
      spatstat.options(npixel=newnpixel)
      on.exit(spatstat.options(npixel=oldnpixel))
    } 
    

    ################   How to evaluate residuals ##########################
  
    ## Data values
    datim <- residualfield(fit, type=type, ..., addtype=TRUE)
    ## type of residuals (partially matched and validated by diagnose.ppm)
    type <- attr(datim, "type")
    typename <- attr(datim, "typename") 

    ## Extract pixel values
    dat <- pixelvalues(datim) #' all pixel values including NA's
    npix <- length(dat)
    nfinite <- if(anyNA(dat)) sum(!is.na(dat)) else npix
      
    ## pixel weights (= length of network in each pixel)
    L <- as.linnet(data.lppm(fit))
    pixwt <- pixellate(L, Window(datim))
    pixwt <- pixelvalues(pixwt)
    pixwt <- pixwt/sum(pixwt, na.rm=TRUE)

    if(length(pixwt) != npix)
      stop("Internal error: pixel weights have wrong length", call.=FALSE)

    ## Quantiles of the residual field will be computed.
    if(is.null(probs))
      probs <- ppoints(min(100, nfinite), 3/8)
    
    ##################  How to perform simulations?  #######################

    ## envir.call <- sys.parent()
    envir.here <- sys.frame(sys.nframe())

    ## extract.from.list <- FALSE
    inext <- 0 # to placate package checker
    dont.complain.about(inext)
    
    if(is.null(expr)) {

      ## We will simulate from the fitted model 'nsim' times
      ## and refit the model to these simulations
      simsource <- "fit"
      how.simulating <- "simulating from fitted model" 

      ## expression to be evaluated each time
      expr <- expression( refit(fit, simulate(fit, nsim=1, drop=TRUE)) )
      envir.expr <- envir.here
      
    } else if(is.expression(expr)) {
      simsource <- "expr"
      how.simulating <- paste("evaluating", sQuote("expr"))  
      if(missing(envir.expr) || is.null(envir.expr))
        envir.expr <- parent.frame()
    } else if(inherits(expr, "envelope")) {
      simpat <- attr(expr, "simpatterns")
      if(!is.null(simpat) && all(sapply(simpat, is.ppp))) {
        expr <- expression(simpat[[inext]])
        envir.expr <- envir.here
        dont.complain.about(simpat)
        simsource <- "list"
        how.simulating <- "extracting point pattern from list"
      } else stop(paste("Argument", sQuote("expr"),
                        "is an envelope object,",
                        "but does not contain point patterns"),
                  call.=FALSE)
    } else if(is.list(expr) && all(sapply(expr, is.ppp))) {
      simpat <- expr
      expr <- expression(simpat[[inext]])
      envir.expr <- envir.here
      dont.complain.about(simpat)
      simsource <- "list"
      how.simulating <- "extracting point pattern from list"
    } else stop(paste(sQuote("expr"),
                      "should be an expression, or an envelope object,",
                      "or a list of point patterns"),
                call.=FALSE)

    exprstring <- if(simsource == "expr") deparse(expr) else NULL

    ######  Perform simulations
    if(verbose) {
      cat(paste("Simulating", nsim, "realisations... "))
      pstate <- list()
    }
    sim <- matrix(, nrow=npix, ncol=nsim)
    simul.sizes <- numeric(nsim)
    isim <- 0
    ierr <- 0
    repeat {
      inext <- isim + 1
      ## protect from randomly-generated crashes in gam
      ##      ei <- try(eval(expr, envir=envir.expr), silent=!verbose)
      ei <- eval(expr, envir=envir.expr)
      if(inherits(ei, "try-error")) {
        ## error encountered in evaluating 'expr'
        ierr <- ierr + 1
        if(ierr > maxerr) 
          stop(paste("Exceeded maximum of", maxerr,
                     "failures in", how.simulating,
                     "after generating only", isim, "realisations"))
        else break
      } else {
        ## simulation successful
        isim <- isim + 1
        fiti <- 
          if(simsource == "fit")
            ei
          else if(is.lppm(ei))
            ei
          else if(is.lpp(ei))
            refit(fit, ei)
          else
            stop("result of eval(expr) is not an lppm or lpp object")
        ## diagnostic info
        simul.sizes[isim] <- npoints(data.lppm(fiti))
        ## compute residual field
        resi <- residualfield(fiti, type=type, ...)
        resi <- pixelvalues(resi)
        if(length(resi) != npix)
          stop("Internal error: wrong number of pixels in simulated residual",
               call.=FALSE)
        sim[,isim] <- resi
        if(verbose) 
          pstate <- progressreport(isim, nsim, state=pstate)
        if(isim >= nsim)
          break
      }
    }

    ###### Report diagnostics
    if(ierr > 0)
      cat(paste("\n\n**Alert:",
                ierr, "failures occurred in", how.simulating, "\n\n"))
    nempty <- sum(simul.sizes == 0)
    if(nempty > 0)
      cat(paste("\n\n**Alert:",
                nempty, "out of", nsim,
                "simulated patterns were empty.\n\n"))
    else
      cat(paste("\nDiagnostic info:\n",
                "simulated patterns contained an average of",
                mean(simul.sizes), "points.\n"))
    if(nempty == nsim)
      warning("All simulated patterns were empty")
    ############ Plot them
    switch(style,
           classical = {
             rr <- range(c(dat,sim))
             datquant <- weighted.quantile(x=dat,
                                           w=pixwt,
                                           probs=probs)
             simquant <- weighted.quantile(x=as.numeric(sim),
                                           w=as.numeric(pixwt[row(sim)]),
                                           probs=probs)
             result <- list(x = simquant, y = datquant)
             if(plot.it) {
               plot(simquant, datquant,
                    xlim=rr, ylim=rr, asp=1.0,
                    xlab="Quantiles of simulation",
                    ylab="Quantiles of data")
               title(sub=typename)
               abline(0,1, lty=2)
             }
             result <- list(x = simquant,
                            y = datquant,
                            data=dat,
                            sim=sim,
                            probs=probs,
                            xlim=rr,
                            ylim=rr,
                            xlab="Quantiles of simulation",
                            ylab="Quantiles of data",
                            rtype=type,
                            typename=typename,
                            nsim=nsim,
                            fit=fit,
                            expr=exprstring,
                            simsource = simsource
                            )
           },
           mean = {
             ## compute quantiles corresponding to probabilities p[i]
             ## separately in each realisation.
             browser()
             datquant <- weighted.quantile(x=dat,
                                           w=pixwt,
                                           probs=probs)
             if(verbose) cat("Calculating quantiles...")
             qsim <- apply(sim, 2,
                           weighted.quantile,
                           w=pixwt, probs=probs)
             if(verbose) cat("averaging...")
             ## sample mean of each quantile
             simquant <- apply(qsim, 1, mean, na.rm=TRUE)
             ## et cetera
             varq <- apply(qsim, 1, var, na.rm=TRUE)
             sdq <- sqrt(varq)
             q.025 <- apply(qsim, 1, quantile, probs=0.025, na.rm=TRUE)
             q.975 <- apply(qsim, 1, quantile, probs=0.975, na.rm=TRUE)
             rr <- range(c(simquant,dat), na.rm=TRUE)
             if(verbose) cat("..Done.\n")
             if(plot.it) {
               plot(simquant, datquant,
                    xlim=rr, ylim=rr, asp=1.0,
                    xlab="Mean quantiles of simulations",
                    ylab="Quantiles of data")
               title(sub=typename)
               abline(0,1)
               lines(simquant, q.025, lty=2, col=limcol)
               lines(simquant, q.975, lty=2, col=limcol)
             }
             result <- list(x     = simquant,
                            y     = datquant,
                            sdq   = sdq,
                            q.025 = q.025,
                            q.975 = q.975,
                            data  = dat,
                            sim   = sim,
                            probs = probs,
                            xlim  = rr,
                            ylim  = rr,
                            xlab  = "Mean quantiles of simulations",
                            ylab  = "Quantiles of data",
                            rtype = type,
                            typename = typename,
                            nsim     = nsim,
                            fit      = fit,
                            expr     = exprstring,
                            simsource = simsource)
           },
           stop(paste("Unrecognised option for", sQuote("style")))
           )

    ## Throw out baggage if not wanted         
    if(!saveall) {
      result$fit <- summary(fit)
      result$sim <- NULL
    }
         
    ## reset npixel is done by on.exit()
    ##
    class(result) <- unique(c("qqlppm", class(result)))
    return(invisible(result))
  }

  qqplot.lppm

})


plot.qqlppm <- local({

  plot.qqlppm <- function(x, ..., limits=TRUE,
                         monochrome=spatstat.options('monochrome'),
                         limcol=if(monochrome) "black" else "red") {
    stopifnot(inherits(x, "qqlppm"))
    default.type <- if(length(x$x) > 150) "l" else "p"
    do.call(myplot,
            resolve.defaults(list(quote(x), ..., type=default.type,
                                  limits=limits, limcol=limcol)))
    return(invisible(x))
  }

  myplot <- function(object,
                     xlab = object$xlab, ylab = object$ylab,
                     xlim = object$xlim, ylim = object$ylim,
                     asp = 1,
                     type = "l",
                     ..., limits=TRUE, limcol="red") {
    plot(object$x, object$y, xlab = xlab, ylab = ylab,
         xlim = xlim, ylim = ylim, asp = asp, type = type, ...)
    abline(0, 1)
    
    if(limits) {
      if(!is.null(object$q.025))
        lines(object$x, object$q.025, lty = 2, col=limcol)
      if(!is.null(object$q.975))
        lines(object$x, object$q.975, lty = 2, col=limcol)
    }
    typename <- object$typename %orifnull% paste(object$rtype, "residuals")
    title(sub=typename)
  }

  plot.qqlppm
})


print.qqlppm <- function(x, ...) {
  stopifnot(inherits(x, "qqlppm"))
  splat("Q-Q plot of point process residuals",
        "of type", sQuote(x$rtype), "\n",
        "based on", x$nsim, "simulations")
  simsource <- x$simsource
  if(is.null(simsource)) # old version
    simsource <- if(x$simulate.from.fit) "fit" else "expr"
  switch(simsource,
         fit = {
           fit  <- x$fit
           sumfit <- if(is.lppm(fit)) summary(fit)
                     else if(inherits(fit, "summary.lppm")) fit
                     else list(name="(unrecognised format)")
           splat("\nSimulations from fitted model:", sumfit$name)
         },
         expr = {
           splat("Simulations obtained by evaluating the following expression:")
           print(x$expr)
         },
         list = {
           splat("Simulated point patterns were provided in a list")
         })
  invisible(NULL)
}

