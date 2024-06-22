#'  exactlppm.R
#'
#'  Twin of 'exactppm.R' for linear networks
#' 
#'  An internal device to represent Poisson point process models on a network
#'  for which the MLE is computable exactly.
#'       - uniform intensity
#'       - intensity proportional to baseline
#'  These models are used mainly as a mathematical device
#'  in nonparametric methods such as 'rhohat'
#'  to represent the null/reference model,
#'  so that the code for nonparametric methods does not depend on 'lppm'
#'
#'  $Revision: 1.4 $ $Date: 2024/06/22 08:16:18 $

exactlppm <- function(X, baseline=NULL, ..., subset=NULL) {
  stopifnot(is.lpp(X))
  marx <- marks(X) # may be null
  lev <- levels(marx) # may be null
  if(is.null(subset)) {
    Xfit <- X
  } else {
    verifyclass(subset, "owin")
    Xfit <- X[subset]
  }
  L <- domain(X)
  if(is.null(baseline)) {
    #' stationary Poisson process
    beta <- intensity(Xfit)
  } else {
    #' Poisson process with intensity proportional to baseline
    if(is.linim(baseline)) {
      denom <- integral(baseline, domain=subset)
    } else if(is.list(baseline) &&
              all(sapply(baseline, is.linim)) &&
              is.multitype(X) &&
              length(baseline) == length(lev)) {
      denom <- sapply(baseline, integral, domain=subset)
    } else if(is.function(baseline)) {
      if(!is.multitype(X) || length(formals(baseline)) == 2) {
        ba <- as.linim(baseline, L=L, ...)
        denom <- integral(ba, domain=subset)
      } else {
        ba <- lapply(lev, function(z) { as.linim(baseline, L=L, z, ...) })
        denom <- sapply(ba, integral, domain=subset)
      }
    } else if(identical(baseline, "x")) {
      ba <- as.linim(function(x,y){x}, L=L, ...)
      denom <- integral(ba, domain=subset)
    } else if(identical(baseline, "y")) {
      ba <- as.linim(function(x,y){y}, L=L, ...)
      denom <- integral(ba, domain=subset)
    } else if(is.numeric(baseline) &&
              (length(baseline) == 1 ||
               is.multitype(X) && length(baseline) == length(lev))) {
      denom <- baseline * volume(L)
    } else stop("Format of 'baseline' is not understood")
    numer <- if(is.multitype(Xfit)) as.integer(table(marks(Xfit))) else npoints(Xfit)
    beta <- numer/denom
    if(length(beta) == length(lev))
      names(beta) <- lev
  }
  model <- list(X=X, baseline=baseline, subset=subset, beta=beta)
  class(model) <- c("exactlppm", class(model))
  return(model)
}

print.exactlppm <- function(x, ...) {
  with(x, {
    splat("Exactly-fitted point process model on a network")
    if(!is.multitype(X)) {
      if(is.null(baseline)) splat("Homogeneous intensity", signif(beta, 4)) else
                            splat("Intensity proportional to baseline", 
                                  paren(paste("proportionality constant",
                                              signif(beta, 4))))
    } else {
      lab <- levels(marks(X))
      if(is.null(baseline)) {
        splat("Homogeneous intensities:")
        splat(paste(paste(lab, signif(beta, 4), sep=": "), collapse=", "))
      } else {
        splat("Intensities proportional to baseline")
        splat("Proportionality constants:")
        splat(paste(paste(lab, signif(beta, 4), sep=": "), collapse=", "))
      }
    }
  })
  return(invisible(NULL))
}

is.poisson.exactlppm <- function(x) { TRUE }

is.stationary.exactlppm <- function(x) {
  is.null(x$baseline) || is.numeric(x$baseline)
}

predict.exactlppm <- function(object, locations=NULL, ...) {
  X        <- object$X
  beta     <- object$beta # numeric
  baseline <- object$baseline # covariate or NULL
  if(is.null(locations))
    locations <- domain(X)
  if(length(beta) > 1) {
    ## Intensities for different types
    ## Syntax of 'evaluateCovariate' requires a list in this case
    beta <- as.list(beta)
  }
  ## evaluate at desired locations
  Beta <- evaluateNetCovariate(beta, locations, ...)
  if(is.null(baseline)) {
    Lambda <- Beta
  } else {
    Baseline <- evaluateNetCovariate(baseline, locations, ...)
    if(is.im(Beta) || is.imlist(Beta) || is.im(Baseline) || is.imlist(Baseline)) {
      Lambda <- LinimListOp(Beta, Baseline, "*")
    } else {
      Lambda <- Beta * Baseline
    }
  }
  ## tidy
  if(is.imlist(Lambda)) {
    if(length(Lambda) == 1) { 
      Lambda <- Lambda[[1]]
    } else if(length(Lambda) == length(Beta)) {
      names(Lambda) <- names(beta)
    }
  }
  return(Lambda)
}

response.exactlppm <- function(object) { object$X }
