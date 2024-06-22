#'
#'  evaluateNETcovariate.R
#'
#'  Twin of 'evaluatecovariates.R' for linear networks
#' 
#'  Evaluate a covariate along a network or at points on a network
#'
#'  $Revision: 1.7 $ $Date: 2024/06/22 10:12:53 $

evaluateNetCovariate <- function(covariate, locations, ...) {
  if(is.lpp(locations)) {
    evaluateNetCovariateAtPoints(covariate, locations, ...)
  } else {
    evaluateNetCovariateAlongNetwork(covariate, locations, ...)
  }
}

evaluateNetCovariateAlongNetwork <- function(covariate, locations, ...,
                                             types=NULL) {
  stopifnot(is.linnet(locations))
  L <- locations
  if(is.linim(covariate) || is.im(covariate) || inherits(covariate, "linfun")) {
    value <- as.linim(covariate, L=L, ...)
  } else if(is.imlist(covariate)) {
    value <- solapply(covariate, as.linim, L=L, ...)
  } else if(is.function(covariate)) {
    if(is.null(types) || length(formals(covariate)) <= 2L ||
       !any(c("m", "marks") %in% names(formals(covariate)))) {
      ## function (x,y) or function(x,y,...)  does not use mark value
      value <- as.linim(covariate, L=L, ...)
    } else {
      ## function f(x,y,m) or function(x,y,m, ...) expects mark value
      value <- solapply(types,
                        function(a) { as.linim(covariate, L=L, m=a, ...) })
    }
  } else if(is.list(covariate) && all(sapply(covariate, is.function))) {
    ## list of function(x,y) or list of function(x,y,..)
    value <- solapply(covariate, as.linim, L=L, ...)
  } else if(identical(covariate, "x")) {
    value <- as.linim(function(x,y){x}, L=L)
  } else if(identical(covariate, "y")) {
    value <- as.linim(function(x,y){y}, L=L)
  } else if(is.numeric(covariate) || is.factor(covariate)) {
    value <- solapply(covariate, as.linim, L=L)
  } else if(is.list(covariate) && all(lengths(covariate) == 1L) &&
            all(sapply(covariate, is.numeric) | sapply(covariate, is.factor))) {
    ## list of single values, associated with types
    if(length(types) > 0 && length(covariate) != length(types))
      stop(paste("Length of list of covariate values", paren(length(covariate)),
                 "does not match number of possible types in data",
                 paren(length(types))),
           call.=FALSE)
    value <- solapply(covariate, as.linim, L=L)
  } else stop("Format of covariate is not understood")
  if(!is.null(types)) {
    ## Ensure result is a solist of the right length
    if(is.linim(value)) {
      value <- rep(list(value), length(types))
    } else if(length(value) != length(types)) 
      warning("Mismatch between number of covariates and number of types",
              call.=FALSE)
    names(value) <- types
  }
  return(value)
}

evaluateNetCovariateAtPoints <- local({

  evaluateNetCovariateAtPoints <- function(covariate, locations, ...,
                                         allow.column=TRUE) {
    AvoidNames <- c("eps", "dimyx", "rule.eps", "delta", "nd", "types")
    stopifnot(is.lpp(locations))
    n <- npoints(locations)
    marx <- marks(locations) # may be null
    lev <- levels(marx) # may be null
    if(is.linim(covariate) || is.im(covariate)) {
      ## single pixel image 
      values <- safelooknet(covariate, locations)
    } else if(is.imlist(covariate)) {
      ## list of images for each type of point
      if(length(covariate) != length(lev)) 
        stop(paste("Number of covariate images", paren(length(covariate)),
                   "does not match number of possible types in data",
                   paren(length(lev))),
             call.=FALSE)
      values <- vector(mode="list", length=n)
      mm <- as.integer(marx)
      for(k in which(as.integer(table(mm)) > 0)) {
        relevant <- which(mm == k)
        values[relevant] <- safelooknet(covariate[[k]], locations[relevant])
      }
      values <- unlist(values)
    } else if(is.function(covariate)) {
      ## function(x,y) or function(x,y,m)
      xy <- coords(locations)
      if(length(formals(covariate)) <= 2L ||
         !any(c("m", "marks") %in% names(formals(covariate)))) {
        ## function does not use mark value
        values <- do.call.without(covariate,
                                  xy$x, xy$y, ...,
                                  avoid=AvoidNames)
      } else {
        ## function expects the mark values
        values <- do.call.without(covariate,
                                  xy$x, xy$y, marx, ...,
                                  avoid=AvoidNames)
      }
    } else if(is.list(covariate) && all(sapply(covariate, is.function))) {
      ## list of functions for each type of point
      if(length(covariate) != length(lev))
        stop(paste("Number of covariate functions", paren(length(covariate)),
                   "does not match number of possible types in data",
                   paren(length(lev))),
             call.=FALSE)
      values <- vector(mode="list", length=n)
      xy <- coords(locations)
      xx <- xy$x
      yy <- xy$y
      mm <- as.integer(marx)
      for(k in which(as.integer(table(mm)) > 0)) {
        relevant <- which(mm == k)
        values[relevant] <- do.call.without(covariate[[k]],
                                            xx[relevant], yy[relevant], ...,
                                            avoid=AvoidNames)
      }
      values <- unlist(values)
    } else if(is.numeric(covariate) || is.factor(covariate)) {
      ## numeric/categorical value or vector
      if(length(covariate) == 1L) {
        values <- rep.int(covariate, n)
      } else if(allow.column && length(covariate) == n) {
        ## NOTE
        values <- covariate 
      } else stop("Inappropriate length for covariate vector")
    } else if(is.list(covariate) && all(lengths(covariate) == 1L) &&
              (all(sapply(covariate, is.numeric)) ||
               all(sapply(covariate, is.factor)))) {
      ## list of single values, assumed to correspond to types
      if(length(covariate) != length(lev)) 
        stop(paste("Length of list of covariate values",
                   paren(length(covariate)),
                   "does not match number of possible types in data",
                   paren(length(lev))),
             call.=FALSE)
      values <- unlist(covariate[as.integer(marx)])
    } else if(identical(covariate, "x")) {
      values <- coords(locations)$x
    } else if(identical(covariate, "y")) {
      values <- coords(locations)$y
    } else 
      stop("Covariate should be an image, a function or a factor/numeric vector")
    return(values)
  }
  
  safelooknet <- function(Z, X) {
    if(is.linim(Z) && is.lpp(X)) {
      values <- Z[X, drop=FALSE]
      if(anyNA(values)) {
        bad <- is.na(values)
        values[bad] <- safelookup(as.im(Z), as.ppp(X))
      }
    } else {
      values <- safelookup(as.im(Z), as.ppp(X))
    }
    return(values)
  }

  evaluateNetCovariateAtPoints
})


