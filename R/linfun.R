#
#   linfun.R
#
#   Class of functions of location on a linear network
#
#   $Revision: 1.29 $   $Date: 2026/01/21 06:26:39 $
#

linfun <- function(f, L) {
  stopifnot(is.function(f))
  stopifnot(inherits(L, "linnet"))
  fargs <- names(formals(f))
  needargs <- c("x", "y", "seg", "tp")
  if(!all(needargs %in% fargs))
    stop(paste("Function must have formal arguments",
               commasep(sQuote(needargs))),
         call.=FALSE)
  otherfargs <- setdiff(fargs, needargs)
  g <- function(...) {
    argh <- list(...)
    extra <- names(argh) %in% otherfargs
    if(!any(extra)) {
      X <- as.lpp(..., L=L)
      value <- do.call(f, as.list(coords(X)))
    } else {
      extrargs <- argh[extra]
      mainargs <- argh[!extra]
      X <- do.call(as.lpp, append(mainargs, list(L=L)))
      value <- do.call(f, append(as.list(coords(X)), extrargs))
    }
    return(value)
  }
  class(g) <- c("linfun", class(g))
  attr(g, "L") <- L
  attr(g, "f") <- f
  unitname(g) <- unitname(L)
  return(g)
}

print.linfun <- function(x, ...) {
  L <- as.linnet(x)
  if(!is.null(explain <- attr(x, "explain"))) {
    explain(x)
  } else {
    splat("Function on linear network:")
    print(attr(x, "f"), ...)
    splat("Function domain:")
    print(L)
  }
  invisible(NULL)
}

summary.linfun <- function(object, ...) { print(object, ...) }

as.linim.linfun <- function(X, L=domain(X),
                            ...,
                            eps = NULL, dimyx = NULL, xy = NULL,
                            rule.eps=c("adjust.eps",
                                        "grow.frame", "shrink.frame"),
                            delta=NULL, nd=NULL) {
  if(is.linim(L) && !is.null(attr(L, "df"))) {
    #' use L as template object
    Y <- L
    L <- as.linnet(Y)
  } else {
    #' create template object
    typical <- X(runiflpp(1, L), ...)
    if(length(typical) != 1)
      stop(paste("The function must return a single value",
                 "when applied to a single point"))
    rule.eps <- match.arg(rule.eps)
    Y <- as.linim(typical, L, eps=eps, dimyx=dimyx, xy=xy,
                  rule.eps=rule.eps,
                  delta=delta, nd=nd)
  }
  # extract coordinates of sample points along network in template 
  df <- attr(Y, "df")
  coo <- df[, c("x", "y", "mapXY", "tp")]
  colnames(coo)[3L] <- "seg"
  # evaluate function at sample points
  vals <- do.call(X, append(as.list(coo), list(...)))
  # write values in data frame
  df$values <- vals
  attr(Y, "df") <- df
  #' overwrite values in pixel array
  Y$v[] <- NA
  pix <- nearest.raster.point(df$xc, df$yc, Y)
  Y$v[cbind(pix$row, pix$col)] <- vals
  #'
  return(Y)
}

as.data.frame.linfun <- function(x, ...) {
  as.data.frame(as.linim(x, ...))
}

as.function.linim <- function(x, ...) { as.linfun.linim(x, ...) }

as.linfun.linim <- function(X, ...) {
  trap.extra.arguments(..., .Context="as.linfun.linim")
  ## extract info
  L <- as.linnet(X)
  df <- attr(X, "df")
  ## function values and corresponding locations
  values <- df$values
  locations <- with(df, as.lpp(x=x, y=y, seg=mapXY, tp=tp, L=L))
  ## Function that maps any spatial location to the nearest data location 
  nearestloc <- nnfun(locations)
  ## Function that reads value at nearest data location
  f <- function(x, y, seg, tp) {
    values[nearestloc(x,y,seg,tp)]
  }
  g <- linfun(f, L)
  return(g)
}

plot.linfun <- function(x, ..., L=NULL, main) {
  if(missing(main)) main <- short.deparse(substitute(x))
  if(is.null(L)) L <- as.linnet(x)
  argh <- list(...)
  fargnames <- get("otherfargs", envir=environment(x))
  resolution <- c("eps", "dimyx", "xy", "rule.eps", "delta", "nd")
  convert <- names(argh) %in% c(fargnames, resolution)
  Z <- do.call(as.linim, append(list(x, L=L), argh[convert]))
  rslt <- do.call(plot.linim, append(list(Z, main=main), argh[!convert]))
  return(invisible(rslt))
}

as.owin.linfun <- function(W, ...) {
  as.owin(as.linnet(W))
}

domain.linfun <- as.linnet.linfun <- function(X, ...) {
  attr(X, "L")
}

as.function.linfun <- function(x, ...) {
  nax <- names(attributes(x))
  if(!is.null(nax)) {
    retain <- (nax == "srcref")
    attributes(x)[!retain] <- NULL
  }
  return(x)
}

integral.linfun <- function(f, domain=NULL, weight=NULL, ..., exact=FALSE,
                            delta, nd) {
  if(exact) {
    ## decide whether exact calculation is feasible
    weightok <- is.null(weight) || is.function(weight)
    if(!weightok) {
      warning("exact=TRUE ignored because 'weight' must be discretised",
              call.=FALSE)
      exact <- FALSE
    }
    domainok <- is.null(domain) || inherits(domain, c("lintess","tess","owin"))
    if(!domainok) {
      warning("exact=TRUE ignored because 'domain' must be discretised",
              call.=FALSE)
      exact <- FALSE
    }
  }
  if(!exact) {
    ## discrete calculation using integral.linim
    if(missing(delta)) delta <- NULL
    if(missing(nd)) nd <- NULL
    fim <- as.linim(f, delta=delta, nd=nd)
    result <- integral(fim, domain=domain, weight=weight, ...)
    return(result)
  }
  ## 'exact' integration using integrate()
  L <- domain(f)
  if(is.function(weight)) {
    ## absorb weight into f and continue
    f.orig      <- f
    weight.orig <- weight
    wtf <- if(inherits(weight,"linfun")) weight.orig else linfun(weight.orig, L)
    fw <- function(x, y, seg, tp, ...) {
      f.orig(x, y, seg, tp, ...) * wtf(x, y, seg, tp, ...)
    }
    f <- linfun(fw, L)
    weight <- NULL
  }
  ## map network to an interval [0, T] on the real line
  m <- mapLinnetToLine(L)
  totlen    <- m$totlen
  netcoords <- m$netcoords
  ## convert the integrand to a function on [0, T]
  g <- function(x, ...) { do.call(f, append(netcoords(x), list(...))) }
  eps <- .Machine$double.eps
  ## preprocess different kinds of 'domain'
  if(waswin <- is.owin(domain)) {
    A <- domain
    B <- Frame(L)
    notA <- setminus.owin(B, A)
    domain <- tess(tiles=list(inside=A, outside=notA), window=B)
  }
  if(is.tess(domain))
    domain <- intersect.lintess(domain, L)
  ## ----------------- start computing --------------------------
  if(is.null(domain)) {
    ## use integrate() on each segment for accuracy
    nseg <- nsegments(L)
    if(nseg == 0) return(0)
    cumlen <- m$cumlen
    for(seg in seq_len(nseg)) {
      z <- integrate(g,
                     lower = if(seg == 1L) 0 else cumlen[seg-1],
                     upper = cumlen[seg] - eps,
                     ...)$value
      results <- if(seg == 1L) z else c(results, z)
    }
    result <- sum(results)
    return(result)
  } else if(inherits(domain, "lintess")) {
    ## use integrate() on each tile fragment
    df <- domain$df
    npieces <- nrow(df)
    ## Each row of df represents a tile fragment
    ## Map each fragment to an interval of the real line
    flatcoord <- m$flatcoord
    df$lower <- with(df, flatcoord(seg=seg, tp=t0))
    df$upper <- with(df, flatcoord(seg=seg, tp=t1))
    ## use integrate() on each fragment
    integ <- numeric(npieces)
    for(i in seq_len(npieces)) 
      integ[i] <- integrate(g,
                            lower=df$lower[i]+eps,
                            upper=df$upper[i]-eps,
                            ...)$value
    ## sum contributions from fragments for each tile
    result <- tapplysum(integ, list(tile=df$tile))
    if(waswin) {
      result <- result[[1L]]
    } else {
      names(result) <- tilenames(domain)
    }
    return(result)
  } else stop("Internal error - unrecognised kind of domain", call.=FALSE)
}

as.linfun <- function(X, ...) {
  UseMethod("as.linfun")
}

as.linfun.linfun <- function(X, ...) { return(X) }

as.linfun.linnet <- function(X, ..., values=marks(X)) {
  if(missing(values) || is.null(values)) {
    values <- marks(X) %orifnull% seq_len(nsegments(X))
    if(!is.null(dim(values))) values <- values[,1]
  } 
  stopifnot(length(values) == nsegments(X))
  f <- function(x, y, seg, tp) { values[seg] }
  g <- linfun(f, unmark(X))
  return(g)
}

as.function.linfun <- function(x, ...) { return(x) }

as.function.linnet <- function(x, ...) { as.linfun.linnet(x, ...) }

levelset.linfun <- function(X, thresh, compare="<=", ...) {
  levelset(as.linim(X), thresh=thresh, compare=compare, ...)
}
