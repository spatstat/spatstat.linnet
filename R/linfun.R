#
#   linfun.R
#
#   Class of functions of location on a linear network
#
#   $Revision: 1.18 $   $Date: 2023/05/02 08:02:49 $
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
  if(is.null(L))
    L <- domain(X)
  #' create template
  typical <- X(runiflpp(1, L), ...)
  if(length(typical) != 1)
    stop(paste("The function must return a single value",
               "when applied to a single point"))
  rule.eps <- match.arg(rule.eps)
  Y <- as.linim(typical, L, eps=eps, dimyx=dimyx, xy=xy,
                rule.eps=rule.eps,
                delta=delta, nd=nd)
  # extract coordinates of sample points along network
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

integral.linfun <- function(f, domain=NULL, weight=NULL, ..., delta, nd) {
  if(missing(delta)) delta <- NULL
  if(missing(nd)) nd <- NULL
  integral(as.linim(f, delta=delta, nd=nd), domain=domain, weight=weight, ...)
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
