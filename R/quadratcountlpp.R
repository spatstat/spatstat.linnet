#'
#' quadratcountlpp.R
#'
#' Quadrat counting on a network
#'
#' $Revision: 1.5 $ $Date: 2025/11/16 09:23:22 $

quadratcount.lpp <- function(X, ...,
                             nx=5, ny=nx,
                             xbreaks=NULL, ybreaks=NULL, left.open=TRUE, 
                             tess=NULL)  {
  stopifnot(is.lpp(X))
  L <- domain(X)
  cooX <- coords(X)
  if(missing(nx) && missing(ny) && is.null(xbreaks) && is.null(tess)) {
    ## default: each network segment is a tile
    tess <- as.lintess(L)
    segid <- factor(coords(X)$seg,
                    levels=seq_len(nsegments(L)),
                    labels=tilenames(tess))
    Xcount <- table(segid)
  } else if(is.null(tess)) {
    ## rectangular tiles
    if(!is.numeric(nx))
      stop("nx should be numeric")
    ## start with rectangular tessellation
    tess <- quadrats(Frame(L),
                     nx=nx, ny=ny, xbreaks=xbreaks, ybreaks=ybreaks)
    ## fast code for counting points in rectangular grid
    ## (yields two-dimensional table)
    Xcount <- rectquadrat.countEngine(cooX$x, cooX$y, tess$xgrid, tess$ygrid,
                                      left.open=left.open)
    ## now intersect rectangles with network
    tess <- intersect.lintess(tess, L)
  } else if(is.tess(tess)) {
    tess <- intersect.lintess(tess, L)
    ## classify points of X by tiles
    j <- lineartileindex(cooX$seg, cooX$tp, tess)
    ## count them (one-dimensional table)
    Xcount <- table(j)
  } else if(!inherits(tess, "lintess")) {
    stop(paste("Unrecognised format of tessellation argument", sQuote("tess")),
         call.=FALSE)
  }
  ## wrap up
  attr(Xcount, "tess") <- tess
  class(Xcount) <- unique(c("linearquadratcount", class(Xcount)))
  return(Xcount)
}

as.lintess.linearquadratcount <- function(X, ...) { attr(X, "tess") }

domain.linearquadratcount <- function(X, ...) { domain(as.lintess(X), ...) }

as.owin.linearquadratcount <- function(W, ..., fatal=TRUE) {
  as.owin(as.lintess(W), ..., fatal=fatal)
}

Window.linearquadratcount <- function(X, ...) { as.owin(X, ...) }

intensity.linearquadratcount <- function(X, ..., fun=FALSE) {
  tes <- as.lintess(X)
  len <- tile.lengths(tes)
  lam <- as.numeric(X)/len
  if(fun) lam <- as.linfun(tes, values=lam)
  return(lam)
}

shift.linearquadratcount <- function(X, ...) {
  Y <- X
  Y$domain <- L <- shift(X$domain, ...)
  Y <- putlastshift(Y, getlastshift(L))
  return(Y)
}

plot.linearquadratcount <- function(x, ..., main, add=FALSE,
                                    do.plot=TRUE, textargs=list()) {
  if(missing(main)) main <- short.deparse(substitute(x))
  ## plot tessellation
  tess <- as.lintess(x)
  cmap <- plot(tess, add=add, main=main, ..., do.plot=do.plot)
  if(!do.plot) return(invisible(cmap))
  ## annotate with counts
  Xcount <- as.integer(t(x))
  ## first find centroids of tiles
  len <- tile.lengths(tess)
  delta <- max(min(len)/10, sum(len)/1000)
  L <- as.linnet(tess)
  P <- pointsAlongNetwork(L, delta)
  ti <- lineartileindex(P$seg, P$tp, tess)
  x <- tapply(P$x, ti, mean)
  y <- tapply(P$y, ti, mean)
  do.call(text.default, append(list(x=x, y=y, labels=Xcount), textargs))
  return(invisible(cmap))
}

  



  
    
