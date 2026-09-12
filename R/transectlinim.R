#'
#' transectlinim.R
#'
#' Transect of image on a network
#'
#' $Revision: 1.6 $ $Date: 2026/09/09 05:01:33 $
#'

transect.linim <- function(X, ..., path=NULL, click=FALSE, add=FALSE,
                           nsample=512, Xname=NULL) {
  if(is.null(Xname)) {
    Xname <- short.deparse(substitute(X))
    Xname <- sensiblevarname(Xname, "X")
  }
  stopifnot(is.linim(X))
  click <- isTRUE(click)
  check.1.integer(nsample)
  L <- unmark(as.linnet(X))
  V <- unmark(vertices(L))
  #' determine path through network
  co <- getNetTransectCoords(L, path=path,
                             click=click, add=add, nsample=nsample)
  ## look up values
  with(co, {
    SamplePoints <- as.lpp(x=x, y=y, seg=seg, tp=tp, L=L)
    values <- X[SamplePoints, drop=FALSE]
    ## package into fv object
    df <- data.frame(t=t, v=values)
    colnames(df) <- c(tname, Xname)
    result <- fv(df,
                 argu = tname,
                 ylab = substitute(Xname(tname),
                                   list(Xname=as.name(Xname),
                                        tname=as.name(tname))),
                 valu = Xname,
                 labl = c(tname, paste0("%s", paren(tname))),
                 desc = c(tdescrip, "pixel value of %s"),
                 unitname = tunits,
                 fname = Xname)
    ## save the path
    attr(result, "path") <- journey
    attr(result, "lines") <- L$lines[segmap]
    attr(result, "SamplePoints") <- SamplePoints
    return(result)
  })
}

transect.linfun <- function(X, ..., path=NULL, click=FALSE, add=FALSE,
                           nsample=512, Xname=NULL) {
  if(is.null(Xname)) {
    Xname <- short.deparse(substitute(X))
    Xname <- sensiblevarname(Xname, "X")
  }
  stopifnot(inherits(X, "linfun"))
  click <- isTRUE(click)
  check.1.integer(nsample)
  L <- unmark(as.linnet(X))
  V <- unmark(vertices(L))
  #' determine path through network
  co <- getNetTransectCoords(L, path=path,
                             click=click, add=add, nsample=nsample)
  ## evaluate function at sample points
  with(co, {
    SamplePoints <- as.lpp(x=x, y=y, seg=seg, tp=tp, L=L)
    values <- X(SamplePoints)
    ## package into fv object
    df <- data.frame(t=t, v=values)
    colnames(df) <- c(tname, Xname)
    result <- fv(df,
                 argu = tname,
                 ylab = substitute(Xname(tname),
                                   list(Xname=as.name(Xname),
                                        tname=as.name(tname))),
                 valu = Xname,
                 labl = c(tname, paste0("%s", paren(tname))),
                 desc = c(tdescrip, "function value of %s"),
                 unitname = tunits,
                 fname = Xname)
    ## save the path
    attr(result, "path") <- journey
    attr(result, "lines") <- L$lines[segmap]
    attr(result, "SamplePoints") <- SamplePoints
    return(result)
  })
}

getNetTransectCoords <- function(L, ...,
                                 path=NULL, click=FALSE, add=FALSE,
                                 nsample=512) {
  stopifnot(is.linnet(L))
  V <- unmark(vertices(L))
  nV <- npoints(V)
  if(!is.null(path)) {
    ## 'path' is given -  should be a sequence of vertices
    stopifnot(is.numeric(path))
    stopifnot(length(path) >= 2)
    path <- as.integer(path)
    if(min(path) < 1 || max(path) > nV)
      stop(paste("Argument 'path' should contain indices between 1 and",
                 nV, "(= number of vertices of network)"),
           call.=FALSE)
  } else if(click) {
    ## interactive selection of a sequence of vertices
    if(add) 
      plot(L, main=c("select vertices along path",
                     "tap ESC to end"))
    path <- integer(0)
    repeat {
      p <- identify(V, n=1, plot=FALSE)
      if(length(p)) {
        plot(V[p], add=TRUE, pch=16)
        path <- c(path, p)
      } else break
    }
    if(length(path) < 2)
      stop("At least two vertices must be selected", call.=FALSE)
  } else {
    ## no path information --- choose two maximally distant vertices
    L <- as.linnet(L, sparse=FALSE)
    d <- L$dpath
    ok <- is.finite(d)
    if(is.null(d) || !any(ok))
      stop("Unable to extract path distances", call.=FALSE)
    dmax <- max(d[ok])
    ij <- which(ok & (d == max(d)), arr.ind=TRUE)
    path <- as.integer(ij[1,])
  }
  #' fill in missing steps 
  VP <- lpp(V, L)
  nP <- length(path)
  journey <- path[1L]
  for(i in 2:nP) {
    vprev <- path[i-1L]
    vnext <- path[i]
    ## find vertices along shortest path from previous vertex to next vertex
    m <- suppressWarnings(shortestpath(VP, vprev, vnext))
    midsteps <- attr(m, "steps")
    journey <- c(journey, midsteps, vnext)
  }
  ## remove repeats
  journey <- journey[c(TRUE, diff(journey) != 0)]
  ## construct sample points
  nJ <- length(journey)
  vfrom  <- journey[-nJ]
  vto    <- journey[-1L]
  leglen <- L$dpath[cbind(vfrom, vto)]
  nlegs  <- length(leglen)
  dJ     <- sum(leglen)
  cumlen <- cumsum(c(0,leglen))
  ## arc length coordinate
  s <- seq(0, dJ, length.out=nsample)
  ## which leg of journey
  legid <- findInterval(s,
                        cumlen,
                        rightmost.closed=TRUE, all.inside=TRUE,
                        checkNA=FALSE, checkSorted=FALSE)
  ## local coordinates on network (up to tp <-> 1-tp )
  tp <- (s-cumlen[legid])/leglen[legid]
  tp[!is.finite(tp) | tp < 0 | tp > 1] <- 0.5
  ## map legs of journey to segments of network
  segmap <- integer(nlegs)
  revmap <- logical(nlegs)
  for(leg in seq_len(nlegs)) {
    candidates <- which(L$from == vfrom[leg] & L$to == vto[leg])
    if(length(candidates)) {
      segmap[leg] <- min(candidates)
    } else {
      candidates <- which(L$to == vfrom[leg] & L$from == vto[leg])
      if(length(candidates)) {
        segmap[leg] <- min(candidates)
        revmap[leg] <- TRUE
      } else stop("Internal error: cannot find segment", call.=FALSE)
    }
  }
  seg <- segmap[legid]
  ## fix coordinates tp <-> 1 - tp
  if(any(revmap)) {
    reverse <- revmap[legid]
    tp <- ifelse(reverse, 1 - tp, tp)
  }
  ## spatial coordinates
  x <- V$x[vfrom[legid]] * (1-tp) + V$x[vto[legid]] * tp
  y <- V$y[vfrom[legid]] * (1-tp) + V$y[vto[legid]] * tp
  ## pack up
  result <- list(t=s, # cumulative arc length
                 x=x, y=y, seg=seg, tp=tp, # local coordinates
                 tname="s",
                 tdescrip="distance along transect",
                 tunits=unitname(L),
                 journey = journey,
                 segmap = segmap)
  return(result)
}
