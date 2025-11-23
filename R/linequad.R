#
# linequad.R
#
#  $Revision: 1.24 $ $Date: 2025/11/23 03:49:17 $
#
# create quadscheme for a pattern of points lying *on* line segments

linequad <- function(X, Y, ..., eps=NULL, nd=1000, random=FALSE) {
  epsgiven <- !is.null(eps)
  if(is.lpp(X)) {
    # extract local coordinates from lpp object
    coo <- coords(X)
    mapXY <- coo$seg
    tp    <- coo$tp
    Xproj <- as.ppp(X)
    if(!missing(Y) && !is.null(Y))
      warning("Argument Y ignored when X is an lpp object")
    Y <- as.linnet(X)
  } else if(is.ppp(X)) {
    ## project data points onto segments
    if(missing(Y))
      stop("Argument Y is required when X is a point pattern", call.=FALSE)
    stopifnot(is.psp(Y) || is.linnet(Y))
    v <- project2segment(X, as.psp(Y))
    Xproj <- v$Xproj
    mapXY <- v$mapXY
    tp    <- v$tp
  } else stop("X should be a point pattern object of class lpp or ppp")
  
  # handle multitype
  ismulti <- is.multitype(X)
  if(is.marked(X) && !ismulti)
    stop("Not implemented for marked patterns")
  if(ismulti) {
    marx <- marks(X)
    flev <- factor(levels(marx))
  }
  #
  win <- as.owin(Y)
  len <- lengths_psp(as.psp(Y))
  nseg <- length(len)
  if(is.null(eps)) {
    stopifnot(is.numeric(nd) && length(nd) == 1L & is.finite(nd) && nd > 0)
    eps <- sum(len)/nd
  } else
  stopifnot(is.numeric(eps) && length(eps) == 1L && is.finite(eps) && eps > 0)
  ##
  if(is.lpp(X) && spatstat.options('Clinequad')) {
    L <- as.linnet(X)
    W <- Frame(L)
    V <- vertices(L)
    nV <- npoints(V)
    coordsV <- coords(V)
    coordsX <- coords(X)
    nX <- npoints(X)
    ooX <- order(coordsX$seg)
    ndumeach <- ceiling(len/eps) + 1L
    ndummax <- sum(ndumeach)
    maxdataperseg <- max(table(factor(coordsX$seg, levels=1:nsegments(L))))
    maxscratch <- max(ndumeach) + maxdataperseg
    if(!ismulti) {
      if(!random) {
        z <- .C(SL_Clinequad,
                ns    = as.integer(nseg),
                from  = as.integer(L$from-1L),
                to    = as.integer(L$to-1L), 
                nv    = as.integer(nV),
                xv    = as.double(coordsV$x),
                yv    = as.double(coordsV$y), 
                eps   = as.double(eps),
                ndat  = as.integer(nX),
                sdat  = as.integer(coordsX$seg[ooX]-1L),
                tdat  = as.double(coordsX$tp[ooX]),
                wdat  = as.double(numeric(nX)),
                ndum  = as.integer(integer(1L)),
                xdum  = as.double(numeric(ndummax)),
                ydum  = as.double(numeric(ndummax)),
                sdum  = as.integer(integer(ndummax)),
                tdum  = as.double(numeric(ndummax)),
                wdum  = as.double(numeric(ndummax)),
                maxscratch = as.integer(maxscratch),
                PACKAGE="spatstat.linnet")
      } else {
        z <- .C(SL_ClineRquad,
                ns    = as.integer(nseg),
                from  = as.integer(L$from-1L),
                to    = as.integer(L$to-1L), 
                nv    = as.integer(nV),
                xv    = as.double(coordsV$x),
                yv    = as.double(coordsV$y), 
                eps   = as.double(eps),
                ndat  = as.integer(nX),
                sdat  = as.integer(coordsX$seg[ooX]-1L),
                tdat  = as.double(coordsX$tp[ooX]),
                wdat  = as.double(numeric(nX)),
                ndum  = as.integer(integer(1L)),
                xdum  = as.double(numeric(ndummax)),
                ydum  = as.double(numeric(ndummax)),
                sdum  = as.integer(integer(ndummax)),
                tdum  = as.double(numeric(ndummax)),
                wdum  = as.double(numeric(ndummax)),
                maxscratch = as.integer(maxscratch),
                PACKAGE="spatstat.linnet")
      }
      ## dummy points
      ndum <- z$ndum
      seqdum <- seq_len(ndum)
      wdum <- z$wdum[seqdum]
      sdum <- z$sdum[seqdum] + 1L
      tdum <- z$tdum[seqdum]
      xdum <- z$xdum[seqdum]
      ydum <- z$ydum[seqdum]
      dum <- ppp(xdum, ydum, window=W, check=FALSE)
      ## data points (original ordering)
      wdat <- numeric(nX)
      wdat[ooX] <- z$wdat
      sdat <- coordsX$seg
      tdat <- coordsX$tp
      dat <- as.ppp(X)
    } else {
      ntypes <- length(flev)
      ndummax <- ntypes * (ndummax + nX)
      maxscratch <- ntypes * maxscratch
      if(!random) {
        z <- .C(SL_ClineMquad,
                ns    = as.integer(nseg),
                from  = as.integer(L$from-1L),
                to    = as.integer(L$to-1L), 
                nv    = as.integer(nV),
                xv    = as.double(coordsV$x),
                yv    = as.double(coordsV$y), 
                eps   = as.double(eps),
                ntypes = as.integer(ntypes),
                ndat  = as.integer(nX),
                xdat  = as.double(coordsX$x),
                ydat  = as.double(coordsX$y),
                mdat  = as.integer(as.integer(marx)-1L),
                sdat  = as.integer(coordsX$seg[ooX]-1L),
                tdat  = as.double(coordsX$tp[ooX]),
                wdat  = as.double(numeric(nX)),
                ndum  = as.integer(integer(1L)),
                xdum  = as.double(numeric(ndummax)),
                ydum  = as.double(numeric(ndummax)),
                mdum  = as.integer(integer(ndummax)),
                sdum  = as.integer(integer(ndummax)),
                tdum  = as.double(numeric(ndummax)),
                wdum  = as.double(numeric(ndummax)),
                maxscratch = as.integer(maxscratch),
                PACKAGE="spatstat.linnet")
      } else {
        z <- .C(SL_ClineRMquad,
                ns    = as.integer(nseg),
                from  = as.integer(L$from-1L),
                to    = as.integer(L$to-1L), 
                nv    = as.integer(nV),
                xv    = as.double(coordsV$x),
                yv    = as.double(coordsV$y), 
                eps   = as.double(eps),
                ntypes = as.integer(ntypes),
                ndat  = as.integer(nX),
                xdat  = as.double(coordsX$x),
                ydat  = as.double(coordsX$y),
                mdat  = as.integer(as.integer(marx)-1L),
                sdat  = as.integer(coordsX$seg[ooX]-1L),
                tdat  = as.double(coordsX$tp[ooX]),
                wdat  = as.double(numeric(nX)),
                ndum  = as.integer(integer(1L)),
                xdum  = as.double(numeric(ndummax)),
                ydum  = as.double(numeric(ndummax)),
                mdum  = as.integer(integer(ndummax)),
                sdum  = as.integer(integer(ndummax)),
                tdum  = as.double(numeric(ndummax)),
                wdum  = as.double(numeric(ndummax)),
                maxscratch = as.integer(maxscratch),
                PACKAGE="spatstat.linnet")
      }
      ## dummy points
      ndum <- z$ndum
      seqdum <- seq_len(ndum)
      wdum <- z$wdum[seqdum]
      sdum <- z$sdum[seqdum] + 1L
      tdum <- z$tdum[seqdum]
      mdum <- factor(z$mdum[seqdum] + 1L, labels=flev)
      xdum <- z$xdum[seqdum]
      ydum <- z$ydum[seqdum]
      dum <- ppp(xdum, ydum, marks=mdum, window=W, check=FALSE)
      ## data points (original order)
      wdat <- numeric(nX)
      wdat[ooX] <- z$wdat
      sdat <- coordsX$seg
      tdat <- coordsX$tp
      dat <- as.ppp(X)
    }      
  } else {
    ## Interpreted code
    W <- Frame(X)
    xX <- coords(X)$x
    yX <- coords(X)$y
    mX <- if(ismulti) marx else NULL
    xdat <- xdum <- ydat <- ydum <- wdat <- wdum <- tdat <- tdum <- numeric(0)
    jdat <- sdat <- sdum <- integer(0)
    mdat <- mdum <- if(ismulti) marx[integer(0)] else NULL
    ## consider each segment in turn
    YY    <- as.data.frame(as.psp(Y))
    for(i in 1:nseg) {
      ## divide segment into pieces of length eps
      ## with shorter bits at each end
      leni <- len[i]
      nwhole <- floor(leni/eps)
      if(leni/eps - nwhole < 0.5 && nwhole > 2)
        nwhole <- nwhole - 1
      epsfrac <- eps/leni
      rumpfrac <- (1 - nwhole * epsfrac)/2
      brks <- c(0, rumpfrac + (0:nwhole) * epsfrac, 1)
      nbrks <- length(brks)
      ## create a dummy point at the middle of each piece
      tdumi <- (brks[-1L] + brks[-nbrks])/2
      ndumi <- length(tdumi)
      xdumi <- with(YY, x0[i] + tdumi * (x1[i]-x0[i]))
      ydumi <- with(YY, y0[i] + tdumi * (y1[i]-y0[i]))
      sdumi <- rep(i, ndumi)
      IDdumi <- 1:ndumi
      ## data points on this segment
      relevant <- (mapXY == i)
      jdati   <- which(relevant)
      xdati   <- xX[relevant]
      ydati   <- yX[relevant]
      tdati   <- tp[relevant]
      ndati   <- length(tdati)
      sdati   <- rep(i, ndati)
      IDdati  <- findInterval(tdati, brks,
                              rightmost.closed=TRUE, all.inside=TRUE)
      ## determine weights
      w <- countingweights(id=c(IDdumi, IDdati), areas=diff(brks) * leni)
      wdumi <- w[1:ndumi]
      wdati <- w[-(1:ndumi)]
      ##
      if(ismulti) {
        ## attach correct marks to data points
        mdati <- marx[relevant]
        ## replicate dummy points with each mark
        ## also add points at data locations with other marks
        xdumiU <- xdumi
        ydumiU <- ydumi
        sdumiU <- sdumi
        tdumiU <- tdumi
        wdumiU <- wdumi
        for(k in seq_len(length(flev))) {
          le <- flev[k]
          avoid <- (mdati != le)
          mdumi <- c(mdumi, rep(le, ndumi), rep(le, sum(avoid)))
          xdumi <- c(xdumi, xdumiU,          xdati[avoid])
          ydumi <- c(ydumi, ydumiU,          ydati[avoid])
          sdumi <- c(sdumi, sdumiU,          sdati[avoid])
          tdumi <- c(tdumi, tdumiU,          tdati[avoid])
          wdumi <- c(wdumi, wdumiU,          wdati[avoid])
        }
      }
      ## concatenate data points
      xdat <- c(xdat, xdati)
      ydat <- c(ydat, ydati)
      sdat <- c(sdat, sdati)
      tdat <- c(tdat, tdati)
      wdat <- c(wdat, wdati)
      jdat <- c(jdat, jdati)
      ## concatenate dummy points
      xdum <- c(xdum, xdumi)
      ydum <- c(ydum, ydumi)
      sdum <- c(sdum, sdumi)
      tdum <- c(tdum, tdumi)
      wdum <- c(wdum, wdumi)
      if(ismulti) {
        mdat <- c(mdat, mdati)
        mdum <- c(mdum, mdumi)
      }
    }
    ## map data points back to their original order
    odat <- order(jdat)
    xdat <- xdat[odat]
    ydat <- ydat[odat]
    wdat <- wdat[odat]
    sdat <- sdat[odat]
    tdat <- tdat[odat]
    ## Now create point patterns
    dat <- ppp(xdat, ydat, marks=mdat, window=W, check=FALSE)
    dum <- ppp(xdum, ydum, marks=mdum, window=W, check=FALSE)
  }
  ## save parameters
  dmethod <- paste("Equally spaced along each segment at spacing eps =",
                    signif(eps, 4),
                    summary(unitname(X))$plural)
  if(!epsgiven)
    dmethod <- paste0(dmethod, "\nOriginal parameter nd = ", nd)
  wmethod <- "Counting weights based on segment length"
  param <- list(dummy = list(method=dmethod),
                weight = list(method=wmethod))
  ## make quad scheme
  Qout <- quad(dat, dum, c(wdat, wdum), param=param)
  ## silently attach lines/network and local coordinates
  attr(Qout, "domain") <- Y
  attr(Qout, "localcoords") <- list(seg = c(sdat, sdum),
                                    tp  = c(tdat, tdum))
  return(Qout)
}

