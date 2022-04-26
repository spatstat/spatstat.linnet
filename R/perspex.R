#'  perspex.R
#'
#'  Perspective plots for linim and linfun objects
#'
#'  Copyright (C) Adrian Baddeley 2015
#'  GPL Public Licence >= 2.0


persp.linfun <- function(x, ..., main, eps=NULL, dimyx=NULL, xy=NULL) {
  if(missing(main)) main <- short.deparse(substitute(x))
  #' convert linfun to linim
  y <- as.linim(x, L=as.linnet(x), eps=eps, dimyx=dimyx, xy=xy)
  persp(y, ..., main=main)
}

persp.linim <- local({
  
  persp.linim <- function(x, ..., main, grid=TRUE, ngrid=10,
                          col.grid="grey", col.base="white",
                          neg.args=list(), warncross=FALSE) {
    xname <- short.deparse(substitute(x))
    if(missing(main)) main <- xname
    dotargs <- list(...)
    #' 
    L <- as.linnet(x)
    R <- Frame(L)
    zlim <- range(x, 0)
    #' set up perspective transformation and plot horizontal plane
    Z <- as.im(0, W=R, dimyx=ngrid)
    col.grid.used <- if(grid && (zlim[1] >= 0)) col.grid else NA
    argh <- resolve.defaults(list(x=Z, main=main,
                                  border=col.grid.used,
                                  col=col.base),
                             dotargs,
                             list(axes=FALSE, box=FALSE,
                                  zlim=zlim, zlab=xname, 
                                  scale=TRUE, expand=0.1))
    M <- do.call.matched(persp.im, argh,
                         funargs=graphicsPars("persp"))
    #' compute the projection of the linear network
    S <- as.psp(L)
    E <- S$ends
    x0y0 <- trans3dz(E[,1], E[,2], rep(0, nrow(E)), M)
    x1y1 <- trans3dz(E[,3], E[,4], rep(0, nrow(E)), M)
    #' order the segments by their depth in the picture (back to front)
    segmentsequence <- NULL
    ## new algorithm
    po <- depthrelations(x0y0$x, x0y0$z, x1y1$x, x1y1$z, verbose=warncross)
    if(!is.null(po))
      segmentsequence <- partial2total(po, nrow(E))
    ## backup
    if(is.null(segmentsequence)) {
      ## old algorithm
      depf <- pmin(x0y0$z, x1y1$z)
      segmentsequence <- order(depf)
    }
    #' extract function data
    df <- attr(x, "df")
    #' handle negative values separately if a grid is shown
    if(grid && zlim[1] < 0) {
      dfneg <- df
      dfneg$values <- pmin(0, df$values)
      df$values    <- pmax(0, df$values)
      #' plot negative part
      neg.args <- resolve.defaults(neg.args, dotargs)
      spectiveWalls(dfneg, segmentsequence, E, M, neg.args)
      #' plot baseline grid on top again
      spectiveGrid(R, ngrid, M, col=col.grid)
    }
    #' plot network
    do.call.matched(segments,
                    list(x0=x0y0$x,
                         y0=x0y0$y,
                         x1=x1y1$x,
                         y1=x1y1$y, ...))
    #' plot function above grid (or entire function if no grid)
    spectiveWalls(df, segmentsequence, E, M, dotargs)
    #'   
    invisible(M)
  }

  trans3dz <- function(x,y,z,pmat) {
    tr <- cbind(x, y, z, 1) %*% pmat
    list(x = tr[, 1]/tr[, 4],
         y = tr[, 2]/tr[, 4],
         z = tr[, 3]/tr[, 4])
  }

  spectiveWalls <- function(df, segmentsequence, E, M, pargs) {
    #' split by segment
    for(i in segmentsequence) {
      dfi <- df[df$mapXY == i, , drop=FALSE]
      if((nn <- nrow(dfi)) > 0) {
        #' order by position along segment
        ord <- order(dfi$tp)
        dfi <- dfi[ord, , drop=FALSE]
        #' extrapolate to segment endpoints
        Ei <- E[i,, drop=FALSE]
        xx <- c(Ei$x0, dfi$x, Ei$x1)
        yy <- c(Ei$y0, dfi$y, Ei$y1)
        zz <- with(dfi, c(values[1], values, values[nn]))
        #' make polygon in 3D      
        xx <- c(xx, rev(xx))
        yy <- c(yy, rev(yy))
        zz <- c(zz, rep(0, nn+2))
        #' project polygon
        xy <- trans3d(xx, yy, zz, M)
        do.call.matched(polygon,
                        resolve.defaults(list(x=xy$x, y=xy$y),
                                         pargs,
                                         list(col="grey")))
      }
    }
    invisible(NULL)
  }

  spectiveGrid <- function(B, ngrid, M, ...) {
    B <- Frame(B)
    xr <- B$xrange
    yr <- B$yrange
    xx <- seq(xr[1], xr[2], length.out=ngrid+1)
    yy <- seq(yr[1], yr[2], length.out=ngrid+1)
    spectiveSegments(xr[1], yy, xr[2], yy, ..., M=M)
    spectiveSegments(xx, yr[1], xx, yr[2], ..., M=M)
    invisible(NULL)
  }

  spectiveSegments <- function(x0, y0, x1, y1, ..., M) {
    a0 <- trans3dz(x0, y0, 0, M)
    a1 <- trans3dz(x1, y1, 0, M)
    do.call.matched(segments,
                    list(x0=a0$x,
                         y0=a0$y,
                         x1=a1$x,
                         y1=a1$y, ...))
    invisible(NULL)
  }

    
  depthrelations <- function(x0, z0, x1, z1, verbose=TRUE, useC=TRUE) {
    ## Find which segments lie in front of/behind other segments
    ## x is the projected x-coordinate
    ## z is the negative depth (z increases as we move closer to the viewer)
    front <- back <- integer(0)
    n <- length(x0)
    if(n >= 2) {
      ## enforce x0 <= x1
      if(any(swap <- (x0 > x1))) {
        tmp <- x1[swap]
        x1[swap] <- x0[swap]
        x0[swap] <- tmp
        tmp <- z1[swap]
        z1[swap] <- z0[swap]
        z0[swap] <- tmp
      }
      if(useC) {
        stuff <- .Call(SL_depthrel,
                       x0, z0, x1, z1, Verb=verbose,
                       PACKAGE="spatstat.linnet")
        front <- stuff[[1]]
        back  <- stuff[[2]]
        status <- stuff[[3]]
        if(status != 0) return(NULL)
        return(data.frame(front=front, back=back))
      }
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          if(x1[i] > x0[j] && x1[j] > x0[i]) {
            ## overlap occurs
            ## consider left side
            z0i <- z0[i]
            z0j <- z0[j]
            if(x0[j] < x0[i]) {
              xleft <- x0[i]
              tval <- (xleft - x0[j])/(x1[j]-x0[j])
              if(!is.finite(tval)) tval <- 0
              z0j <- z0j + tval * (z1[j] - z0j)
            } else {
              xleft <- x0[j]
              tval <- (xleft - x0[i])/(x1[i]-x0[i])
              if(!is.finite(tval)) tval <- 0
              z0i <- z0i + tval * (z1[i] - z0i)
            }
            ## consider right side
            z1i <- z1[i]
            z1j <- z1[j]
            if(x1[j] > x1[i]) {
              xright <- x1[i]
              tval <- (xright - xleft)/(x1[j]-xleft)
              if(!is.finite(tval)) tval <- 0
              z1j <- z0j + tval * (z1j - z0j)
            } else {
              xright <- x1[j]
              tval <- (xright - xleft)/(x1[i]-xleft)
              if(!is.finite(tval)) tval <- 0
              z1i <- z0i + tval * (z1i - z0i)
            }
            ## now determine which is in front
            if(z0i >= z0j && z1i >= z1j) {
              ## 'i' is in front
              front <- c(front, i)
              back  <- c(back, j)
            } else if(z0i <= z0j && z1i <= z1j){
              ## 'j' is in front
              front <- c(front, j)
              back  <- c(back, i)
            } else {
              if(verbose) {
                warning(paste("segments", i, "and", j, "cross:",
                              "z0i=", z0i, "z0j=", z0j,
                              "z1i=", z1i, "z1j=", z1j))
              }
              return(NULL)
            }
          }
        }
      }
    }
    return(data.frame(front=front, back=back))
  }

  partial2total <- function(po, n=max(po), verbose=TRUE) {
    ##  'po' represents a partial order amongst 'n' items.
    ##   Find an ordering of the 'n' items that respects the partial order.
    front <- po$front
    back  <- po$back
    unresolved <- rep(TRUE, n)
    sequence <- integer(0)
    front <- factor(front, levels=1:n)
    while(any(unresolved)) {
      bad <- sapply(split(unresolved[back], front), 
                    any)
      found <- unresolved & !bad
      if(!any(found)) {
        if(verbose)
          warning("Could not resolve some entries in the partial order")
        return(NULL)
      }
      sequence <- c(sequence, which(found))
      unresolved[found] <- FALSE
    }
    return(sequence)
  }

  persp.linim
})
