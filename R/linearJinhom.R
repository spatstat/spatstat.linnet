#'
#' linearJinhom.R
#'
#' Author: Mehdi Moradi
#' Rewritten by: Adrian Baddeley
#' 
#' Copyright (c) Mehdi Moradi and Adrian Baddeley 2020-2022
#' GNU Public Licence >= 2.0
#' 
#' $Revision: 1.14 $ $Date: 2022/06/16 06:12:10 $
#'

linearJinhom <- function(X, lambda=NULL, lmin=NULL,
                         ...,
                         r=NULL, rmax=NULL,
                         distance=c("path","euclidean"),
                         densitymethod=c("kernel", "Voronoi"),
                         sigma=bw.scott.iso,
                         f=0.2, nrep=200,
                         ngrid=256) {
  verifyclass(X, "lpp")
  nX <- npoints(X)
  L <- domain(X)
  distance <- match.arg(distance)
  
  ## sample the network L on a grid
  gridDF <- pointsAlongNetwork(L, volume(L)/ngrid)
  nG <- nrow(gridDF)
  G <- lpp(gridDF, L)

  ## 'border' 
  B <- terminalvertices(L)

  ## Determine intensity values
  if(is.null(lambda)){
    ## Estimate intensity 
    densitymethod <- match.arg(densitymethod)
    switch(densitymethod,
           kernel = {
             if(is.function(sigma)) sigma <- sigma(X)
             lambda <- density(X, sigma = sigma, distance = "euclidean",...)
           },
           Voronoi = {
             lambda <- densityVoronoi.lpp(X,f=f,nrep = nrep,...)
           })
    ## Evaluate intensity at data/grid points
    lambdaX <- lambda[X, drop=FALSE]  # intensity at data points
    if(anyNA(lambdaX)) {
      bad <- is.na(lambdaX)
      lambdaX[bad] <- safelookup(as.im(lambda), as.ppp(X[bad]))
    }
  } else {
    ## Validate intensity and evaluate at data points
    lambdaX <- getlambda.lpp(lambda, X, ...,
                             update=TRUE, leaveoneout=TRUE,
                             loo.given=FALSE, lambdaname="lambda")
  }

  ## lower bound on lambda
  if(missing(lmin)) {
    ## lmin <- min(min(lambdaG), min(lambdaX))
    lmin <- min(lambdaX)
  } else {
    check.1.real(lmin)
    if(lmin > min(lambdaX)) 
      stop("lmin exceeds the minimum value of lambda")
  }
  ## measure distances
  switch(distance,
         path = {
           dXX <- pairdist(X)
           dGX <- crossdist(G, X)
           nndXB <- nncross(X, B, what="dist")
           nndGB <- nncross(G, B, what="dist")
         },
         euclidean = {
           PX <- as.ppp(X)
           PG <- as.ppp(G)
           PB <- as.ppp(B)
           dXX <- pairdist(PX)
           dGX <- crossdist(PG, PX)
           nndXB <- nncross(PX, PB, what="dist")
           nndGB <- nncross(PG, PB, what="dist")
         })
  diag(dXX) <- Inf

  ## eliminate grid points which coincide with data points
  nndGX <- apply(dGX, 1, min)
  ok <- (nndGX >= .Machine$double.eps)
  if(!all(ok)) {
    G <- G[ok]
    if(distance == "euclidean") PG <- PG[ok]
    dGX <- dGX[ok,,drop=FALSE]
    nndGB <- nndGB[ok, drop=FALSE]
    nG <- sum(ok)
  }
  
  ## determine sequence of r values
  rmaxdefault <- rmax %orifnull% quantile(nndXB, 0.6)
  brks <- handle.r.b.args(r, NULL, Window(L), rmaxdefault=rmaxdefault)
  r <- brks$r
  rmax <- brks$max
  nr <- length(r)

  ## calculate geometrical edge weights
  switch(distance,
         path = {
           toler <- default.linnet.tolerance(L)
           ## count number m(u,r) of endpoints of disc of centre u and radius r
           ## mXX[i,j] = m(X[i], dXX[i,j])
           mXX <- DoCountEnds(X, dXX, toler)
           ## mGX[i,j] = m(G[i], dXG[j,i])
           U <- superimpose(G, X)
           I <- (seq_len(npoints(U)) <= nG)
           J <- !I
           mGX <- DoCountCrossEnds(U, I, J, dGX, toler)
           ## weights are reciprocals of counts
           wXX <- 1/mXX
           wGX <- 1/mGX
         },
         euclidean = {
           ## intersect circles with segments
           le <- as.data.frame(as.psp(L))
           ## X to X
           stuffXX <- xysegMcircle(PX$x, PX$y, dXX, le$x0, le$y0, le$x1, le$y1)
           ## stuffXX = list(i, j, k, sinalpha)
           ## indices i, j, k specify provenance of each intersection point
           ## i = centre, j = segment, k = radius (column of dXX)
           seqX <- seq_len(nX)
           ii <- factor(stuffXX$i, levels=seqX)
           jj <- factor(stuffXX$k, levels=seqX) # stet 
           wXX <- tapply(stuffXX$sinalpha, list(ii, jj), harmonicsum)
           ## G to X
           stuffGX <- xysegMcircle(PG$x, PG$y, dGX, le$x0, le$y0, le$x1, le$y1)
           seqG <- seq_len(nG)
           ii <- factor(stuffGX$i, levels=seqG)
           jj <- factor(stuffGX$k, levels=seqX)
           wGX <- tapply(stuffGX$sinalpha, list(ii, jj), harmonicsum)
         })
  
  wXX[!is.finite(wXX)] <- 0
  wGX[!is.finite(wGX)] <- 0

  ######## Calculate estimates ##################################

  RlambdaX <- pmax(lambdaX/lmin, 1)

  ## nearest neighbour calculation (survival function 1 - G(r))
  Gsurv <- rep(1,nr)
  for(i in 1:nr) {
    uncensored <- (nndXB > r[i])
    if(any(uncensored)) {
      counted <- (dXX <= r[i])
      diag(counted) <- FALSE
      if(any(counted)) {
        tot <- 0
        for (j in which(uncensored)) {
          p <- which(counted[j,])
          if(length(p)) {
            tot <- tot+prod(1 - wXX[j,p]/RlambdaX[p])
          } else {
            tot <- tot+1
          }
        }
        Gsurv[i] <- tot/(sum(uncensored))
      }
    }
  }

  ## empty space calculation (survival function 1-F(r))
  Fsurv <- rep(1, nr)
  for(i in 1:nr) {
    uncensored <- (nndGB > r[i])
    if(any(uncensored)) {
      counted <- (dGX <= r[i])
      if(any(counted)) {
        tot <- 0
        for (j in which(uncensored)) {
          p <- which(counted[j,])
          if(length(p)) {
            tot <- tot+prod(1 - wGX[j,p]/RlambdaX[p])
          } else {
            tot <- tot+1
          }
        }
        Fsurv[i] <- tot/(sum(uncensored))
      }
    }
  }
  ## Finally, calculate the estimate of J
  Jvalues <- ratiotweak(Gsurv, Fsurv)

  ## .......... pack up ..............................
  ## Make fv objects
  fname <- c("J", "list(L, inhom)")
  desc <- c("distance argument r",
            "theoretical Poisson %s",
            "estimate of %s")
  Z <- fv(data.frame(r=r, theo=1, est=Jvalues),
          "r",
          quote(J[L, inhom](r)),
          "est",
          . ~ r,
          c(0,rmax),
          c("r",
            makefvlabel(NULL, NULL, fname, "pois"),
            makefvlabel(NULL, "hat", fname, "est")),
          desc,
          fname=fname,
          yexp=quote(J[list(L, "inhom")](r)))
  fname[1] <- "G"
  attr(Z, "Ginhom") <- fv(data.frame(r=r, theo=1, est=1-Gsurv),
                          "r",
                          quote(G[L, inhom](r)),
                          "est",
                          . ~ r,
                          c(0,rmax),
                          c("r",
                            makefvlabel(NULL, NULL, fname, "pois"),
                            makefvlabel(NULL, "hat", fname, "est")),
                          desc,
                          fname=fname,
                          yexp=quote(G[list(L, "inhom")](r)))
  fname[1] <- "F"
  attr(Z, "Finhom") <- fv(data.frame(r=r, theo=1, est=1-Fsurv),
                          "r",
                          quote(F[L, inhom](r)),
                          "est",
                          . ~ r,
                          c(0,rmax),
                          c("r",
                            makefvlabel(NULL, NULL, fname, "pois"),
                            makefvlabel(NULL, "hat", fname, "est")),
                          desc,
                          fname=fname,
                          yexp=quote(F[list(L, "inhom")](r)))
  attr(Z, "alim") <- range(c(0, r[Fsurv >= 0.1]))
  return(Z)
}

