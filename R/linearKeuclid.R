#'
#'    linearKeuclid.R
#'
#'    Copyright (C) Suman Rakshit, Gopalan Nair and Adrian Baddeley 2017-2022
#'
#'    GNU Public Licence 2.0
#'
#'    $Revision: 1.8 $ $Date: 2022/06/27 02:16:56 $

linearKEuclid <- function(X, r=NULL, ...) {
  stopifnot(inherits(X, "lpp"))
  # extract info about pattern
  np <- as.double(npoints(X))
  lengthL <- volume(domain(X))
  # compute K
  samplesize <- npairs <- np * (np - 1)
  denom <- npairs/lengthL
  K <- linearEuclidEngine(X, "K", r=r, denom=denom, ..., samplesize=samplesize)
  # set appropriate y axis label
  ylab <- quote(K[euclid](r))
  fname <- c("K", "euclid")
  K <- rebadge.fv(K, new.ylab=ylab, new.fname=fname)
  return(K)
}

linearpcfEuclid <- function(X, r=NULL, ...) {
  stopifnot(inherits(X, "lpp"))
  # extract info about pattern
  np <- as.double(npoints(X))
  lengthL <- volume(domain(X))
  # compute g
  samplesize <- npairs <- np * (np - 1)
  denom <- npairs/lengthL
  g <- linearEuclidEngine(X, "g", r=r, denom=denom, ..., samplesize=samplesize)
  # extract bandwidth
  bw <- attr(g, "bw")
  # set appropriate y axis label
  ylab <- quote(g[euclid](r))
  fname <- c("g", "euclid")
  g <- rebadge.fv(g, new.ylab=ylab, new.fname=fname)
  # reattach bandwidth
  attr(g, "bw") <- bw
  return(g)
}

linearKEuclidInhom <- function(X, lambda=NULL, r=NULL,  ...,
                          normalise=TRUE, normpower=2,
			  update=TRUE, leaveoneout=TRUE) {
  stopifnot(inherits(X, "lpp"))
  if(is.null(lambda))
    linearKEuclid(X, r=r, ...)
  if(normalise) {
    check.1.real(normpower)
    stopifnot(normpower >= 1)
  }
  # extract info about pattern
  lengthL <- volume(domain(X))
  #
  lambdaX <- resolve.lambda.lpp(X, lambda, ...,
               update=update, leaveoneout=leaveoneout)
  #
  invlam <- 1/lambdaX
  invlam2 <- outer(invlam, invlam, "*")
  denom <- if(!normalise) lengthL else
           if(normpower == 1) sum(invlam) else
           lengthL * (sum(invlam)/lengthL)^normpower
  K <- linearEuclidEngine(X, "K", ...,
                          r=r, reweight=invlam2, denom=denom)
  # set appropriate y axis label
  ylab <- quote(K[euclid, inhom](r))
  fname <- c("K", "list(euclid, inhom)")
  K <- rebadge.fv(K, new.fname=fname, new.ylab=ylab)
  attr(K, "dangerous") <- attr(lambdaX, "dangerous")
  return(K)
}

linearpcfEuclidInhom <- function(X, lambda=NULL, r=NULL,  ...,
                          normalise=TRUE, normpower=2,
			  update=TRUE, leaveoneout=TRUE) {
  stopifnot(inherits(X, "lpp"))
  if(is.null(lambda))
    linearpcfEuclid(X, r=r, ...)
  if(normalise) {
    check.1.real(normpower)
    stopifnot(normpower >= 1)
  }
  # extract info about pattern
  lengthL <- volume(domain(X))
  #
  lambdaX <- resolve.lambda.lpp(X, lambda, ...,
               update=update, leaveoneout=leaveoneout)
  #
  invlam <- 1/lambdaX
  invlam2 <- outer(invlam, invlam, "*")
  denom <- if(!normalise) lengthL else
           if(normpower == 1) sum(invlam) else
           lengthL * (sum(invlam)/lengthL)^normpower
  g <- linearEuclidEngine(X, "g", ...,
                          r=r, reweight=invlam2, denom=denom)
  # extract bandwidth
  bw <- attr(g, "bw")
  # set appropriate y axis label
  ylab <- quote(g[euclid, inhom](r))
  fname <- c("g", "list(euclid, inhom)")
  g <- rebadge.fv(g, new.fname=fname, new.ylab=ylab)
  attr(g, "dangerous") <- attr(lambdaX, "dangerous")
  # reattach bandwidth
  attr(g, "bw") <- bw
  return(g)
}

linearEuclidEngine <- function(X, fun=c("K", "g"), ...,
		               r=NULL, reweight=NULL,
                               denom=1, samplesize=NULL,
                               showworking=FALSE,
                               correction=NULL) {
  #                           'correction' is ignored
  ## compute either K-function or pair correlation function
  fun <- match.arg(fun)
  ## extract info about pattern
  np <- npoints(X)
  ## extract linear network
  L <- domain(X)
  W <- Window(L)
  # determine r values
  rmaxdefault <- 0.98 * rmaxEuclidean(L, verbose=FALSE, show=FALSE)
  breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
  ##
  fname <- c(fun, "euclid")
  ylab <- substitute(funky[euclid](r), list(funky=as.name(fun)))
  ##
  if(np < 2) {
    ## no pairs to count: return zero function
    zeroes <- numeric(length(r))
    df <- data.frame(r = r, est = zeroes)
    result <- fv(df, "r", ylab,
                "est", . ~ r, c(0, rmax),
                c("r", makefvlabel(NULL, "hat", fname)), 
                c("distance argument r", "estimated %s"),
                fname = fname)
    return(result)
  }
  #' compute pairwise distances
  Y <- as.ppp(X)
  D <- pairdist(Y)
  diag(D) <- 2 * diameter(Frame(X))
  #' intersect circles with segments
  df <- as.data.frame(as.psp(L))
  stuff <- xysegMcircle(Y$x, Y$y, D, df$x0, df$y0, df$x1, df$y1)
  #' compute edge weights
  seqX <- seq_len(np)
  ii <- factor(stuff$i, levels=seqX)
  jj <- factor(stuff$k, levels=seqX) # stet 
  edgewt <- tapply(stuff$sinalpha, list(ii, jj), harmonicsum)
  edgewt[is.na(edgewt)] <- 0
  ## compute K or g
  wt <- if(!is.null(reweight)) edgewt * reweight else edgewt
  switch(fun,
         K = {
	   result <- compileK(D, r, weights=wt, denom=denom, fname=fname,
                              samplesize=samplesize)
	   bw <- NA
	   theo <- r
	 },
	 g = {
	   result <- compilepcf(D, r, weights=wt, denom=denom, fname=fname,
	                        ..., samplesize=samplesize)
	   bw <- attr(result, "bw")			
           theo <- rep.int(1, length(r))
	 })
  ## tack on theoretical value
  result <- bind.fv(result, data.frame(theo=theo),
                    makefvlabel(NULL, NULL, fname, "theo"),
                   "theoretical Poisson %s")
  result <- rebadge.fv(result, ylab, fname)
  unitname(result) <- unitname(X)
  fvnames(result, ".") <- rev(fvnames(result, "."))
  ## tack on bandwidth again
  if(fun == "g") attr(result, "bw") <- bw
  ## show working
  if(showworking)
    attr(result, "working") <- list(D=D, wt=wt)
  return(result)
}

#' calculate rmax for Euclidean distance on linear network
#'

rmaxEuclidean <- function(L, verbose=TRUE, show=FALSE) {
  L <- as.linnet(L)
  V <- vertices(L)
  # pick out the convex null
  V <- V[chull(V)]
  # extract segments and their endpoints
  S <- as.psp(L)
  ns <- nsegments(S)
  A <- endpoints.psp(S, "first")
  B <- endpoints.psp(S, "second")
  len <- lengths_psp(S)
  # find furthest points
  DA <- crossdist(A, V)
  DB <- crossdist(B, V)
  MA <- apply(DA, 1, max)
  MB <- apply(DB, 1, max)
  # upper and lower bounds on minimum for each segment
  Mup <- pmin(MA, MB)  # min of a subset
  Mlow <- pmax(0, Mup - len)  # triangle inequality
  # upper and lower bounds on rmax 
  rmax1 <- min(Mup)
  rmax0 <- min(Mlow)
  # restrict attention to segments capable of achieving minimum
  ok <- (Mlow <= rmax1)
  if(verbose)
    cat(paste("Eliminated", sum(!ok), "out of", ns, "segments\n"))
  Sok <- S[ok]
  # brute force minimisation along segments
  X <- pointsOnLines(Sok, np=1000)
  D <- crossdist(X, V)
  maxDV <- apply(D, 1, max)
  rmax <- min(maxDV)
  if(verbose)
    cat(paste0("Lower bound: \t", format(rmax0), "\n",
              "Estimate: \t", format(rmax), "\n",
              "Upper bound: \t", format(rmax1), "\n"))
  if(show) {
    ii <- which.min(maxDV)
    jj <- which(D[ii,] == rmax)
    Xcentre <- X[ii]
    Vmax <- V[jj]
    plot(L, main="")
    plot(Xcentre, add=TRUE, pch=16)
    plot(Vmax, add=TRUE)
    plot(disc(centre=Xcentre, radius=rmax), add=TRUE)
    attr(rmax, "centre") <- coords(Xcentre)
  }
  return(rmax)
}

