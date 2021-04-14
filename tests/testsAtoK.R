#'
#'   Header for all (concatenated) test files
#'
#'   Require spatstat.linnet
#'   Obtain environment variable controlling tests.
#'
#'   $Revision: 1.5 $ $Date: 2020/04/30 05:31:37 $

require(spatstat.linnet)
FULLTEST <- (nchar(Sys.getenv("SPATSTAT_TEST", unset="")) > 0)
ALWAYS   <- TRUE
cat(paste("--------- Executing",
          if(FULLTEST) "** ALL **" else "**RESTRICTED** subset of",
          "test code -----------\n"))
#'  tests/aucroc.R
#'
#'  AUC and ROC code
#'
#'  $Revision: 1.6 $ $Date: 2020/11/02 06:26:45 $

local({
  if(FULLTEST) {
    A <- roc(spiders, "x")
    B <- auc(spiders, "y")
    fut <- lppm(spiders ~ I(y-x))
    f <- roc(fut)
    g <- auc(fut)
  }
})
## tests/cdf.test.R


local({
  NSIM <- 9
  op <- spatstat.options(ndummy.min=16, npixel=32)
  if(ALWAYS) {
    ## (3) linear networks
    set.seed(42)
    X <- runiflpp(20, simplenet)
    cdf.test(X, "x")
    if(FULLTEST) {
      cdf.test(X, "x", "cvm")
      cdf.test(X %mark% runif(20), "x")
    }
    fit <- lppm(X ~1)
    cdf.test(fit, "y", "cvm", nsim=NSIM)
    if(FULLTEST) {
      cdf.test(fit, "y", nsim=NSIM)
      cdf.test(fit, "y", "ad", nsim=NSIM)
    }
    if(FULLTEST) {
      ## marked
      cdf.test(chicago, "y")
      cdf.test(subset(chicago, marks != "assault"), "y")
    }
  }
  reset.spatstat.options()
})


#'
#'   tests/cluck.R
#'
#'   Tests of "click*" functions
#'   using queueing feature of spatstatLocator
#'
#'   $Revision: 1.7 $ $Date: 2020/11/02 06:53:30 $

local({
  Y <- coords(runiflpp(6, simplenet))
  if(FULLTEST) {
    #' clicklpp
    spatstat.utils::queueSpatstatLocator(Y)
    XL <- clicklpp(simplenet)
  }
  if(ALWAYS) {
    spatstat.utils::queueSpatstatLocator(Y)
    XM <- clicklpp(simplenet, n=3, types=c("a", "b"))
  }
  if(ALWAYS) {
    #' lineardisc
    plot(simplenet)
    spatstat.utils::queueSpatstatLocator(as.ppp(runiflpp(1, simplenet)))
    V <- lineardisc(simplenet, r=0.3)
  }
})
#'
#'   tests/disconnected.R
#'
#'   disconnected linear networks
#'
#'    $Revision: 1.4 $ $Date: 2020/04/28 12:58:26 $


local({

  #'   disconnected network
  m <- simplenet$m
  m[4,5] <- m[5,4] <- m[6,10] <- m[10,6] <- m[4,6] <- m[6,4] <- FALSE
  L <- linnet(vertices(simplenet), m)
  if(FULLTEST) {
    L
    summary(L)
    is.connected(L)
    Z <- connected(L, what="components")
  }

  #' point pattern with no points in one connected component
  set.seed(42)
  X <- rpoislpp(lambda=function(x,y) { 10 * (x < 0.5)}, L)
  B <- lineardirichlet(X)
  if(FULLTEST) {
    plot(B)
    summary(B)
  }
  if(ALWAYS) {
    D <- pairdist(X)
    A <- nndist(X)
  }
  if(FULLTEST) {
    H <- nnwhich(X)
    Y <- rpoislpp(lambda=function(x,y) { 10 * (x < 0.5)}, L)
    G <- nncross(X, Y)
    J <- crossdist(X, Y)
    plot(distfun(X))  # includes evaluation of nncross(what="dist")
  }
  
  #' K functions in disconnected network
  if(ALWAYS) {
    K <- linearK(X)
    lamX <- intensity(X)
    nX <- npoints(X)
    KI <- linearKinhom(X, lambda=rep(lamX, nX))
    P <- linearpcf(X)
    PJ <- linearpcfinhom(X, lambda=rep(lamX, nX))
  }
  Y <- X %mark% factor(rep(1:2, nX)[1:nX])
  if(FULLTEST) {
    Y1 <- split(Y)[[1]]
    Y2 <- split(Y)[[2]]
    KY <- linearKcross(Y)
    PY <- linearpcfcross(Y)
    KYI <- linearKcross.inhom(Y, lambdaI=rep(intensity(Y1), npoints(Y1)),
                              lambdaJ=rep(intensity(Y2), npoints(Y2)))
    PYI <- linearpcfcross.inhom(Y, lambdaI=rep(intensity(Y1), npoints(Y1)),
                                lambdaJ=rep(intensity(Y2), npoints(Y2)))
  }

  #' internal utilities
  if(FULLTEST) {
    K <- ApplyConnected(X, linearK, rule=function(...) list())
  }
})


#
#  tests/envelopes.R
#
#  Test validity of envelope data
#
#  $Revision: 1.24 $  $Date: 2020/11/02 06:53:20 $
#


if(FULLTEST) {
local({
  X <- runiflpp(10, simplenet)
  Xr <- X %mark% runif(10)
  Xc <- X %mark% factor(letters[c(1:4,3,2,4:1)])
  X2 <- X %mark% data.frame(height=runif(10), width=runif(10))

  E  <- envelope(X,  linearK, nsim=9)
  Er <- envelope(Xr, linearK, nsim=9)
  Ec <- envelope(Xc, linearK, nsim=9)
  E2 <- envelope(X2, linearK, nsim=9)
  
  Erf <- envelope(Xr, linearK, nsim=9, fix.n=TRUE)
  E2f <- envelope(X2, linearK, nsim=9, fix.n=TRUE)

  Ecf <- envelope(Xc, linearK,      nsim=9, fix.n=TRUE)
  Ecm <- envelope(Xc, linearKcross, nsim=9, fix.n=TRUE, fix.marks=TRUE)

  fut <- lppm(Xc ~ marks)
  EEf <- envelope(fut, linearK,      fix.n=TRUE)
  EEm <- envelope(fut, linearKcross, fix.n=TRUE, fix.marks=TRUE)
})
}

#
#  tests/func.R
#
#   $Revision: 1.8 $   $Date: 2020/12/03 03:28:44 $
#
#  Tests of 'funxy' infrastructure etc

if(FULLTEST) {
local({
  ## Check the peculiar function-building code in funxy
  W <- square(1)
  f1a <- function(x, y) sqrt(x^2 + y^2)
  F1a <- funxy(f1a, W)
  Y <- runiflpp(5, simplenet)
  b <- F1a(Y)
})
}


#'     tests/hypotests.R
#'     Hypothesis tests
#' 
#'  $Revision: 1.9 $ $Date: 2020/11/02 06:39:23 $

if(FULLTEST) {
local({
  berman.test(spiders, "x")
  berman.test(lppm(spiders ~ x), "y")

})
}
#
#  tests/imageops.R
#
#   $Revision: 1.32 $   $Date: 2021/04/14 08:57:21 $
#


if(FULLTEST) {
  local({
    d <- distmap(cells, dimyx=32)
    ## linear networks
    ee  <- d[simplenet, drop=FALSE]
    eev <- d[simplenet]
  })
}
