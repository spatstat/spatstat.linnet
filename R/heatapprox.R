#'
#'     heatapprox.R
#' 
#'  Approximation to the heat kernel kappa(u,u) on a network,
#'  using only paths on the current segment.
#'
#'  Copyright (c) Greg McSwiggan and Adrian Baddeley 2017-2020
#'
#'  $Revision: 1.5 $  $Date: 2022/05/21 09:52:11 $
#'

heatkernelapprox <- function(X, sigma, nmax=20, floored=TRUE) {
  stopifnot(is.lpp(X))
  nX <- npoints(X)
  if(nX == 0) return(numeric(0))
  check.nvector(sigma, nX, oneok=TRUE, vname="sigma")
  stopifnot(all(sigma > 0))
  if(length(sigma) == 1) sigma <- rep(sigma, nX)
  check.1.integer(nmax)
  lenf <- lengths_psp(as.psp(domain(X)))
  coo <- coords(X)
  seg <- coo$seg
  len <- lenf[seg]
  pos <- len * coo$tp
  L <- domain(X)
  vv <- vertexdegree(L)
  dleft <- vv[L$from[seg]]
  dright <- vv[L$to[seg]]
  z <- .C(SL_heatApprox,
          n = as.integer(nX),
          a = as.double(len),
          x = as.double(pos),
          y = as.double(pos), #sic
          s = as.double(sigma),
          degl = as.integer(dleft),
          degr = as.integer(dright),
          m = as.integer(nmax),
          z = as.double(numeric(nX)),
          PACKAGE="spatstat.linnet")
  ans <- z$z
  if(floored) ans <- pmax(ans, 1/volume(L))
  return(ans)
}

                             

