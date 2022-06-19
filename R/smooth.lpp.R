#'
#' smooth.lpp.R
#'
#' $Revision: 1.3 $ $Date: 2022/06/19 05:21:40 $

Smooth.lpp <- function(X, sigma,
                       ...,
                       at=c("pixels", "points"),
                       weights=rep(1, npoints(X)),
                       leaveoneout=TRUE) {
  marx <- marks(X)
  at <- match.arg(at)
  if(missing(sigma)) stop("Bandwidth sigma is required")
  check.1.real(sigma)
  if(!is.null(dim(marx)) && ncol(marx) > 1) {
    y <- solapply(as.data.frame(marx), setmarks, x=X)
    z <- lapply(y, Smooth.lpp, sigma=sigma, at=at, weights=weights, ...)
    z <- switch(at,
                pixels = as.solist(z),
                points = simplify2array(z))
    return(z)
  }
  ## Single column of marks
  marx <- as.numeric(marx)
  U <- unmark(X)
  ## Calculate
  Numer <- density(U, sigma=sigma, ..., at=at, leaveoneout=leaveoneout,
                   weights=weights * marx)
  Denom <- density(U, sigma=sigma, ..., at=at, leaveoneout=leaveoneout,
                   weights=weights)
  z <- Numer/Denom
  if(at == "points" && any(bad <- !is.finite(z))) {
    ## underflow
    z[bad] <- marx[nnwhich(U)[bad]]
  }
  attr(z, "sigma") <- sigma
  return(z)
}

