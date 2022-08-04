#'
#'   rThomaslpp.R
#'
#' Thomas process on linear network, using heat kernel
#'
#' e.g. rThomaslpp(4, 5, 0.1, simplenet)
#' 
#'   $Revision: 1.4 $  $Date: 2022/08/04 05:30:48 $

rThomaslpp <- function(kappa, scale, mu, L, ..., nsim=1, drop=TRUE) {
  if(missing(L) || is.null(L)) {
    L <- as.linnet(kappa)
  } else stopifnot(is.linnet(L))
  check.1.integer(nsim)
  simlist <- vector(mode="list", length=nsim)
  parentlist <- rpoislpp(kappa, L, ..., nsim=nsim, drop=FALSE)
  for(i in seq_len(nsim)) {
    X <- parentlist[[i]]
    Y <- density(X, scale)
    Y <- eval.linim(Y * mu)
    Z <- rpoislpp(Y, L)
    attr(Z, "parents") <- X
    attr(Z, "Lambda") <- Y
    simlist[[i]] <- Z
  }
  result <- simulationresult(simlist, nsim, drop)
  return(result)
}
