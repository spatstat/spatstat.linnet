#'
#' unstacklpp.R
#'
#' $Revision: 1.1 $ $Date: 2020/06/17 04:38:38 $

unstack.lpp <- function(x, ...) {
  trap.extra.arguments(...)
  marx <- marks(x)
  d <- dim(marx)
  if(is.null(d)) return(solist(x))
  y <- rep(list(unmark(x)), d[2])
  for(j in seq_along(y))
    marks(y[[j]]) <- marx[,j,drop=FALSE]
  names(y) <- colnames(marx)
  return(as.solist(y))
}

