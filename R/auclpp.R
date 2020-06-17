##
## auclpp.R
##
##  Calculate ROC curve or area under it
##
##  code for linear networks
##
## $Revision: 1.1 $ $Date: 2020/06/17 04:32:13 $

roc.lpp <- function(X, covariate, ..., high=TRUE) {
  nullmodel <- lppm(X)
  result <- rocData(covariate, nullmodel, ..., high=high)
  return(result)
}

roc.lppm <- function(X, ...) {
  stopifnot(is.lppm(X))
  model <- X
  lambda <- predict(model, ...)
  Y <- X$X
  nullmodel <- lppm(Y)
  result <- rocModel(lambda, nullmodel, ...)
  return(result)
}

#    ......................................................

auc.lpp <- function(X, covariate, ..., high=TRUE) {
  d <- spatialCDFframe(lppm(X), covariate, ...)
  U <- d$values$U
  EU <- mean(U)
  result <- if(high) EU else (1 - EU) 
  return(result)
}

auc.lppm <- function(X, ...) {
  stopifnot(inherits(X, "lppm"))
  model <- X
  if(is.multitype(model)) {
    # cheat
    ro <- roc(model, ...)
    aobs <- with(ro, mean(fobs))
    atheo <- with(ro, mean(ftheo))
  } else {
    lambda <- predict(model, ...)
    Fl <- ecdf(lambda[])
    lamX <- lambda[model$X]
    aobs <- mean(Fl(lamX))
    atheo <- mean(lambda[] * Fl(lambda[]))/mean(lambda)
  }
  result <- c(aobs, atheo)
  names(result) <- c("obs", "theo")
  return(result)
}


