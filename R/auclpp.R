##
##
##  auclpp.R
##
##  Calculate Area Under ROC curve
##       (in spatstat.linnet)
##
##
## Copyright (c) 2017-2025 Adrian Baddeley/Ege Rubak/Rolf Turner
##



auc.lpp <- function(X, covariate, ...,
           high=TRUE,  subset=NULL) {
  verifyclass(X, "lpp")
  if(needROC(...)) {
    ro <- roc(X, covariate, ..., high=high, subset=subset)
    result <- auc(ro)
  } else {
    nullmodel <- exactlppm(X) 
    result <- aucData(covariate, nullmodel, ..., high=high, subset=subset)
  }
  return(result)
}



auc.lppm <-
  function(X, ..., subset=NULL) {
    model <- X
    use.roc <- is.multitype(model) || needROC(...)
    if(use.roc) {
      ro <- roc(model, ..., subset=subset)
      aobs <- with(ro, mean(.y))
      atheo <- with(ro, mean(theo))
    } else if(is.ppm(model) && is.stationary(model)) {
      aobs <- atheo <- 1/2
    } else {
      lambda <- predict(model, ..., type="trend")
      X <- if(is.ppm(model)) data.ppm(model) else model$X
      if(!is.null(subset)) {
        #' restrict to subset
        lambda <- lambda[subset, drop=FALSE]
        X <- X[subset]
      }
      lamX <- lambda[X]
      lamW <- lambda[]
      Fl <- ecdf(lamW)
      aobs <- mean(Fl(lamX))
      atheo <- mean(lamW * Fl(lamW))/mean(lamW)
    }
    result <- c(aobs, atheo)
    names(result) <- c("obs", "theo")
    return(result)
  }




