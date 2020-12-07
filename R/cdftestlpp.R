#'
#'  cdftestlpp.R
#'
#'  methods for linear networks
#' 
#'  $Revision: 1.1 $  $Date: 2020/06/17 04:34:21 $

cdf.test.lpp <-
  function(X, covariate, test=c("ks", "cvm", "ad"), ...,
           interpolate=TRUE, jitter=TRUE) {
    Xname <- short.deparse(substitute(X))
    covname <- singlestring(short.deparse(substitute(covariate)))
    test <- match.arg(test)
    if(is.character(covariate)) covname <- covariate
    if(!is.marked(X, dfok=TRUE)) {
      # unmarked
      model <- lppm(X)
      modelname <- "CSR"
    } else if(is.multitype(X)) {
      # multitype
      mf <- table(marks(X))
      if(all(mf > 0)) {
        model <- lppm(X ~ marks)
        modelname <- "CSRI"
      } else {
        warning("Ignoring marks, because some mark values have zero frequency")
        X <- unmark(X)
        model <- lppm(X)
        modelname <- "CSR"
      } 
    } else {
      # marked - general case
      X <- unmark(X)
      warning("marks ignored")
      model <- lppm(X)
      modelname <- "CSR"
    }
    do.call(spatialCDFtest,
            resolve.defaults(list(quote(model), 
				  quote(covariate), 
				  test=test),
                             list(interpolate=interpolate, jitter=jitter),
                             list(...),
                             list(modelname=modelname,
                                  covname=covname, dataname=Xname)))
}

cdf.test.lppm <- function(model, covariate,
                          test=c("ks", "cvm", "ad"),
                          ..., interpolate=TRUE, jitter=TRUE,
                          nsim=99, verbose=TRUE) {
  modelname <- short.deparse(substitute(model))
  covname <- singlestring(short.deparse(substitute(covariate)))
  test <- match.arg(test)
  verifyclass(model, "lppm")
  if(is.character(covariate)) covname <- covariate
  if(is.poisson(model) && is.stationary(model))
    modelname <- "CSR"
  do.call(spatialCDFtest,
          resolve.defaults(list(quote(model), 
				quote(covariate), 
				test=test),
                           list(interpolate=interpolate, jitter=jitter,
                                nsim=nsim, verbose=verbose),
                           list(...),
                           list(modelname=modelname,
                                covname=covname)))
}

