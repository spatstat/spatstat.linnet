##
##
##  roclpp.R
##
##  Calculate ROC curve
##       (in spatstat.linnet)
##
## Copyright (c) 2017-2025 Adrian Baddeley/Ege Rubak/Rolf Turner
##


roc.lpp <-
  function(X, covariate, 
                               ...,
                               baseline=NULL,
                               high = TRUE,
                               weights = NULL,
                               method = "raw",
                               CI = "none", alpha=0.05,
                               subset=NULL) {
  callframe <- parent.frame()
  stopifnot(is.lpp(X))
  observations <- "exact" # presence-absence analysis is not implemented
  ## estimation method
  method <- match.arg(method, c("raw", "monotonic", "smooth", "rhohat", "all"),
                      several.ok = TRUE)
  method <- sub("rhohat", "smooth", method)
  if("all" %in% method)
    method <- c("raw", "monotonic", "smooth", "all")
  ## Resolve null model / baseline
  baseline.is.pattern <- is.ppp(baseline) || is.lpp(baseline)
  if(!baseline.is.pattern) {
    ## baseline is NULL, or a function/image/...  that determines null model
    nullmodel <- resolveNullModel(baseline, X, observations, ...)
  } else {
    ## baseline is a set of dummy points
    if(observations == "presence") {
      ## discretise both patterns
      X <- do.call.matched(discretise, 
                           list(X=quote(X), ..., move.points=TRUE),
                           envir=callframe)
      baseline <-
        do.call.matched(discretise,
                        list(X=quote(baseline), ..., move.points=TRUE),
                        envir=callframe)
    }
  }
  ## confidence intervals using which estimate?
  CI <- match.arg(CI, c("none", "raw", "monotonic", "smooth", "rhohat"))
  CI <- sub("rhohat", "smooth", CI)
  if(CI != "none") {
    method <- union(method, CI)
    if(baseline.is.pattern && CI != "smooth")
      warning(paste("Confidence interval for", sQuote(CI), "method",
                    "treats the baseline ('control') point pattern as fixed"),
              call.=FALSE)
  }
  #' Get covariates
  covariate <- digestCovariates(covariate, W = Window(X))
  #' ensure 'high' is a vector
  ncov <- length(covariate)
  if(length(high) == 1)
    high <- rep(high, ncov)
  stopifnot(length(high) == ncov)
  #' process
  result <- vector(mode="list", length=ncov)
  for(i in seq_len(ncov)) {
    if(baseline.is.pattern){
      result[[i]] <- rocDummy(X, baseline, covariate[[i]],
                              method=method,
                              high = high[i],
                              weights=weights,
                              subset=subset,
                              CI=CI, alpha=alpha,
                              ...)
    } else{
      result[[i]] <- rocEngine(covariate[[i]],
                               nullmodel,
                               method = method,
                               high=high[i],
                               weights=weights,
                               subset=subset,
                               CI=CI, alpha=alpha,
                               ...)
      }
  }
  names(result) <- names(covariate)
  result <- as.anylist(result)
  if(ncov == 1)
    result <- result[[1]]
  return(result)
}



  roc.lppm <-
  function(X, covariate=NULL, ...,
           baseline=NULL, high=TRUE,
           method = "raw", CI = "none", alpha=0.05,
           leaveoneout=FALSE, subset=NULL) {
  model <- X
  X <- response(model)
  traditional <- is.null(covariate)
  lambda <- do.call.matched(predict,
                            list(object=model, ...),
                            c("object", "ngrid", "dimyx", "eps",
                              "correction", "new.coef"))
  if(traditional) {
    #' discriminant is the predicted intensity
    covariate <- lambda
    covtype <- if(is.slrm(model)) "probability" else "intensity"
  } else {
    #' discriminant is a user-specified covariate
    covariate <- digestCovariates(covariate, W=Window(model))
    if(length(covariate) == 1) covariate <- covariate[[1L]]
    covtype <- "covariate"
  }
  if(is.ppp(baseline) || is.lpp(baseline)) {
    if(is.slrm(model) && is.null(model$Data$dataAtPoints)) {
      W <- model$Data$W
      X <- discretise(X, xy=W, move.points=TRUE)
      baseline <- discretise(baseline, xy=W, move.points=TRUE)
    }
    result <- rocDummy(X, baseline, covariate, ..., high=high, method=method,
                       CI=CI, alpha=alpha, subset=subset)
  } else {
    nullmodel <- resolveNullModel(baseline, model)
    leaveoneout <- if(!traditional) FALSE else as.logical(leaveoneout)
    ## leaveoneout can be TRUE, FALSE or c(TRUE, FALSE)
    if(any(leaveoneout)) {
      fittedX <- fitted(model, dataonly=TRUE, leaveoneout=TRUE)
      Rminus <- rocEngine(covariate, nullmodel, method = method,
                          covtype = covtype,
                          fittedmodel = if(traditional) NULL else model,
                          discrimAtPoints = fittedX, # NB
                          leftoneout=TRUE, # for labelling
                          high=high,
                          CI = CI, alpha=alpha,
                          subset=subset,
                          ...)
    } else {
      Rminus <- NULL
    }
    if(any(!leaveoneout)) {
      Rplus <- rocEngine(covariate, nullmodel, method = method,
                         covtype = covtype,
                         fittedmodel = if(traditional) NULL else model,
                         discrimAtPoints = NULL,  # NB
                         high=high,
                         CI = CI, alpha=alpha,
                         subset=subset,
                         ...)
    } else {
      Rplus <- NULL
    }
    if(all(leaveoneout)) {
      ## return only the leave-one-out calculation
      result <- Rminus
    } else if(all(!leaveoneout)) {
      ## return only the vanilla calculation
      result <- Rplus
    } else {
      ## return both results, combined
      common <- c("null", "theo", "fun", "thresh")
      if(leaveoneout[1]) {
        Rfirst <- Rminus
        Rsecond <- Rplus
      } else {
        Rfirst <- Rplus
        Rsecond <- Rminus
      }
      result <- bind.fv(Rfirst,
                        Rsecond[, setdiff(colnames(Rsecond), common)])
      ## preferred column
      fvnames(result, ".y") <- fvnames(Rfirst, ".y")
      ## shade columns if any
      if(!is.null(sn <- fvnames(Rfirst, ".s")))
        fvnames(result, ".s") <- sn
      ## dotnames
      dn <- fvnames(Rplus, ".")
      dnew <- character(0)
      for(d in dn) {
        if(!leaveoneout[1]) dnew <- c(dnew, d)
        if(!(d %in% common))
          dnew <- c(dnew, makeRocTag(d, TRUE))
        if(leaveoneout[1]) dnew <- c(dnew, d)
      }
      fvnames(result, ".") <- dnew
      ## return as roc
      class(result) <- union("roc", class(result))
    }
  }
  return(result)
}









