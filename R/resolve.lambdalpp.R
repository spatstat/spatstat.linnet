#'
#'  resolve.lambdalpp.R
#'
#'  $Revision: 1.3 $ $Date: 2026/02/16 02:01:30 $


resolve.lambda.lpp <- function(X, lambda, subset=NULL, ...,
                               update=TRUE, leaveoneout=TRUE,
                               everywhere=FALSE,
                               loo.given=TRUE,
                               sigma=NULL,
                               lambdaname,
                               subsetlevels=NULL) {
  if(missing(lambdaname)) lambdaname <- short.deparse(substitute(lambda))
  Y <- if(is.null(subset)) X else X[subset]
  danger <- TRUE
  if(is.null(lambda)) {
    ## compute a kernel density estimate at X
    lambdaX <- density(X, sigma=sigma, ..., at="points",
                       leaveoneout=leaveoneout)
    ## ensure it is positive
    lambdaX <- pmax(lambdaX, .Machine$double.eps)
    ## Restrict if required
    lambdaY <- if(is.null(subset)) lambdaX else lambdaX[subset]
    ## Compute at all locations if required
    if(everywhere)
      everylambda <- density(X, sigma=sigma, ..., at="pixels") 
  } else if(is.ppm(lambda) || is.lppm(lambda)) {
    ## fitted model
    model <- lambda
    if(loo.given && leaveoneout && !update) 
      warning("leave-one-out calculation for fitted models is only available when update=TRUE",
              call.=FALSE)
    leaveoneout <- leaveoneout && update
    ##
    if(!is.multitype(model)) {
      ## ............ Unmarked point process model for unmarked subset
      UY <- unmark(Y)
      if(update) {
        ## Refit to unmarked subset of data
        model <- if(is.lppm(model)) update(model, UY) else update(model, as.ppp(UY))
        ## Intensity of unmarked subset
        lambdaY <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
        danger <- FALSE
      } else {
        lambdaY <- predict(model, locations=UY)
      }
    } else {
      ## ...........  Marked point process model for entire marked point pattern dataset
      if(update) {
        ## Refit to entire marked point pattern
        model <- if(is.lppm(model)) update(model, X) else update(model, as.ppp(X))
        danger <- FALSE
      }
      ## predict intensity of each type at each location
      nX <- npoints(X)
      levX <- levels(marks(X))
      nM <- length(levX)
      lamXM <- matrix(, nX, nM)
      for(m in levX) 
        lamXM[,m] <- predict(model, locations=X %mark% m)
      if(update && leaveoneout) {
        ## Use leave-one-out predictions for actual observed (location, mark)
        fitlam <- fitted(model, dataonly=TRUE, leaveoneout=TRUE)
        if(length(fitlam) != nX)
          stop(paste("Internal error: trying to overwrite",
                     nX, "predicted values with",
                     length(fitlam), "fitted values"),
               call.=FALSE)
        lamXM[cbind(1:nX, as.integer(marks(X)))] <- fitlam
      }
      ## Now marginalise to find intensity of subset
      if(!is.null(subset)) {
        if(missing(subsetlevels))
          stop(paste("Argument 'subsetlevels' is required when 'subset' is given",
                     "and 'lambda' is a multitype point process"))
        lamXM <- lamXM[, match(subsetlevels, levX), drop=FALSE]
      }
      lambdaX <- rowSums(lamXM)
      ## Restrict to subset in data, if required
      lambdaY <- if(is.null(subset)) lambdaX else lambdaX[subset]
    }
    ## Also predict everywhere if desired
    if(everywhere) everylambda <- predict(model)
  } else {
    ## lambda is some other kind of object
    lambdaY <-
      if(is.vector(lambda)) lambda  else
      if(inherits(lambda, "linfun")) lambda(Y, ...) else
      if(inherits(lambda, "linim")) lambda[Y, drop=FALSE] else
      if(is.function(lambda)) {
        coo <- coords(Y)
        do.call.matched(lambda, list(x=coo$x, y=coo$y, ...))
      } else if(is.im(lambda)) safelookup(lambda, as.ppp(Y)) else 
      stop(paste(lambdaname, "should be",
                 "a numeric vector, function, pixel image, or fitted model"))
    if(everywhere) {
      L <- domain(X)
      everylambda <-
        if(inherits(lambda, "linim")) lambda else
        if(inherits(lambda, "linfun")) as.linim(lambda, L, ...) else
        if(is.function(lambda)) as.linim(lambda, L, ...) else
        if(is.im(lambda)) as.linim(lambda, L, ...) else NULL
    }
  }
  if(!is.numeric(lambdaY))
    stop(paste("Values of", lambdaname, "are not numeric"))
  if((nv <- length(lambdaY)) != (np <- npoints(Y)))
    stop(paste("Obtained", nv, "values of", lambdaname,
	   "but point pattern contains", np, "points"))
  if(any(lambdaY < 0))
    stop(paste("Negative values of", lambdaname, "obtained"))
  if(any(lambdaY == 0))
    stop(paste("Zero values of", lambdaname, "obtained"))
  if(danger)
    attr(lambdaY, "dangerous") <- lambdaname
  if(everywhere) return(list(lambdaX=lambdaY, lambdaImage=everylambda))
  return(lambdaY)
}

