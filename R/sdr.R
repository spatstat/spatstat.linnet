#'  Sufficient Dimension Reduction 
#'    including linear network case

sdr.lpp <- local({

  sdr.lpp <- function(X, covariates,
                      method=c("DR", "NNIR", "SAVE", "SIR", "TSE"),
                      Dim1=1, Dim2=1, predict=FALSE, ...) {
    stopifnot(is.ppp(X) || is.lpp(X))
    method <- match.arg(method)
    trap.extra.arguments(...)
    #' ensure 'covariates' is a list of compatible images
    if(!inherits(covariates, "imlist") && !all(sapply(covariates, is.im)))
      stop("Argument 'covariates' must be a list of images")
    nc <- length(covariates)
    if(nc == 0)
      stop("Need at least one covariate!")
    if(nc < Dim1 + (method == "TSE") * Dim2)
      stop(paste(if(method == "TSE") "Dim1 + Dim2" else "Dim1",
                 "must not exceed the number of covariates"),
           call.=FALSE)
    #' extract "all" values
    covariates <- as.solist(covariates)
    Ypixval <- pairs(covariates, plot=FALSE)
    Ypixval <- Ypixval[complete.cases(Ypixval), , drop=FALSE]
    #' compute sample mean and covariance matrix
    m <- colMeans(Ypixval)
    V <- cov(Ypixval)
    #' evaluate each image at point data locations
    YX <- sapply(covariates, safelook, Y=X)
    #' apply precomputed standardisation
    Zx <- t(t(YX) - m) %*% matrixinvsqrt(V)
    #' ready
    coordsX <- coords(X)
    result <-
      switch(method,
             DR   =   calc.DR(COV=V, z=Zx,              Dim=Dim1),
             NNIR = calc.NNIR(COV=V, z=Zx, pos=coordsX, Dim=Dim1),
             SAVE = calc.SAVE(COV=V, z=Zx,              Dim=Dim1),
             SIR  =  calc.SIR(COV=V, z=Zx                      ),
             TSE  =  calc.TSE(COV=V, z=Zx, pos=coordsX, Dim1=Dim1, Dim2=Dim2)
             )
    #'
    covnames <- names(covariates) %orifnull% paste0("Y", 1:nc)
    dimnames(result$B) <- list(covnames, paste0("B", 1:ncol(result$B)))
    if(method == "TSE") {
      result$M1 <- namez(result$M1)
      result$M2 <- namez(result$M2)
    } else {
      result$M <- namez(result$M)
    }
    if(predict) result$Y <- sdrPredict(covariates, result$B)
    return(result)
  }

  safelook <- function(Z, Y, ...) {
    if(is.im(Z)) safelookup(Z, Y, ...) else Z[Y]
  }

  namez <- function(M, prefix="Z") {
    dimnames(M) <- list(paste0(prefix, 1:nrow(M)),
                     paste0(prefix, 1:ncol(M)))
    return(M)
  }

  sdr.lpp
})
