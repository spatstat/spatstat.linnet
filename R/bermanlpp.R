#
# bermanlpp.R
#
# Tests from Berman (1986) adapted to linear networks
#
#  $Revision: 1.3 $  $Date: 2022/01/04 05:30:06 $
#
#

berman.test.lpp <-
  function(X, covariate,
           which=c("Z1", "Z2"),
           alternative=c("two.sided", "less", "greater"),
           ...) {
    Xname <- short.deparse(substitute(X))
    covname <- short.deparse(substitute(covariate))
    force(covariate)
    if(is.character(covariate)) covname <- covariate
    which <- match.arg(which)
    alternative <- match.arg(alternative)

    model <- lppm(X)
    dont.complain.about(model)
    
    do.call(bermantestEngine,
            resolve.defaults(list(quote(model),
                                  quote(covariate),
                                  which, alternative),
                             list(...),
                             list(modelname="CSR",
                                  covname=covname, dataname=Xname)))
}

berman.test.lppm <- function(model, covariate,
                           which=c("Z1", "Z2"),
                           alternative=c("two.sided", "less", "greater"),
                           ...) {
  modelname <- short.deparse(substitute(model))
  covname <- short.deparse(substitute(covariate))
  force(model)
  force(covariate)
  if(is.character(covariate)) covname <- covariate
  verifyclass(model, "lppm")
  which <- match.arg(which)
  alternative <- match.arg(alternative)
  if(is.poisson(model) && is.stationary(model))
    modelname <- "CSR"
  do.call(bermantestEngine,
          resolve.defaults(list(quote(model),
                                quote(covariate),
                                which, alternative),
                           list(...),
                           list(modelname=modelname,
                                covname=covname,
                                dataname=model$Xname)))
}
