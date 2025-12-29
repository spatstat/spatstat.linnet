#'
#'   Header for all (concatenated) test files
#'
#'   Require spatstat.linnet
#'   Obtain environment variable controlling tests.
#'
#'   $Revision: 1.5 $ $Date: 2020/04/30 05:31:37 $

require(spatstat.linnet)
FULLTEST <- (nchar(Sys.getenv("SPATSTAT_TEST", unset="")) > 0)
ALWAYS   <- TRUE
cat(paste("--------- Executing",
          if(FULLTEST) "** ALL **" else "**RESTRICTED** subset of",
          "test code -----------\n"))
##
##  tests/segments.R
##   Tests of psp class and related code
##                      [SEE ALSO: tests/xysegment.R]
##
##  $Revision: 1.33 $  $Date: 2022/05/22 08:39:47 $




if(FULLTEST) {
  local({
    ## more tests of lppm code
    fit <- lppm(unmark(chicago) ~ polynom(x,y,2))
    Z <- predict(fit)
  })
}

#'
#'   spatstat.linnet/tests/sidefx.R
#'
#' Test whether plot(do.plot=FALSE) has no side effects on graphics system
#'
#'  $Revision: 1.6 $  $Date: 2025/12/23 00:24:24 $

local({
  if(FULLTEST) {
    ## test whether a graphics device has been started
    deviceExists <- function() { length(dev.list()) != 0 }
    ## check whether executing 'expr' causes creation of a graphics device
    chk <- function(expr) {
      ename <- sQuote(deparse(substitute(expr)))
      if(deviceExists()) {
        ## try switching off the graphics
        graphics.off()
        if(deviceExists()) {
          ## Dang
          warning(paste("Cannot check", ename, 
                        "as a graphics device already exists"),
                  call.=FALSE)
          return(FALSE)
        }
      }
      eval(expr)
      if(deviceExists()) {
        stop(paste("Evaluating", ename, 
                   "caused a graphics device to be started"),
             call.=FALSE)
      }
      return(TRUE)
    }
    




    ## linnet
    chk(plot(simplenet, do.plot=FALSE))
    ## lpp
    chk(plot(spiders, do.plot=FALSE))
    chk(plot(chicago, do.plot=FALSE))
    ## lintess
    chk(plot(as.lintess(simplenet), do.plot=FALSE))
    ## linfun, linim
    f <- function(x,y,seg,tp) { x - y }
    aLinfun <- Y <- linfun(f, simplenet)
    aLinim  <- Z <- as.linim(Y)
    chk(plot(aLinfun, do.plot=FALSE))    
    chk(plot(aLinim, do.plot=FALSE))
    ## linearquadratcount
    chk(plot(quadratcount(spiders, nx=3), do.plot=FALSE))

  }
})
#
#  tests/undoc.R
#
#   $Revision: 1.17 $   $Date: 2025/07/21 07:35:03 $
#
#  Test undocumented hacks, experimental code, etc


if(FULLTEST) {
  local({
    ## linim helper functions
    df <- pointsAlongNetwork(simplenet, 0.2)
  })
}
  
##
##  tests/updateppm.R
##
##  Check validity of update.ppm
##
##  $Revision: 1.9 $ $Date: 2025/12/19 00:37:21 $

local({
  
  if(ALWAYS) {
    ## test update.lppm
    X <- runiflpp(20, simplenet)
    fit0 <- lppm(X ~ 1)
    fit1 <- update(fit0, ~ x)
    anova(fit0, fit1, test="LR")
    cat("update.lppm(fit, ~trend) is OK\n")
    fit2 <- update(fit0, . ~ x)
    anova(fit0, fit2, test="LR")
    cat("update.lppm(fit, . ~ trend) is OK\n")
  }
})
##
## tests/xysegment.R
##                      [SEE ALSO tests/segments.R]
##
##    Test weird problems and boundary cases for line segment code
##
##    $Version$ $Date: 2022/10/23 01:21:09 $ 
##


if(FULLTEST) local({ S <- as.psp(simplenet) })

