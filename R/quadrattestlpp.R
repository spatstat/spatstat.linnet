#'
#'  quadrattestlpp.R
#'
#'  Quadrat counting tests for lpp and lppm
#'
#'  $Revision: 1.1 $ $Date: 2025/11/17 09:34:34 $

quadrat.test.lpp <-
  function(X, ...,
           tess=NULL,
           nx=5, ny=nx,
           xbreaks=NULL, ybreaks=NULL,
           alternative = c("two.sided", "regular", "clustered"),
           method = c("Chisq", "MonteCarlo"),
           conditional=TRUE, CR=1,
           lambda=NULL, df.est=NULL,
           nsim=1999)
{
   Xname <- short.deparse(substitute(X))
   method <- match.arg(method)
   alternative <- match.arg(alternative)
   do.call(quadrat.testEngine,
           resolve.defaults(list(quote(X),
                                tess=tess,
                                nx=nx, ny=ny,
                                xbreaks=xbreaks, ybreaks=ybreaks,
                                alternative=alternative,
                                method=method,
                                conditional=conditional,
                                CR=CR,
                                fit=lambda,
                                df.est=df.est,
                                nsim=nsim),
                           list(...), 
                           list(Xname=Xname, fitname="CSR on a network")))
}

quadrat.test.lppm <- 
  function(X, ...,
           tess=NULL,
           nx=5, ny=nx,
           xbreaks=NULL, ybreaks=NULL,
           alternative = c("two.sided", "regular", "clustered"),      
           method=c("Chisq", "MonteCarlo"),
           conditional=TRUE, CR=1, df.est=NULL, 
           nsim=1999)
{
   fitname <- short.deparse(substitute(X))
   dataname <- paste("data from", fitname)
   method <- match.arg(method)
   alternative <- match.arg(alternative)
   if(!is.poisson(X))
    stop("Test is only defined for Poisson point process models")
   if(is.marked(X))
    stop("Sorry, not yet implemented for marked point process models")
   Xdata <- response(X)
   dont.complain.about(Xdata)
   do.call(quadrat.testEngine,
          resolve.defaults(list(quote(Xdata), 
                                tess=tess,
                                nx=nx, ny=ny,
                                xbreaks=xbreaks, ybreaks=ybreaks,
                                alternative=alternative,
                                method=method,
                                conditional=conditional, CR=CR,
                                nsim=nsim, 
                                fit=X,
                                df.est=df.est),
                           list(...),
                           list(Xname=dataname, fitname=fitname)))
}

quadrat.test.linearquadratcount <-
  function(X, ...,
           alternative = c("two.sided", "regular", "clustered"),
           method=c("Chisq", "MonteCarlo"),
           conditional=TRUE, CR=1,
           lambda=NULL, df.est=NULL,
           nsim=1999) {
   trap.extra.arguments(...)
   method <- match.arg(method)
   alternative <- match.arg(alternative)
   quadrat.testEngine(Xcount=X,
                      alternative=alternative,
                      fit=lambda, df.est=df.est,
                      method=method, conditional=conditional, CR=CR, nsim=nsim)
}
