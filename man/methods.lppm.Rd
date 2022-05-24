\name{methods.lppm}
\alias{methods.lppm} %DoNotExport
\alias{coef.lppm}
\alias{emend.lppm}
\alias{extractAIC.lppm}
\alias{formula.lppm}
\alias{logLik.lppm}
\alias{deviance.lppm}
\alias{nobs.lppm}
\alias{print.lppm}
\alias{summary.lppm}
\alias{terms.lppm}
\alias{update.lppm}
\alias{valid.lppm}
\alias{vcov.lppm}
\alias{as.linnet.lppm}
\alias{response.lppm}
\title{
  Methods for Fitted Point Process Models on a Linear Network
}
\description{
  These are methods for the class \code{"lppm"} of fitted point process
  models on a linear network.
}
\usage{
  \method{coef}{lppm}(object, ...)

  \method{emend}{lppm}(object, \dots)

  \method{extractAIC}{lppm}(fit, ...)

  \method{formula}{lppm}(x, ...)

  \method{logLik}{lppm}(object, ...)

  \method{deviance}{lppm}(object, ...)

  \method{nobs}{lppm}(object, ...)

  \method{print}{lppm}(x, ...)

  \method{summary}{lppm}(object, ...)

  \method{terms}{lppm}(x, ...)

  \method{update}{lppm}(object, ...)

  \method{valid}{lppm}(object, ...)

  \method{vcov}{lppm}(object, ...)

  \method{as.linnet}{lppm}(X, ...)

  \method{response}{lppm}(object)
}
\arguments{
  \item{object,fit,x,X}{
    An object of class \code{"lppm"} representing a fitted point process
    model on a linear network.
  }
  \item{\dots}{
    Arguments passed to other methods, usually the
    method for the class \code{"ppm"}.
  }
}
\details{
  These are methods for the \R generic commands
  \code{\link[stats]{coef}},
  \code{\link[stats]{extractAIC}},
  \code{\link[stats]{formula}},
  \code{\link[stats]{logLik}},
  \code{\link[stats]{deviance}},
  \code{\link[stats]{nobs}},
  \code{\link[base]{print}}, 
  \code{\link[base]{summary}}, 
  \code{\link[stats]{terms}},
  \code{\link[stats]{update}} and
  \code{\link[stats]{vcov}},
  and the \pkg{spatstat} generic commands
  \code{\link{as.linnet}},
  \code{\link[spatstat.model]{emend}},
  \code{\link[spatstat.model]{response}} and 
  \code{\link[spatstat.model]{valid}},
  for the class \code{"lppm"}. 
}
\value{
  For \code{as.linnet.lppm} a linear network (object of class \code{"linnet"}).
  For \code{emend.lppm} another fitted model of the same class \code{"lppm"}.
  For \code{response.lppm} a spatial point pattern on a linear network
  (object of class \code{"lpp"}).
  For \code{valid.lppm} a logical value.

  For the other methods, see the help for the default methods.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{lppm}},
  \code{\link{plot.lppm}}.
}
\examples{
  X <- runiflpp(15, simplenet)
  fit <- lppm(X ~ x)
  print(fit)
  coef(fit)
  formula(fit)
  terms(fit)
  logLik(fit)
  deviance(fit)
  nobs(fit)
  extractAIC(fit)
  update(fit, ~1)
  valid(fit)
  vcov(fit)
  response(fit)
}
\keyword{spatial}
\keyword{models}
