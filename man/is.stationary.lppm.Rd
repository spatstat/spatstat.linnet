\name{is.stationary.lppm}
\alias{is.stationary.lppm}
\alias{is.poisson.lppm}
\title{
  Recognise Stationary and Poisson Point Process Models on a Network
}
\description{
  Given a point process model that has been fitted to data on a network,
  determine whether the model is a stationary point process,
  and whether it is a Poisson point process.
}
\usage{
\method{is.stationary}{lppm}(x)

\method{is.poisson}{lppm}(x)
}
\arguments{
  \item{x}{
    A fitted spatial point process model on a linear network
    (object of class \code{"lppm"}).
  }
}
\details{
  The argument \code{x} represents a fitted spatial point process model
  on a linear network.

  \code{is.stationary(x)} returns \code{TRUE} if \code{x} represents
  a stationary point process, and \code{FALSE} if not.

  \code{is.poisson(x)} returns \code{TRUE} if \code{x} represents
  a Poisson point process, and \code{FALSE} if not.

  The functions \code{\link[spatstat.core]{is.stationary}}
  and \code{\link[spatstat.core]{is.poisson}} are generic,
  with methods for many classes of models.
}
\value{
  A logical value.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[spatstat.core]{is.marked}}
  to determine whether a model is a marked
  point process. 
  
  \code{\link{summary.lppm}} for detailed information.

  Model-fitting function
  \code{\link{lppm}}.
}
\examples{
  fit <- lppm(spiders ~ x)
  is.stationary(fit)
  is.poisson(fit)
}
\keyword{spatial}
\keyword{models}
