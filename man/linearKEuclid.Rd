\name{linearKEuclid}
\alias{linearKEuclid}
\title{
  Linear K Function Using Euclidean Distance
}
\description{
  Computes an estimate of the linear \eqn{K} function
  based on Euclidean distances,
  for a point pattern on a linear network.
}
\usage{
linearKEuclid(X, r = NULL, ...)
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{r}{
    Optional. Numeric vector of values of the function argument \eqn{r}.
    There is a sensible default.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This command computes an estimate of the 
  linear \eqn{K} function based on Euclidean distances
  between the points, as described by
  Rakshit, Nair and Baddeley (2017).

  This is different from the linear \code{K} function
  based on shortest-path distances, which is computed by
  \code{\link{linearK}}.

  The linear \eqn{K} function based on Euclidean distances
  is defined in equation (20) of Rakshit, Nair and Baddeley (2017).
  The estimate is computed from the point pattern as described in equation (25).
}
\value{
  Function value table (object of class \code{"fv"}).
}
\references{
  Rakshit. S., Nair, G. and Baddeley, A. (2017)
  Second-order analysis of point patterns on a network
  using any distance metric. \emph{Spatial Statistics} \bold{22} (1) 129--154.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{linearpcfEuclid}}, \code{\link{linearKEuclidInhom}}.

  See \code{\link{linearK}} for the corresponding function
  based on shortest-path distances.
}
\examples{
  X <- rpoislpp(5, simplenet)
  K <- linearKEuclid(X)
}
\keyword{spatial}
\keyword{nonparametric}