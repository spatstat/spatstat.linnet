\name{linearpcfEuclid}
\alias{linearpcfEuclid}
\title{
  Linear Pair Correlation Function Using Euclidean Distance
}
\description{
  Computes an estimate of the pair correlation function
  based on Euclidean distances,
  for a point pattern on a linear network.
}
\usage{
linearpcfEuclid(X, r = NULL, ...)
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
  This command computes an estimate of the pair correlation function
  based on Euclidean distances between the points, as described by
  Rakshit, Nair and Baddeley (2017).

  This is different from the linear pair correlation function
  based on shortest-path distances, which is computed by \code{\link{linearpcf}}.

  The linear pair correlation function based on Euclidean distances
  is defined in equation (15) of Rakshit, Nair and Baddeley (2017).
  The estimate is computed from the point pattern as described in equation (31).
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
  \code{\link{linearKEuclid}}, \code{\link{linearpcfEuclidInhom}}.

  See \code{\link{linearpcf}} for the corresponding function
  based on shortest-path distances.
}
\examples{
  X <- rpoislpp(5, simplenet)
  g <- linearpcfEuclid(X)
}
\keyword{spatial}
\keyword{nonparametric}