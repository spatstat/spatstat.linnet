\name{rjitter.lpp}
\alias{rjitter.lpp}
\title{Random Perturbation of a Point Pattern on a Network}
\description{
  Applies independent random displacements to each point
  in a point pattern on a network.
}
\usage{
\method{rjitter}{lpp}(X, radius, \dots, nsim = 1, drop = TRUE)
}
\arguments{
  \item{X}{
    A point pattern on a linear network
    (object of class \code{"lpp"}).
  }
  \item{radius}{
    Scale of perturbations. A positive numerical value.
    Each point will be displaced by a random distance,
    with maximum displacement equal to this value.
  }
  \item{\dots}{
    Ignored.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\details{
  The function \code{\link[spatstat.geom]{rjitter}} is generic. This function is
  the method for the class \code{"lpp"} of point patterns on a linear network.

  Each of the points in \code{X} will be displaced along the network
  by a random amount, independently of other points.
  The maximum displacement distance is specified by \code{radius}.
  Each point remains on the same line segment of the network
  as it originally was.
}
\value{
  A point pattern on a linear network
  (object of class \code{"lpp"})
  or a list of such point patterns.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[spatstat.geom]{rjitter}} for point patterns in two dimensions.
}
\examples{
   X <- runiflpp(3, simplenet)
   plot(X, pch=16)
   Y <- rjitter(X, 0.1)
   plot(Y, add=TRUE, cols=3)
}
\keyword{spatial}
\keyword{datagen}
\concept{Linear network}
