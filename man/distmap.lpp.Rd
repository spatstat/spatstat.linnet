\name{distmap.lpp}  
\alias{distmap.lpp}
\title{
  Distance Map of Point Pattern on Linear Network
}
\description{
  Computes the distance from each pixel to the nearest
  point in the given point pattern on a linear network.
}
\usage{
  \method{distmap}{lpp}(X, \dots, k=1)
}
\arguments{
  \item{X}{A point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{k}{
    Integer. The distance to the \code{k}-th nearest data point
    will be computed.
  }
  \item{\dots}{Arguments passed to \code{\link{as.linim.linfun}}
    to control pixel resolution.
  }
}
\value{
  A pixel image on the network (object of class \code{"linim"})
  whose greyscale values are the values of the distance map.
}
\details{
  This is a method for the generic function \code{\link{distmap}}.
  It computes the distance map of the point pattern \code{X}
  as a pixel image on the network.

  At a pixel \eqn{u}, the greyscale value
  equals the distance from \eqn{u} to the
  nearest point of the pattern \code{X} (or the \eqn{k}-th nearest point
  of \code{X}).
}
\seealso{
  Generic function \code{\link{distmap}} and methods.
}
\examples{
  plot(distmap(spiders))
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{math}
