\name{diameter.linnet}
\alias{boundingradius.linnet}
\alias{diameter.linnet}
\title{
  Diameter and Bounding Radius of a Linear Network
}
\description{
  Compute the diameter or bounding radius of a linear network
  measured using the shortest path distance.
}
\usage{
\method{diameter}{linnet}(x)

\method{boundingradius}{linnet}(x, \dots)
}
\arguments{
  \item{x}{
    Linear network (object of class \code{"linnet"}).
  }
  \item{\dots}{Ignored.}
}
\details{
  The diameter of a linear network (in the shortest path distance)
  is the maximum value of the shortest-path distance between any
  two points \eqn{u} and \eqn{v} on the network.

  The bounding radius of a linear network (in the shortest path distance)
  is the minimum value, over all points \eqn{u} on the network,
  of the maximum shortest-path distance from \eqn{u} to another point
  \eqn{v} on the network.
  
  The functions \code{\link[spatstat.geom]{boundingradius}} and
  \code{\link[spatstat.geom]{diameter}} are generic;
  the functions \code{boundingradius.linnet} and \code{diameter.linnet}
  are the methods for objects of class \code{linnet}.
}
\value{
  A single numeric value.
}
\author{
  \adrian
}
\seealso{
 \code{\link[spatstat.geom]{boundingradius}}, \code{\link[spatstat.geom]{diameter}}, 
 \code{\link{linnet}}
}
\examples{
   diameter(simplenet)
   boundingradius(simplenet)
}
\keyword{spatial}
\keyword{math}
\concept{Linear network}
