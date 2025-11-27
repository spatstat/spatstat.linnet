\name{lineardirichlet}
\alias{lineardirichlet}
\title{
  Dirichlet Tessellation on a Linear Network
}
\description{
  Given a point pattern on a linear network, compute the Dirichlet
  (or Voronoi or Thiessen) tessellation induced by the points.
}
\usage{
lineardirichlet(X, metric=c("shortestpath", "Euclidean"))
}
\arguments{
  \item{X}{
    Point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{metric}{
    Character string (partially matched) specifying the distance
    metric used to define the Dirichlet tessellation.
  }
}
\details{
  The Dirichlet tessellation induced by a point pattern \code{X}
  on a linear network \code{L}
  is a partition of \code{L} into subsets. The subset \code{L[i]}
  associated with the data point \code{X[i]} is the part of \code{L}
  lying closer to \code{X[i]} than to any other data point \code{X[j]}.

  If \code{metric="shortestpath"} (the default),
  distance between points on the network
  is measured by the shortest path in the network.
  If \code{metric="Euclidean"}, distance is measured by the
  Euclidean distance in two dimensions.
}
\section{Missing tiles}{
  If the linear network is not connected, and if one of the connected
  components contains no data points, then the Dirichlet tessellation
  using \code{metric="shortestpath"}
  is mathematically undefined inside this component.
  The resulting tessellation object includes
  a tile with label \code{NA}, which contains this component of the network.
  A plot of the tessellation will not show this tile. 
}
\value{
  A tessellation on a linear network
  (object of class \code{"lintess"}).
}
\author{
  \adrian.
}
\seealso{
  \code{\link{lintess}}.

  For the Dirichlet tessellation in two-dimensional space,
  see \code{\link[spatstat.geom]{dirichlet}}.
}
\examples{
  X <- runiflpp(5, simplenet)
  plot(lineardirichlet(X), lwd=3)
  points(X)
  plot(lineardirichlet(X, metric="E"), lwd=3)
  points(X)
  plot(dirichlet(as.ppp(X)), add=TRUE, lty=2)
}
\keyword{spatial}
\keyword{manip}
\concept{Dirichlet tessellation}
\concept{Linear network}
