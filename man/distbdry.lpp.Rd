\name{distbdry.lpp}
\alias{distbdry.lpp}
\title{
  Distance To Boundary in a Linear Network
}
\description{
  For a point pattern on a linear network,
  compute the distance from each point
  to the nearest terminal vertex of the network.
}
\usage{
\method{distbdry}{lpp}(X)
}
\arguments{
  \item{X}{
    A point pattern on a linear network (class \code{"lpp"}).
  }
}
\details{
  This function computes the distance from each data point \code{X[i]}
  to the nearest terminal vertex of the network.
  A terminal vertex is a vertex
  which is the endpoint of exactly one segment of the network.
  Distance is measured by the length of the shortest path in the network.

  The function \code{distbdry} is generic, with methods for
  class \code{"ppp"}, \code{"pp3"} and other classes.
}
\value{
  A numeric vector with one entry for each point of \code{X}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{bdist.points}}
}
\examples{
  X <- if(interactive()) chicago else spiders
  b <- distbdry(X)
  plot(X \%mark\% b, main=expression(distbdry(X)),
       pch=16, size=1, cols=topo.colors(128))
  L <- domain(X)
  V <- vertices(L)
  d <- vertexdegree(L)
  Terminal <- V[d == 1]
  points(Terminal, pch=2, col="red")
}
\keyword{spatial}
\keyword{math}
