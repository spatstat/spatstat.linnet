\name{shortestpath}
\alias{shortestpath}
\title{
Shortest Path Between Two Points on a Linear Network
}
\description{
  Find the shortest path between two given points on a linear network.
}
\usage{
shortestpath(X, i=1, j=2)
}
\arguments{
  \item{X}{
    Point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{i}{Integer index of the start point}
  \item{j}{Integer index of the end point}
}
\details{
  The shortest path in the network between
  the two specified points \code{X[i]} and \code{X[j]}
  is determined.

  The result is a line segment pattern
  (object of class \code{"psp"}) consisting of (in order) a line segment
  joining \code{X[i]} to a vertex of the network, then a series of
  segments joining adjacent vertices
  of the network, then a line segment joining a vertex to \code{X[j]}.

  The result has an attribute \code{"steps"} which is an integer vector
  giving the sequence of vertices of the network
  through which the shortest path passes.
  This vector will have length zero if \code{X[i]} and
  \code{X[j]} lie on the same segment of the network.
  Otherwise it will contain a sequence of integers
  which index the vertices as given in \code{vertices(domain(X))}.
}
\value{
  Line segment pattern (object of class \code{"psp"}) with an 
  attribute \code{"steps"} which is an integer vector.
}
\author{
  \adrian.
}
\seealso{
  The \emph{length} of the shortest path is computed by
  \code{\link{pairdist.lpp}}.
}
\examples{
  X <- if(interactive()) chicago[c(20, 39)] else runiflpp(2, simplenet)
  P <- shortestpath(X)
  if(interactive()) {
    plot(X, pch=16, main="shortest path")
    ## draw the path
    plot(P, add=TRUE, col=2, lwd=2)
    ## draw the vertices on the path
    V <- vertices(L)
    steps <- attr(P, "steps")
    plot(V[steps], add=TRUE, col=3)
  }
}
\keyword{spatial}
\keyword{math}
\concept{Linear network}
