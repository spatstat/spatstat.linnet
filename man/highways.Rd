\name{highways}
\alias{highways}
\title{
  Identify Uninterrupted Highways in a Linear Network
}
\description{
  Given a linear network, find the highways.
  A highway is a set of connected edges which do not pass through
  any forks in the network.
}
\usage{
highways(X, what = c("tessellation", "labels"))
}
\arguments{
  \item{X}{
    A linear network (object of class \code{"linnet"}).
  }
  \item{what}{
    String (partially matched) specifying the kind of result.
  }
}
\details{
  A \sQuote{highway} in a network is a subset of the edges of the network
  which does not pass through any forks in the network.

  Given a linear network \code{X},
  this function identifies the highways in \code{X} and assigns a label to each
  edge (segment) of the network indicating which highway it belongs to.
  
  Formally a highway is a set of edges joining successive
  vertices \eqn{v_1, \ldots, v_n}{v[1], ..., v[n]} of the network
  such that
  \itemize{
    \item there is an edge of the network joining
    \eqn{v_i}{v[i]} and \eqn{v_{i+1}}{v[i+1]} for each \eqn{i}
    \item each of the intermediate vertices
    vertices \eqn{v_2, \ldots, v_{n-1}}{v[2], ..., v[n-1]}
    has degree 2 (i.e. there are exactly 2 edges of the network which
    end at the vertex)
    \item the terminal vertices \eqn{v_1}{v[1]} and \eqn{v_n}{v[n]}
    do not have degree 2, unless they are the same vertex.
  }

  Every edge (segment) of the network belomgs to a unique highway.
  This algorithm assigns a label to each edge indicating which highway
  it belongs to. If \code{what="labels"} (the default), the result is a
  factor of length equal to \code{nsegments(X)} giving the
  classification of segments.
  If \code{what="tessellation"}, the result is a tessellation of
  \code{X} (object of class \code{"lintess"}) in which the tiles of the
  tessellation are the highways.
}
\value{
  A factor, or a tessellation on the network (object of class \code{"lintess"}).
}
\author{
  \adrian
}
\seealso{
  \code{\link{connected.linnet}}
}
\examples{
  A <- edges(letterR)
  B <- edges(disc(centre=centroid.owin(letterR)))
  L <- as.linnet(superimpose(A,B))
  plot(highways(L, what="tess"),
       col=rainbow, scramble.cols=TRUE,
       lwd=2, show.window=FALSE)
}
\keyword{spatial}
\keyword{manip}
\concept{Linear network}
