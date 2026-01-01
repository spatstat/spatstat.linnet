\name{roadways}
\alias{roadways}
\title{
  Identify Uninterrupted Roadways in a Linear Network
}
\description{
  Given a linear network, find the roadways.
  A roadway is a path in the network which does not pass through
  any forks in the network.
}
\usage{
roadways(X, what = c("labels", "segments", "function"))
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
  A \sQuote{roadway} in a network is a path which does not
  pass through any forks in the network, that is, a path which only
  passes through vertices of degree 2.

  Given a linear network \code{X},
  this function identifies the roadways in \code{X} and assigns a label to each
  edge (segment) of the network indicating which roadway it belongs to.
  
  Formally a roadway is a path joining successive
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

  Every edge (segment) of the network belomgs to a unique roadway.
  This algorithm assigns a label to each edge indicating which roadway
  it belongs to. If \code{what="labels"} (the default), the result is a
  factor of length equal to \code{nsegments(X)} giving the
  classification of segments. If \code{what="segments"}, the result is
  a line segment pattern (object of class \code{"psp"}) with marks
  indicating the classification. If \code{what="function"}, the result
  is a function of class \code{"linfun"} which can be applied to any
  location on the network to determine which roadway it belongs to.
}
\value{
  A factor, or an object of class \code{"psp"}, or a function
  of class \code{"linfun"}.
}
\author{
  \adrian
}
\seealso{
  \code{\link{connected.linnet}}
}
\examples{
  A <- edges(letterR)
  B <- shift(rotate(A, 0.15), c(-0.13, 0))
  L <- as.linnet(superimpose(A,B))
  plot(roadways(L, what="s"),
       col=rainbow, scramble.cols=TRUE,
       lwd=2, show.window=FALSE)
}
\keyword{spatial}
\keyword{manip}
\concept{Linear network}
