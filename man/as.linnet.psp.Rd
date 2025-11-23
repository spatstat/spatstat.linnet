\name{as.linnet.psp}
\alias{as.linnet.psp}
\title{
  Convert Line Segment Pattern to Linear Network
}
\description{
  Converts a line segment pattern to a linear network.
}
\usage{
 \method{as.linnet}{psp}(X, \dots, eps, sparse=FALSE, chop=TRUE, fuse=TRUE)
}
\arguments{
  \item{X}{
    Line segment pattern (object of class \code{"psp"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{eps}{
    Optional. Distance threshold. If two segment endpoints
    are closer than \code{eps} units apart, they will be treated
    as the same point, and will become a single vertex in the
    linear network. Ignored if \code{fuse=FALSE}.
  }
  \item{sparse}{
    Logical value indicating whether to use a sparse matrix
    representation, as explained in \code{\link{linnet}}.
  }
  \item{chop}{
    Logical value specifying whether segments which cross each other
    should be subdivided into separate segments of the network which
    meet at the intersection point.
  }
  \item{fuse}{
    Logical value specifying whether two vertices lyin closer than \code{eps}
    should be treated as a single vertex.
  }
}
\details{
  This command converts any collection of line segments into a linear
  network by guessing the connectivity of the network.

  If \code{chop=TRUE} (the default), then
  if any segments in \code{X} cross over each other, they are first
  cut into pieces using \code{\link[spatstat.geom]{selfcut.psp}}.

  If \code{fuse=TRUE} (the default), 
  then any pair of segment endpoints lying closer than \code{eps}
  units apart, is treated as a single vertex.

  After these modifications, the linear network
  is constructed using \code{\link{linnet}}.

  If \code{chop=FALSE} and \code{fuse=FALSE}, each segment in the
  segment pattern \code{X} becomes an edge in the resulting network,
  and no other edges or vertices are created.

  It would be wise to check the result by plotting the degree
  of each vertex, as shown in the Examples.

  If \code{X} has marks, then these are stored in
  the resulting linear network \code{Y <- as.linnet(X)},
  and can be extracted as \code{marks(as.psp(Y))} or \code{marks(Y$lines)}.
}
\value{
A linear network (object of class \code{"linnet"}).

The result also has an attribute \code{"camefrom"} indicating
the provenance of each line in the resulting network.
For example \code{camefrom[3]=2} means that the third line segment
in the result is a piece of the second segment of \code{X}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{linnet}}, \code{\link[spatstat.geom]{selfcut.psp}},
  \code{\link{methods.linnet}}.
}
\examples{
  # make some data
  A <- psp(0.09, 0.55, 0.79, 0.80, window=owin())
  B <- superimpose(A, as.psp(simplenet))

  # convert to a linear network
  L <- as.linnet(B)

  # check validity
  L
  plot(L)
  text(vertices(L), labels=vertexdegree(L))

  # show the pieces that came from original segment number 1
  S <- as.psp(L)
  (camefrom <- attr(L, "camefrom"))
  parts <- which(camefrom == 1)
  plot(S[parts], add=TRUE, col="green", lwd=2)

  # convert to a network without changing the geometry
  H <- as.linnet(B, chop=FALSE, fuse=FALSE)
}
\keyword{spatial}
\keyword{manip}
\concept{Linear network}
