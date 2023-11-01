\name{intersect.lintess}
\alias{intersect.lintess}
\title{
  Intersection of Tessellations on a Linear Network
}
\description{
  Yields the intersection (common refinement) of two tessellations
  on a linear network.
}
\usage{
intersect.lintess(X, Y)
}
\arguments{
  \item{X,Y}{
    Tessellations (objects of class \code{"lintess"})
    on the same linear network, or data that define such tessellations.
    See Details.
  }
}
\value{
  Another tessellation (object of class \code{"lintess"})
  on the same linear network as \code{X} and \code{Y}.
}
\details{
  \code{X} and \code{Y} should be tessellations on a linear network
  (objects of class \code{"lintess"}) and should be defined on the
  same network. The algorithm finds the common refinement of the
  two tessellations. Each tile in the resulting tessellation
  is the intersection of a tile of \code{X} with a tile of \code{Y}.

  Alternatively, one of the arguments \code{X} or \code{Y} 
  can be a two-dimensional tessellation (object of class
  \code{"tess"}) while the other argument is a network or a tessellation
  on a network. The two-dimensional tessellation will be intersected
  with the network to produce a tessellation on the network, then
  intersected with the other tessellation on the network.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{lintess}},
  \code{\link{divide.linnet}},
  \code{\link{chop.linnet}}
}
\examples{
  X <- divide.linnet(runiflpp(4, simplenet))
  Y <- divide.linnet(runiflpp(3, simplenet))
  opa <- par(mfrow=c(1,3))
  plot(X)
  plot(Y)
  plot(intersect.lintess(X,Y))
  par(opa)
}
\keyword{spatial}
\keyword{manip}
\concept{Tessellation}
\concept{Linear network}
