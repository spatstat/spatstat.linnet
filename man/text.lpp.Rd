\name{text.lpp}
\alias{text.lpp}
\title{
  Add Text Labels to Point Pattern on a Network
}
\description{
  Plots a text label at the location of each point,
  for a point pattern on a linear network.
}
\usage{
\method{text}{lpp}(x, \dots)
}
\arguments{
  \item{x}{
    A point pattern on a linear network (class \code{"lpp"}).
  }
  \item{\dots}{
   Additional arguments passed to \code{\link[graphics]{text.default}}.
  }
}
\details{
  This function is a method for the generic \code{\link[graphics]{text}}.
  A text label is added to the existing plot,
  at the location of each point in the point pattern \code{x},
  or near the location of the midpoint of each segment
  in the segment pattern \code{x}.

  Additional arguments \code{\dots} are passed to
  \code{\link[graphics]{text.default}} and may be used to control the
  placement of the labels relative to the point locations, and
  the size and colour of the labels.

  By default, the labels are the serial numbers 1 to \code{n}, where
  \code{n} is the number of points or segments in \code{x}. This can
  be changed by specifying the argument \code{labels}, which should be
  a vector of length \code{n}.
}
\value{
  Null.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[graphics]{text.default}},
  \code{\link[spatstat.geom]{text.ppp}}
}
\examples{
  X <- runiflpp(5, simplenet)
  plot(X)
  text(X, pos=2, col="blue")
}
\keyword{spatial}
\keyword{hplot}
\concept{Linear network}
