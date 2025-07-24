\name{identify.linnet}
\alias{identify.linnet}
\title{Interactively Identify Segments of a Linear Network}
\description{
  If a linear network is plotted in the graphics window,
  then each time the left mouse button is pressed,
  this function will find the network segment which is closest
  the mouse position, and print its serial number.
}
\usage{
  \method{identify}{linnet}(x, \dots)
}
\arguments{
  \item{x}{
    A linear network
    (object of class \code{"linnet"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{identify.psp}} and ultimately to
    \code{\link[graphics]{identify.default}}.
  }
}
\value{
  A vector containing the serial numbers of the network segments of
  \code{x} that were identified.
}
\details{
  This is a method for the generic function \code{\link[graphics]{identify}}
  for linear networks. 

  The network \code{x} should first be plotted
  using \code{\link{plot.linnet}}, \code{\link{plot.lpp}}
  or \code{\link{plot.linim}}.
  Then \code{identify(x)}
  reads the position of the graphics pointer each time the
  left mouse button is pressed.  It then determines which network segment
  lies closest to the mouse position.
  The index of this segment (and its mark if any) 
  will be returned as part of the value of the call.

  Each time a segment is identified,
  text will be displayed at the midpoint of the segment,
  showing its serial number.
}
\seealso{
  \code{\link[graphics]{identify}},
  \code{\link[spatstat.geom]{identify.psp}},
  \code{\link{plot.linnet}}.
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{iplot}
