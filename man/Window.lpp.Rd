\name{Window.lpp}
\alias{Window.lpp}
\alias{Window.lppm}
\title{Extract Window of Spatial Object on a Network}
\description{
  Given a spatial object on a network,
  these functions extract the window
  in which the network is defined.
}
\usage{
 \method{Window}{lpp}(X, \dots)

 \method{Window}{lppm}(X, \dots)
}
\arguments{
  \item{X}{A spatial object.}
  \item{\dots}{Ignored.}
}
\value{
  An object of class \code{"owin"}
  (see \code{\link[spatstat.geom]{owin.object}})
  specifying an observation window.
}
\details{
  These are methods for the generic function \code{\link[spatstat.geom]{Window}}
  which extract the spatial window in which the object \code{X}
  is defined. 

  For the methods defined here, \code{X} should be a spatial object
  on a linear network (object of class \code{"lpp"} or \code{"lppm"}).
}
\seealso{
  \code{\link[spatstat.geom]{Window}}.
}
\examples{
   Window(spiders)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
 
 
\concept{Linear network}
