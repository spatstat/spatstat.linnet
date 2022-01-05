\name{domain.lpp}
\alias{domain.lpp}
\alias{domain.lppm}
\alias{domain.linfun}
\alias{domain.lintess}
\title{
  Extract the Linear Network on which Spatial Data are Defined
}
\description{
  Given a spatial object representing data on a linear network,
  extract the network.
}
\usage{
 \method{domain}{lpp}(X, \dots)

 \method{domain}{lppm}(X, \dots)

 \method{domain}{linfun}(X, \dots)

 \method{domain}{lintess}(X, \dots)
}
\arguments{
  \item{X}{
    A spatial object representing data on a linear network.
    An object of class \code{"lpp"}, \code{"lppm"}, \code{"linfun"}
    or \code{"lintess"}.
  }
  \item{\dots}{
    Extra arguments. They are ignored by all the methods listed here.
  }
}
\details{
  The function \code{\link[spatstat.geom]{domain}} is generic,
  with methods for many classes.

  For a spatial object \code{X} 
  \code{domain(X)} extracts the spatial domain in which \code{X} is
  defined.

  For a two-dimensional object \code{X}, typically \code{domain(X)}
  is the same as \code{Window(X)}.

  The exception is that, if \code{X} is a point pattern on a linear network
  (class \code{"lpp"}) or a point process model on a linear network
  (class \code{"lppm"}), then \code{domain(X)} is the linear network
  on which the points lie, while \code{Window(X)} is the two-dimensional
  window containing the linear network.
}
\value{
  A linear network (object of class \code{"linnet"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[spatstat.geom]{domain}},
  \code{\link[spatstat.random]{domain.rmhmodel}},
  \code{\link[spatstat.core]{domain.ppm}}.
  
  \code{\link[spatstat.geom]{Window}},
  \code{\link[spatstat.geom]{Frame}}
}
\examples{
  domain(chicago)
}
\keyword{spatial}
\keyword{manip}
