\name{intensity.lpp}  
\alias{intensity.lpp}
\title{
  Empirical Intensity of Point Pattern on Linear Network
}
\description{
  Computes the average number of points per unit length
  in a point pattern on a linear network.
}
\usage{
\method{intensity}{lpp}(X, ...)
}
\arguments{
  \item{X}{
    A point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic function \code{\link[spatstat.geom]{intensity}}
  It computes the empirical intensity of a point pattern
  on a linear network (object of class \code{"lpp"}),
  i.e. the average density of points per unit length.

  If the point pattern is multitype, the intensities of the
  different types are computed separately.
}
\value{
  A numeric value (giving the intensity) or numeric vector
  (giving the intensity for each possible type).
}
\seealso{
  \code{\link[spatstat.geom]{intensity}}, 
  \code{\link[spatstat.geom]{intensity.ppp}}
}
\examples{
  intensity(chicago)
}
\author{\adrian
  and \rolf
}
\keyword{spatial}
\keyword{nonparametric}
\concept{Linear network}
