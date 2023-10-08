\name{is.multitype.lpp}
\alias{is.multitype.lpp}
\title{Test Whether A Point Pattern on a Network is Multitype}
\description{
  Tests whether a point pattern on a network has ``marks''
  attached to the points which classify the points into several types.
}
\usage{
  \method{is.multitype}{lpp}(X, na.action="warn", \dots) 
}
\arguments{
  \item{X}{
    Point pattern on a linear networl
    (object of class \code{"lpp"}).
  }
  \item{na.action}{
    String indicating what to do if \code{NA} values are
    encountered amongst the marks.
    Options are \code{"warn"}, \code{"fatal"} and \code{"ignore"}.
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  Logical value, equal to \code{TRUE} if \code{X} is a multitype point pattern.
}
\details{
  ``Marks'' are observations attached to each point of a point pattern.
  For example the \code{\link[spatstat.data]{chicago}} dataset contains
  the locations of crimes, each crime location
  being marked by the type of crime.

  This function tests whether the point pattern \code{X}
  contains or involves marked points, \bold{and} that the
  marks are a factor.
  It is a method for the generic function
  \code{\link[spatstat.geom]{is.multitype}}.
  
  The argument \code{na.action} determines what action will be taken
  if the point pattern has a vector of marks but some or all of the
  marks are \code{NA}. Options are   \code{"fatal"} to cause a fatal
  error; \code{"warn"} to issue a warning and then return \code{TRUE};
  and \code{"ignore"} to take no action except returning \code{TRUE}.
}
\seealso{
  \code{\link[spatstat.geom]{is.multitype}},
  \code{\link{is.multitype.lppm}}
}
\examples{
   is.multitype(chicago)
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{manip}
\concept{Linear network}
