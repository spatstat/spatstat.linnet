\name{linequad}
\alias{linequad}
\title{
  Quadrature Scheme on a Linear Network
}
\description{
  Generates a quadrature scheme (an object of class \code{"quad"})
  on a linear network.
}
\usage{
linequad(X, Y, \dots, eps = NULL, nd = 1000, random = FALSE)
}
\arguments{
  \item{X}{
    Data points. An object of class \code{"lpp"} or \code{"ppp"}.
  }
  \item{Y}{
    Line segments on which the points of \code{X} lie.
    An object of class \code{"psp"} or \code{"linnet"}.
    Required only when \code{X} is a \code{"ppp"} object.
  }
  \item{\dots}{
    Ignored.
  }
  \item{eps}{
    Optional. Spacing between successive dummy points along each
    segment. (This is the maximum spacing; some spacings will be shorter.)
  }
  \item{nd}{
    Optional. Total number of dummy locations to be generated.
    (Actual number may be larger.)
  }
  \item{random}{
    Logical value indicating whether the sequence of dummy points
    should start at a randomly-chosen position along each segment.
  }
}
\details{
  This command generates a quadrature scheme (object of class
  \code{"quad"}) from a pattern of points on a linear network.

  Normally the user does not need to call \code{linequad} explicitly.
  It is invoked by \pkg{spatstat} functions when needed.
  A quadrature scheme is required by \code{\link{lppm}}
  in order to fit point process models to point pattern data on a linear
  network. A quadrature scheme is also used by \code{\link{rhohat.lpp}}
  and other functions.

  In order to create the quadrature scheme, dummy points are placed
  along each line segment of the network. The dummy points are 
  evenly-spaced with spacing \code{eps}. The default is
  \code{eps = totlen/nd} where \code{totlen} is the total length of
  all line segments in the network.

  Every line segment of the network will contain at least one dummy
  location. Consequently the actual number of dummy location generated
  will typically be greater than \code{nd},
  especially when \code{nd} is small.
  If \code{eps} is specified, the
  number of dummy locations will be greater than \code{totlen/eps},
  especially when \code{eps} is large.

  If \code{X} is a multitype point pattern with \eqn{m} possible types,
  the dummy points will also be a marked point pattern.
  At each dummy location, \eqn{m} marked dummy points will be placed,
  one dummy point of each possible type.
  Additionally at each data location, a further \eqn{m - 1} dummy
  points will be placed, one dummy point of each possible type
  other than the type of the data point. The total number of dummy
  points will be \eqn{mk + (m-1) n = m(k+n) - n}
  and the total number of quadrature points will be \eqn{m(k+n)},
  where \eqn{k} is the number of dummy locations
  and \eqn{n} is the number of data points in \code{X}.
}
\value{
  A quadrature scheme (object of class \code{"quad"}).
}
\author{
  \adrian,
  Greg McSwiggan and Suman Rakshit.
}
\seealso{
  \code{\link{lppm}}
}
\keyword{datagen}
\keyword{spatial}
\concept{Linear network}
