\name{spatstat.linnet-deprecated}
\alias{circumradius.linnet}
\alias{rjitterlpp}
\alias{bw.relrisklpp}
\title{Deprecated spatstat.linnet functions}
\description{
  Deprecated spatstat.linnet functions.
}
\usage{
\method{circumradius}{linnet}(x, \dots)
rjitterlpp(X, \dots)
bw.relrisklpp(X, \dots)
}

\details{
  These functions are deprecated, and will eventually be deleted from
  the \pkg{spatstat.linnet} package.
  

  
  \code{rjitterlpp} is replaced by \code{rjitter.lpp}, a method for the
  generic \code{rjitter}.

  \code{circumradius.linnet} is replaced by the more appropriately named
  \code{boundingradius.linnet}.

  \code{bw.relrisklpp} is replaced by the method \code{bw.relrisk.lpp}
  for the generic \code{bw.relrisk}.
}
\value{
\code{circumradius.linnet} returns a numeric value.

\code{rjitterlpp} returns a point pattern on a linear network (object of
class \code{"lpp"}).

\code{bw.relrisklpp} returns an object of class \code{"bw.optim"} which
is a numeric value with additional attributes.
}
\keyword{internal}
