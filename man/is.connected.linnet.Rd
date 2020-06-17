\name{is.connected.linnet}  
\Rdversion{1.1}
\alias{is.connected.linnet}
\title{
  Determine Whether a Linear Network is Connected
}
\description{
  Determine whether a linear network is 
  topologically connected.
}
\usage{
\method{is.connected}{linnet}(X, \dots)
}
\arguments{
  \item{X}{
    A linear network (object of class \code{"linnet"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{connected.linnet}} to determine the
    connected components.
  }
}
\details{
  The command \code{is.connected(X)} returns \code{TRUE} if the network
  \code{X} consists of a single, topologically-connected piece,
  and returns \code{FALSE} if \code{X} consists of several pieces
  which are not joined together.
  
  The function \code{\link[spatstat.core]{is.connected}} is generic,
  with methods for several classes.
  This help file documents the 
  method for linear networks, \code{is.connected.linnet}.
}
\value{
  A logical value.
}
\seealso{
  \code{\link[spatstat.core]{is.connected}},
  \code{\link[spatstat.core]{connected}},
  \code{\link{connected.lpp}}.
}
\examples{
  is.connected(simplenet)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{math}
