\name{terminalvertices}
\alias{terminalvertices}
\title{Terminal Vertices of a Linear Network}
\description{
  Finds the terminal vertices of a linear network.
}
\usage{
 terminalvertices(L)
}
\arguments{
  \item{L}{A linear network (object of class \code{"linnet"}).}
}
\value{
  A point pattern on the same linear network
  (object of class \code{"lpp"}).
}
\details{
  Given the linear network \code{L}, this function
  examines the vertices (segment endpoints) of the network
  and determines which of them are \sQuote{terminal} vertices
  (i.e. the endpoint of only one segment). These terminal vertices
  are returned as a point pattern on the network.
}
\seealso{
  \code{\link{vertices.linnet}}.
}
\examples{
  B <- terminalvertices(simplenet)
  plot(simplenet, main="")
  plot(B, add=TRUE, pch=16, cex=2)
}
\author{
  \spatstatAuthors and \mehdi.
}
\keyword{spatial}
\keyword{graphs}
