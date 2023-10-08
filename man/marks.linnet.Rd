\name{marks.linnet}
\alias{marks.linnet}
\alias{marks<-.linnet}
\title{Marks of a Network}
\description{
  Extract or change the marks attached to
  vertices or segments of a linear network.
}
\usage{
\method{marks}{linnet}(x, of=c("segments", "vertices"), \dots)           

\method{marks}{linnet}(x, of=c("segments", "vertices"), \dots) <- value
}
\arguments{
  \item{x}{
    Linear network (object of class \code{"linnet"}).
  }
  \item{of}{
    Character string (partially matched)
    specifying whether the marks are attached
    to the vertices of the network (\code{of="vertices"})
    or to the line segments of the network (\code{of="segments"}, the default).
  }
  \item{\dots}{
    Ignored.
  }
  \item{value}{
    Vector or data frame of mark values,
    or \code{NULL}.
  }
}
\value{
  For \code{marks(x)}, the result is a vector, factor or data frame,
  containing the mark values attached to the vertices or the
  segments of \code{x}.
  If there are no marks, the result is \code{NULL}.

  For \code{marks(x) <- value}, the result is the updated network
  \code{x} (with the side-effect that the dataset \code{x} is updated in
  the current environment).
}
\details{
  These functions extract or change the marks
  attached to the network \code{x}.
  They are methods for the generic functions
  \code{\link[spatstat.geom]{marks}} and 
  \code{\link[spatstat.geom]{marks<-}}
  for the class \code{"linnet"} of linear networks.

  A linear network may include a set of marks attached to the line segments,
  and a separate set of marks attached to the vertices.
  Each set of marks can be a vector, a factor, or a data frame.

  The expression \code{marks(x, of)} extracts the marks from \code{x}.
  The assignment \code{marks(x, of) <- value} assigns new marks to the
  dataset \code{x}, and updates the dataset \code{x} in the current
  environment. 
  The argument \code{of} specifies whether we are referring to
  the segments or the vertices.
  
  For the assignment \code{marks(x, "segments") <- value}, the \code{value}
  should be a vector or factor of length equal to the number of
  segments in \code{x}, or a data frame with as many rows
  as there are segments in \code{x}. 
  If \code{value} is a single value,
  or a data frame with one row, then it will be replicated
  so that the same marks will be attached to each segment.
  Similarly for \code{marks(x, "vertices") <- value} the number of
  marks must match the number of vertices.
  
  To remove marks, use \code{marks(x, of) <- NULL}.

  To extract the vertices (including their marks) as a point pattern,
  use \code{vertices(x)}. To extract the segments (including their
  marks) as a line segment pattern, use \code{as.psp(x)}.
}
\seealso{
  \code{\link{linnet}},
  \code{\link[spatstat.geom]{marks}},
  \code{\link[spatstat.geom]{marks<-}}
}
\examples{
  L <- simplenet
  marks(L, "vertices") <- letters[1:nvertices(L)]
  marks(L, "segments") <- runif(nsegments(L))
  L
  marks(L, "v")
  marks(L, "s")
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
\concept{Linear network}
