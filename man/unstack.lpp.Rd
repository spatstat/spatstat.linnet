\name{unstack.lpp}
\alias{unstack.lpp}
\alias{unstack.lintess}
\title{
  Separate Multiple Columns of Marks 
}
\description{
  Given a spatial pattern on a network, with several columns of marks,
  take one column at a time, and return a list of spatial patterns
  each having only one column of marks.
}
\usage{
\method{unstack}{lpp}(x, \dots)

\method{unstack}{lintess}(x, \dots)
}
\arguments{
  \item{x}{
    A spatial point pattern
    (object of class \code{"lpp"})
    or a tessellation on a linear network (object of class \code{"lintess"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  The functions defined here are methods for the generic
  \code{\link[utils]{unstack}}. The functions expect a spatial object
  \code{x} which has several columns of marks; they separate the columns,
  and return a list of spatial objects, each having only one column of marks.
  
  If \code{x} has several columns of marks (i.e. \code{marks(x)} is a
  matrix, data frame or hyperframe with several columns),
  then \code{y <- unstack(x)} is a list of spatial objects, each of the same
  kind as \code{x}. The \code{j}th entry \code{y[[j]]} is equivalent to
  \code{x} except that it only includes
  the \code{j}th column of \code{marks(x)}.
  
  If \code{x} has no marks, or has only a single column of marks,
  the result is a list consisting of one entry, which is \code{x}.
}
\value{
  A list, of class \code{"solist"}, whose entries are objects of the
  same type as \code{x}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[utils]{unstack}}

  \code{\link[spatstat.geom]{unstack.ppp}},
  \code{\link[spatstat.core]{unstack.msr}}.

  See also methods for the generic \code{\link[base]{split}} such as
  \code{\link[spatstat.geom]{split.ppx}} which applies to \code{"lpp"} objects.
}
\examples{
   X <- runiflpp(5, simplenet)
   marks(X) <- data.frame(id=1:5, code=factor(letters[1:5]))
   unstack(X)
}
\keyword{spatial}
\keyword{manip}
