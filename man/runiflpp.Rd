\name{runiflpp}
\alias{runiflpp}
\title{
  Uniform Random Points on a Linear Network
}
\description{
  Generates \eqn{n} random points, independently and
  uniformly distributed, on a linear network.
}
\usage{
runiflpp(n, L, nsim=1, drop=TRUE, ex=NULL)
}
\arguments{
  \item{n}{
    Number of random points to generate.
    A nonnegative integer, or a vector of integers
    specifying the number of points of each type.
  }
  \item{L}{
    A linear network (object of class \code{"linnet"},
    see \code{\link{linnet}}).
  }
  \item{nsim}{Number of simulated realisations to generate.}
  \item{drop}{
    Logical value indicating what to do when \code{nsim=1}.
    If \code{drop=TRUE} (the default), the result is a point pattern.
    If \code{drop=FALSE}, the result is a list with one entry which is a
    point pattern.
  }
  \item{ex}{
    Optional. A point pattern on a network
    (object of class \code{"lpp"}) which serves as an example
    to determine the default values of \code{n} and \code{L}.
    See Details.
  }
}
\details{
  The specified number \code{n} of random points is
  generated with uniform distribution on the network \code{L}.
  The random points are generated using \code{\link[spatstat.random]{runifpointOnLines}}.

  If \code{n} is an integer vector, then a multitype point pattern is
  generated, with \code{n[i]} random points of type \code{i}.
  
  If \code{ex} is given, then it serves as an example for determining
  \code{n} and \code{L}. The default value of \code{n}
  will be the number of points in \code{ex}
  (or the number of points of each type in \code{ex} if it is
  multitype).
  The default value of \code{L} will be 
  the network on which \code{ex} is defined. 
}
\value{
  If \code{nsim = 1} and \code{drop=TRUE},
  a point pattern on a linear network
  (object of class \code{"lpp"}).
  Otherwise, a list of such point patterns.
}
\author{
  \wei and 
  \adrian 
}
\seealso{
 \code{\link{rlpp}} for non-uniform random points;
 \code{\link{rpoislpp}} for Poisson point process;
 
 \code{\link{lpp}},
 \code{\link{linnet}}
}
\examples{
   X <- runiflpp(10, simplenet)
   plot(X)
  # marked
   Z <- runiflpp(c(a=10, b=3), simplenet)
  # using 'ex' 
   U <- runiflpp(ex=Z)
}
\keyword{spatial}
\keyword{datagen}
\concept{Linear network}
