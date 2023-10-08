\name{rThomaslpp}
\alias{rThomaslpp}
\title{
  Simulate Thomas Process on Linear Network
}
\description{
  Generate a random point pattern, a realisation of the
  Thomas cluster process, on a linear network.
}
\usage{
  rThomaslpp(kappa, scale, mu, L, \dots, nsim=1, drop=TRUE)
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres.
    A single positive number, a \code{function(x,y)}, or a pixel image
    (object of class \code{"im"} or \code{"linim"}).
  }
  \item{scale}{
    Standard deviation of random displacement (along the network)
    of a point from its cluster centre.
  }
  \item{mu}{
    Mean number of points per cluster (a single positive number)
    or reference intensity for the cluster points (a function or
    a pixel image).
  }
  \item{L}{
    Linear network (object of class \code{"linnet"})
    on which the point pattern should be generated.
  }
  \item{\dots}{
    Arguments passed to \code{\link{rpoisppOnLines}}.
  }
  \item{nsim}{
    Number of simulated realisations to generate.
  }
  \item{drop}{
    Logical value indicating what to do when \code{nsim=1}.
    If \code{drop=TRUE} (the default), the result is a point pattern.
    If \code{drop=FALSE}, the result is a list with one entry which is a
    point pattern.
  }
}
\details{
  This function generates realisations of the Thomas cluster process
  on a linear network, described by Baddeley et al (2017).
  
  Argument \code{L} can be omitted, and defaults to \code{as.linnet(kappa)},
  when \code{kappa} is a function on a linear network (class
  \code{"linfun"}) or a pixel image on a linear network (\code{"linim"}).
}
\value{
  A point pattern on a network (object of class \code{"lpp"})
  or a list of point patterns on the network.
}
\references{
  Baddeley, A., Nair, G., Rakshit, S. and McSwiggan, G. (2017)
  `Stationary' point processes are uncommon on linear networks.
  \emph{STAT} \bold{6} (1) 68--78.

  Baddeley, A., Nair, G., Rakshit, S., McSwiggan, G. and Davies, T.M. (2021)
  Analysing point patterns on networks --- a review.
  \emph{Spatial Statistics} \bold{42}, 100435,
  DOI 10.1016/j.spasta.2020.100435.
}
\author{
  Greg McSwiggan and \adrian.
}
\seealso{
  \code{\link{rpoislpp}}
}
\examples{
  plot(rThomaslpp(4, 0.07, 5, simplenet))
}
\keyword{datagen}
\keyword{spatial}
\concept{Linear network}
