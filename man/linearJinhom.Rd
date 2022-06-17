\name{linearJinhom}
\alias{linearJinhom}
\title{
  Inhomogeneous Linear J-function for Point Processes on Linear Networks}
\description{
  Computes an estimate of the inhomogeneous linear \eqn{J}-function
  for a point pattern on a linear network.
}
\usage{
linearJinhom(X, lambda = NULL, lmin=NULL,
             ...,
             r=NULL, rmax=NULL,
             distance=c("path","euclidean"),
             densitymethod=c("kernel", "Voronoi"),
             sigma=bw.scott.iso,
             f=0.2, nrep=200, ngrid=256) 
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{lambda}{
    Intensity values for the point pattern. Either a numeric vector,
    a \code{function}, a pixel image
    (object of class \code{"im"} or \code{"linim"}) or
    a fitted point process model (object of class \code{"ppm"}
    or \code{"lppm"}).
  }
  \item{lmin}{
    Optional. The minimum possible value of the intensity
    over the network. A positive numerical value.
  }
  \item{r}{
    Optional. Numeric vector of values of the function argument \eqn{r}.
    There is a sensible default.
  }
  \item{rmax}{
    Optional. Numeric value specifying the largest desired value of \eqn{r}.
    There is a sensible default.
  }
  \item{distance}{
    A string (partially matched) specifying the
    metric that will be used to measure distances between
    points on the network: \code{distance="path"} is the shortest-path
    distance, and \code{distance="euclidean"} is the Euclidean distance.
  }
  \item{densitymethod}{
    String (partially matched) specifying the method that will be used to
    estimate the intensity \code{lambda}, if \code{lambda} is not given:
    \code{densitymethod="kernel"} specifies kernel smoothing
    and \code{densitymethod="Voronoi"} specifies Voronoi estimation.
    See Details.
  }
  \item{sigma}{
    Smoothing bandwidth used to estimate \code{lambda} by kernel smoothing,
    if \code{lambda} is not given and \code{densitymethod="kernel"}.
    Either a numeric value, or a function that can be applied to
    \code{X} to compute the bandwidth.
  }
  \item{f,nrep}{
    Arguments passed to the algorithm for estimating the intensity
    by Voronoi estimation, if \code{lambda} is not given
    and \code{densitymethod="Voronoi"}. 
  }
  \item{\dots}{
    Additional arguments passed to the algorithms that estimate
    the intensity, if \code{lambda} is not given.
  }
  \item{ngrid}{
    Integer specifying the number of sample points on the network that
    will be used to estimate the inhomogeneous empty space function
    \eqn{F}.
  }
}
\details{
  This function computes the geometrically corrected inhomogeneous
  linear \eqn{J}-function for point processes on linear networks
  defined by Cronie et al (2020).

  The argument \code{lambda} is the (estimated) intensity of the
  underlying point process. It should be either a numeric vector
  (giving intensity values at the points of \code{X}), 
  a \code{function}, a pixel image
  (object of class \code{"im"} or \code{"linim"}) or
  a fitted point process model (object of class \code{"ppm"}
  or \code{"lppm"}).

  If \code{lambda} is not given, it will be estimated from the observed
  point pattern \code{X} as follows:
  \itemize{
    \item If \code{densitymethod="kernel"}, the intensity will be
    estimated by kernel smoothing, using the fast estimator
    \code{\link[spatstat.linnet]{densityQuick.lpp}}
    introduced by Rakshit et al (2019). The smoothing bandwidth
    \code{sigma} is required. It may be specified as a numeric value,
    or as a function that can be applied to \code{X} to obtain a
    bandwidth value. Examples of the latter include
    \code{\link[spatstat.explore]{bw.scott.iso}}
    and 
    \code{\link[spatstat.linnet]{bw.lppl}}.
    Additional arguments \code{\dots} will be passed to \code{sigma}
    and to \code{\link[spatstat.linnet]{densityQuick.lpp}}.
    \item If \code{densitymethod = "Voronoi"}, the intensity will be
    estimated using the resample-smoothed Voronoi estimator
    \code{\link[spatstat.linnet]{densityVoronoi.lpp}}
    introduced by Moradi et al (2019). The arguments \code{f} and \code{nrep} 
    are passed to \code{\link[spatstat.linnet]{densityVoronoi.lpp}}
    and determine the retention probability and the number of
    replicates, respectively. Additional arguments \code{\dots} will be
    passed to \code{\link[spatstat.linnet]{densityVoronoi.lpp}}.
  }
}
\author{
  \mehdi and \adrian.
}
\value{
  Function value table (object of class \code{"fv"}).
}
\examples{
  if(interactive()) {
    plot(linearJinhom(spiders))
  } else {
    bottomhalf <- owin(c(0, 1125), c(0, 500))
    plot(linearJinhom(spiders[bottomhalf]))
  }
}
\references{
Cronie, O., Moradi, M., and Mateu, J. (2020) Inhomogeneous higher-order
summary statistics for point processes on linear
networks. \emph{Statistics and Computing} \bold{30} (6) 1221--1239.

Moradi, M., Cronie, 0., Rubak, E., Lachieze-Rey, R.,
Mateu, J. and Baddeley, A. (2019)
Resample-smoothing of Voronoi intensity estimators.
\emph{Statistics and Computing} \bold{29} (5) 995--1010.

Rakshit, S., Davies, T., Moradi, M., 
McSwiggan, G., Nair, G., Mateu, J. and Baddeley, A. (2019)
Fast kernel smoothing of point patterns on a large
network using 2D convolution.
\emph{International Statistical Review} \bold{87} (3) 531--556.
DOI: 10.1111/insr.12327.
}
\seealso{
  \code{\link[spatstat.explore]{bw.scott.iso}},
  \code{\link[spatstat.linnet]{bw.lppl}},
  \code{\link[spatstat.linnet]{densityVoronoi.lpp}},
  \code{\link[spatstat.linnet]{densityQuick.lpp}}

  \code{\link[spatstat.linnet]{linearKinhom}}

  \code{\link[spatstat.explore]{Jinhom}}
}
\keyword{spatial}
\keyword{nonparametric}
