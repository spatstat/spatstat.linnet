\name{lurking.lppm}
\alias{lurking.lppm}               
\alias{lurking.lpp}               
\title{Lurking Variable Plot on a Linear Network}
\description{
  Plot point process residuals against a covariate,
  for a point process model on a linear network,
  or a point pattern on a linear network.
}
\usage{
\method{lurking}{lppm}(object, covariate,
                    type="raw",
                    \dots,
                    covname)
}
\arguments{
  \item{object}{
    The fitted point process model on a linear network
    (an object of class \code{"lppm"})
    for which diagnostics should be produced. This object
    is usually obtained from \code{\link[spatstat.linnet]{lppm}}.
    Alternatively, \code{object} may be a point pattern
    on a linear network (object of class \code{"lpp"}).
  }
  \item{covariate}{
    The covariate against which residuals should be plotted.
    Either a numeric vector, a pixel image, or an \code{expression}.
    See \emph{Details} below.
  }
  \item{type}{
    String indicating the type of residuals or weights to be computed.
    Choices include \code{"eem"},
    \code{"raw"}, \code{"inverse"} and \code{"pearson"}.
  }
  \item{\dots}{
    Further arguments passed to \code{\link[spatstat.model]{lurking.ppm}}
    to control the computation and the plot behaviour.
  }
  \item{covname}{
    A string name for the covariate, to be used in axis labels of plots.
  }
}
\value{
  The (invisible) return value is an object 
  belonging to the class \code{"lurk"}, for which there
  are methods for \code{plot} and \code{print}.

  This object is a list containing two dataframes
  \code{empirical} and \code{theoretical}. 
  The first dataframe \code{empirical} contains columns
  \code{covariate} and \code{value} giving the coordinates of the
  lurking variable plot. The second dataframe \code{theoretical}
  contains columns \code{covariate}, \code{mean} and \code{sd}
  giving the coordinates of the plot of the theoretical mean
  and standard deviation.
}
\details{
  This function generates a \sQuote{lurking variable} plot for a
  fitted point process model on a linear network.
  Residuals from the model represented by \code{object}
  are plotted against the covariate specified by \code{covariate}.
  This plot can be used to reveal departures from the fitted model,
  in particular, to reveal that the point pattern depends on the covariate.

  The function \code{lurking} is generic, with methods for
  \code{lppm} and \code{lpp} documented here.

  The argument \code{object} would usually be a fitted point process model
  on a network, obtained from the model-fitting algorithm
  \code{\link[spatstat.linnet]{lppm}}). 
  If \code{object} is a point pattern on a network
  (object of class \code{"lpp"}) then
  the model is taken to be the uniform Poisson process on the network
  (Complete Spatial Randomness on the Network).
  
  First the residuals from the fitted model (Baddeley et al, 2004)
  are computed at each pixel,
  or alternatively the `exponential energy marks' (Stoyan and Grabarnik,
  1991) are computed at each data point.
  The argument \code{type} selects the type of
  residual or weight. See \code{\link[spatstat.model]{diagnose.ppm}} for options
  and explanation.

  A lurking variable plot for point processes (Baddeley et al, 2004)
  displays either the cumulative sum of residuals/weights
  (if \code{cumulative = TRUE}) or a kernel-weighted average of the
  residuals/weights (if \code{cumulative = FALSE}) plotted against
  the covariate. The empirical plot (solid lines) is shown
  together with its expected value assuming the model is true
  (dashed lines) and optionally also the pointwise
  two-standard-deviation limits (grey shading).
  
  To be more precise, let \eqn{Z(u)} denote the value of the covariate
  at a spatial location \eqn{u} on the network.
  \itemize{
    \item
    If \code{cumulative=TRUE} then we plot \eqn{H(z)} against \eqn{z},
    where \eqn{H(z)} is the sum of the residuals 
    over all quadrature points where the covariate takes
    a value less than or equal to \eqn{z}, or the sum of the
    exponential energy weights over all data points where the covariate
    takes a value less than or equal to \eqn{z}.
    \item
    If \code{cumulative=FALSE} then we plot \eqn{h(z)} against \eqn{z},
    where \eqn{h(z)} is the derivative of \eqn{H(z)},
    computed approximately by spline smoothing.
  }
  For the point process residuals \eqn{E(H(z)) = 0},
  while for the exponential energy weights
  \eqn{E(H(z)) = } length of the subset of the network
  satisfying \eqn{Z(u) <= z}{Z(u) \le z}. 

  If the empirical and theoretical curves deviate substantially
  from one another, the interpretation is that the fitted model does
  not correctly account for dependence on the covariate.
  The correct form (of the spatial trend part of the model)
  may be suggested by the shape of the plot.
  
  If \code{plot.sd = TRUE}, then superimposed on the lurking variable
  plot are the pointwise
  two-standard-deviation error limits for \eqn{H(x)} calculated for the
  inhomogeneous Poisson process. The default is \code{plot.sd = TRUE}
  for Poisson models and \code{plot.sd = FALSE} for non-Poisson
  models.

  By default, the two-standard-deviation limits are calculated
  from the exact formula for the asymptotic variance
  of the residuals under the asymptotic normal approximation,
  equation (37) of Baddeley et al (2006).
  However, for compatibility with the original paper
  of Baddeley et al (2005), if \code{oldstyle=TRUE},
  the two-standard-deviation limits are calculated
  using the innovation variance, an over-estimate of the true
  variance of the residuals.

  The argument \code{covariate} is either a numeric vector, a pixel
  image, or an R language expression.
  If it is a numeric vector, it is assumed to contain
  the values of the covariate for each of the quadrature points
  in the fitted model. The quadrature points can be extracted by
  \code{\link[spatstat.model]{quad.ppm}(object)}.

  If \code{covariate} is a pixel image, it is assumed to contain the
  values of the covariate at each location in the window. The values of
  this image at the quadrature points will be extracted.

  Alternatively, if \code{covariate}
  is an \code{expression}, it will be evaluated in the same environment
  as the model formula used in fitting the model \code{object}. It must
  yield a vector of the same length as the number of quadrature points.
  The expression may contain the terms \code{x} and \code{y} representing the
  cartesian coordinates, and may also contain other variables that were
  available when the model was fitted. Certain variable names are
  reserved words; see \code{\link[spatstat.linnet]{lppm}}.
}
\seealso{
 \code{\link[spatstat.linnet]{lppm}}
}
\references{
  Baddeley, A., Turner, R., \Moller, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Baddeley, A., \Moller, J. and Pakes, A.G. (2006)
  Properties of residuals for spatial point processes.
  \emph{Annals of the Institute of Statistical Mathematics}
  \bold{60}, 627--649.
  
  Stoyan, D. and Grabarnik, P. (1991)
  Second-order characteristics for stochastic structures connected with
  Gibbs point processes.
  \emph{Mathematische Nachrichten}, 151:95--100.
}
\examples{
  fit <- lppm(spiders ~ y)
  (b <- lurking(fit, expression(x), type="raw"))
  lurking(fit, expression(x), type="raw", cumulative=FALSE)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{models}
\keyword{hplot}
\concept{Linear network}

