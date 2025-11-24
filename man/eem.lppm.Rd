\name{eem.lppm}
\alias{eem.lppm}
\title{
  Exponential Energy Marks on a Linear Network
}
\description{
  Given a point process model fitted to a point pattern on a linear network,
  compute the Stoyan-Grabarnik diagnostic ``exponential energy marks''
  for the data points. 
}
\usage{
  \method{eem}{lppm}(fit, \dots)
}
\arguments{
  \item{fit}{
    The fitted point process model. An object of class \code{"lppm"}.
  }
  \item{\dots}{Ignored.}
}
\value{
  A vector containing the values of the exponential energy mark
  for each point in the pattern.
}
\details{
  Stoyan and Grabarnik (1991) proposed a diagnostic
  tool for point process models fitted to spatial point pattern data.
  Each point \eqn{x_i}{x[i]} of the data pattern \eqn{X}
  is given a `mark' or `weight'
  \deqn{m_i = \frac 1 {\hat\lambda(x_i,X)}}{m[i] = 1/\lambda(x[i],X)}
  where \eqn{\hat\lambda(x_i,X)}{\lambda(x[i],X)}
  is the conditional intensity of the fitted model.
  If the fitted model is correct, then the sum of these marks
  for all points in a region \eqn{B} has expected value equal to the
  area of \eqn{B}.

  The function \code{\link[spatstat.model]{eem}} is generic, with methods
  for various classes of models. This page documents the method
  \code{eem.lppm} for the class \code{"lppm"}.
  
  The argument \code{fit} must be a fitted point process model
  on a linear network (object of class \code{"lppm"}).
  Such objects are produced by the fitting algorithm \code{\link{lppm}}).
  This fitted model object contains complete
  information about the original data pattern and the model that was
  fitted to it.

  The value returned by \code{eem} is the vector
  of weights \eqn{m[i]}{m_i} associated with the points \eqn{x[i]}{x_i}
  of the original data pattern. The original data pattern
  (in corresponding order) can be
  extracted from \code{fit} using \code{\link{response}}.
  
  The function \code{\link{diagnose.lppm}}
  produces a set of sensible diagnostic plots based on these weights.
}
\references{
  Stoyan, D. and Grabarnik, P. (1991)
  Second-order characteristics for stochastic structures connected with
  Gibbs point processes.
  \emph{Mathematische Nachrichten}, 151:95--100.
}
\seealso{
 \code{\link[spatstat.model]{diagnose.ppm}},
 \code{\link{residuals.lppm}},
 \code{\link{lppm}}
}
\examples{
    fit <- lppm(spiders ~ x + y)
    ee <- eem(fit)
    sum(ee)/volume(domain(spiders)) # should be about 1 if model is correct
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{models}
