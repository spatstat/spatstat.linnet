\name{residuals.lppm}
\alias{residuals.lppm}
\title{
  Residuals for Fitted Point Process Model on a Network
}
\description{
  Given a point process model fitted to a point pattern on a linear network,
  compute residuals of the fitted model.
}
\usage{
  \method{residuals}{lppm}(object, type="raw", \dots)
}
\arguments{
  \item{object}{
    The fitted point process model (an object of class \code{"ppm"})
    for which residuals should be calculated.
  }
  \item{type}{
    String indicating the type of residuals to be calculated.
    Current options are
    \code{"raw"}, \code{"inverse"}, \code{"pearson"} and \code{"score"}.
    A partial match is adequate.
  }
  \item{\dots}{
    Other arguments are currently ignored.
  }
}
\value{
  An object of class \code{"msr"} 
  representing a signed measure or vector-valued measure
  (see \code{\link[spatstat.model]{msr}}). This object can be plotted.
}
\details{
  This function computes several kinds of residuals for the fit of
  a point process model to a spatial point pattern on a linear network.
  It is an extension of the method of Baddeley et al (2005)
  to point process models on a network.
  Use \code{\link[spatstat.model]{plot.msr}} to plot the residuals directly.

  The argument \code{object} must be a fitted point process model
  on a network (object of class \code{"lppm"}).
  Such objects are produced by the model-fitting algorithm \code{\link{lppm}}.
  This fitted model object contains complete
  information about the original data pattern. 

  Residuals are attached both to the data points and to some
  other points in the window of observation (namely, to the dummy
  points of the quadrature scheme used to fit the model).
  If the fitted model is correct, then the sum of the 
  residuals over all (data and dummy) points in a spatial region \eqn{B}
  has mean zero. For further explanation, see Baddeley et al (2005).
  
  The type of residual
  is chosen by the argument \code{type}. Current options are

  \describe{
    \item{\code{"raw"}:}{
      the raw residuals
      \deqn{
	r_j = z_j - w_j \lambda_j
      }{
	r[j] = z[j] - w[j] lambda[j]
      }
      at the quadrature points \eqn{u_j}{u[j]},
      where \eqn{z_j}{z[j]} is the indicator equal to 1 if \eqn{u_j}{u[j]}
      is a data point and 0 if \eqn{u_j}{u[j]} is a dummy point;
      \eqn{w_j}{w[j]} is the quadrature weight attached to
      \eqn{u_j}{u[j]}; and
      \deqn{\lambda_j = \hat\lambda(u_j,x)}{lambda[j] = lambda(u[j],x)}
      is the conditional intensity of the fitted model at \eqn{u_j}{u[j]}.
      These are the spatial analogue of the martingale residuals
      of a one-dimensional counting process.
    }
    \item{\code{"inverse"}:}{
      the `inverse-lambda' residuals (Baddeley et al, 2005)
      \deqn{
	r^{(I)}_j = \frac{r_j}{\lambda_j}
        = \frac{z_j}{\lambda_j} - w_j
      }{
	rI[j] = r[j]/lambda[j] = z[j]/lambda[j] - w[j]
      }
      obtained by dividing the raw residuals by 
      the fitted conditional intensity. These are 
      a counterpart of the exponential energy marks
      (see \code{\link[spatstat.model]{eem}}).
    }
    \item{\code{"pearson"}:}{
      the Pearson residuals (Baddeley et al, 2005)
      \deqn{
	r^{(P)}_j = \frac{r_j}{\sqrt{\lambda_j}}
        = \frac{z_j}{\sqrt{\lambda_j}}
        - w_j \sqrt{\lambda_j}
      }{
	rP[j] = r[j]/sqrt(lambda[j])
	= z[j]/sqrt(lambda[j]) - w[j] sqrt(lambda[j])
      }
      obtained by dividing the raw residuals by the
      square root of the fitted conditional intensity.
      The Pearson residuals are standardised, in the sense
      that if the model (true and fitted) is Poisson,
      then the sum of the Pearson residuals in a spatial region \eqn{B}
      has variance equal to the area of \eqn{B}.
    }
    \item{\code{"score"}:}{
      the score residuals (Baddeley et al, 2005)
      \deqn{
	r_j = (z_j - w_j \lambda_j) x_j
      }{
	r[j] = (z[j] - w[j] lambda[j]) * x[j,]
      }
      obtained by multiplying the raw residuals \eqn{r_j}{r[j]}
      by the covariates \eqn{x_j}{x[j,]} for quadrature point \eqn{j}.
      The score residuals always sum to zero.
    }
  }
  The result of \code{residuals.ppm} is a measure
  (object of class \code{"msr"}).
  Use \code{\link[spatstat.model]{plot.msr}} to plot the residuals directly.
  Use \code{\link[spatstat.model]{integral.msr}} to compute the total residual.
}
\references{
  Baddeley, A., Turner, R., \Moller, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.
}
\seealso{
 \code{\link{lppm}},
 \code{\link[spatstat.model]{msr}},
 \code{\link[spatstat.model]{residuals.ppm}}.
}
\examples{
   fit <- lppm(unmark(chicago) ~ x + y)

   # raw residuals
   rr <- residuals(fit)
   rr

   # Pearson residuals
   rp <- residuals(fit, type="pe")
   rp
   plot(rp, main="Pearson residuals")

   ## multitype data
   fitm <- lppm(chicago ~ (x+y) * marks, eps=100)
   rpm <- residuals(fitm, type="pe")
   ## plot(rpm) would display 7 panels, one for each crime type
   ## Select residuals for crime type = Assault
   plot(split(rpm)[["assault"]], markscale=2)
}
\author{\spatstatAuthors.}
\keyword{spatial}
\keyword{models}
\keyword{methods}
