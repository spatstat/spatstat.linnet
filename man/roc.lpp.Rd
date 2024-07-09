\name{roc.lpp}
\alias{roc.lpp}
\alias{roc.lppm}
\title{
  Receiver Operating Characteristic for Data on a Network
}
\description{
  Computes the Receiver Operating Characteristic curve
  for a point pattern or a fitted point process model
  on a linear network.
}
\usage{
\method{roc}{lpp}(X, covariate, \dots, high = TRUE)

\method{roc}{lppm}(X, \dots)
}
\arguments{
  \item{X}{
    Point pattern on a network (object of class \code{"lpp"})
    or fitted point process model on a network
    (object of class \code{"lppm"}).
  }
  \item{covariate}{
    Spatial covariate. Either a \code{function(x,y)},
    a pixel image (object of class \code{"im"} or \code{"linim"}), or
    one of the strings \code{"x"} or \code{"y"} indicating the
    Cartesian coordinates.    
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}} controlling the
    pixel resolution for calculations.
  }
  \item{high}{
    Logical value indicating whether the threshold operation
    should favour high or low values of the covariate.
  }
}
\details{
  The command \code{roc}
  computes the Receiver Operating
  Characteristic curve.
  The area under the ROC is computed by \code{\link[spatstat.explore]{auc}}.

  The function \code{\link[spatstat.explore]{roc}} is generic,
  with methods for \code{"ppp"} and \code{"ppm"} described in the help file
  for \code{\link[spatstat.explore]{roc}}.

  This help file describes the methods for classes \code{"lpp"} and
  \code{"lppm"}.

  For a point pattern \code{X} and a covariate \code{Z}, the
  ROC is a plot showing the ability of the 
  covariate to separate the spatial domain
  into areas of high and low density of points.
  For each possible threshold \eqn{z}, the algorithm calculates
  the fraction \eqn{a(z)} of area in the study region where the
  covariate takes a value greater than \eqn{z}, and the
  fraction \eqn{b(z)} of data points for which the covariate value
  is greater than \eqn{z}. The ROC is a plot of \eqn{b(z)} against
  \eqn{a(z)} for all thresholds \eqn{z}. 
  
  For a fitted point process model, 
  the ROC shows the ability of the
  fitted model intensity to separate the spatial domain
  into areas of high and low density of points.
  The ROC is \bold{not} a diagnostic for the goodness-of-fit of the model
  (Lobo et al, 2007).
}
\value{
  Function value table (object of class \code{"fv"})
  which can be plotted to show the ROC curve.
}
\references{
  Lobo, J.M.,
  \ifelse{latex}{\out{Jim{\'e}nez}}{Jimenez}-Valverde, A.
  and Real, R. (2007)
  AUC: a misleading measure of the performance of predictive
  distribution models.
  \emph{Global Ecology and Biogeography} \bold{17}(2) 145--151.

  Nam, B.-H. and D'Agostino, R. (2002)
  Discrimination index, the area under the {ROC} curve.
  Pages 267--279 in 
  Huber-Carol, C., Balakrishnan, N., Nikulin, M.S. 
  and Mesbah, M., \emph{Goodness-of-fit tests and model validity},
  \ifelse{latex}{\out{Birkh{\"a}user}}{Birkhauser}, Basel.
}
\author{
\spatstatAuthors.
}
\seealso{
\code{\link{auc.lpp}}
}
\examples{
  plot(roc(spiders, "x"))
  fit <- lppm(spiders ~ x)
  plot(roc(fit))
}
\keyword{spatial}
\concept{Linear network}
