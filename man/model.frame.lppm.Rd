\name{model.frame.lppm}
\alias{model.frame.lppm}
\title{
Extract the Variables in a Point Process Model on a Network
}
\description{
  Given a fitted point process model on a network, this function
  returns a data frame containing all the variables needed to
  fit the model using the Berman-Turner device.
}
\usage{
 \method{model.frame}{lppm}(formula, ...)
}
\arguments{
  \item{formula}{
    A fitted point process model on a linear network.
    An object of class \code{"lppm"}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{model.frame.glm}}.
  }
}
\details{
  The function \code{\link[stats]{model.frame}} is generic.
  This function is a method for \code{\link[stats]{model.frame}}
  for fitted point process models
  on a linear network (objects of class \code{"lppm"}).

  The first argument should be a fitted point process model;
  it has to be named \code{formula} for consistency with the generic
  function.

  The result is a data frame containing all the variables used in
  fitting the model. The data frame has one row for each quadrature point
  used in fitting the model. The quadrature scheme can be extracted using
  \code{\link{quad.ppm}}.  
}
\value{
  A \code{data.frame} containing all the variables used in the
  fitted model, plus additional variables specified in \code{\dots}.
  It has an additional attribute \code{"terms"} containing information
  about the model formula. For details see \code{\link{model.frame.glm}}.
}
\references{
  Baddeley, A. and Turner, R. (2000)
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42}, 283--322.
}
\seealso{
  \code{\link{lppm}},
  \code{\link[stats]{model.frame}},
  \code{\link[spatstat.model]{model.matrix.ppm}}
}
\examples{
  fit <- lppm(spiders ~ x)
  mf <- model.frame(fit)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{models}

\concept{Linear network}
