\name{is.marked.lppm}
\alias{is.marked.lppm}
\title{Test Whether A Point Process Model is Marked}
\description{
  Tests whether a fitted point process model on a network involves ``marks''
  attached to the points.
}
\usage{
  \method{is.marked}{lppm}(X, \dots) 
}
\arguments{
  \item{X}{
    Fitted point process model on a linear networ
    (object of class \code{"lppm"})
    usually obtained from \code{\link{lppm}}.
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  Logical value, equal to \code{TRUE} if
  \code{X} is a model that was fitted to a marked point pattern dataset.
}
\details{
  ``Marks'' are observations attached to each point of a point pattern.
  For example the \code{\link[spatstat.data]{chicago}} dataset contains
  the locations of crimes, each crime location
  being marked by the type of crime.

  The argument \code{X} is a fitted point process model on a network
  (an object of class \code{"lppm"}) typically obtained
  by fitting a model to point pattern data using \code{\link{lppm}}.

  This function returns \code{TRUE} if the \emph{original data}
  (to which the model \code{X} was fitted) were a marked point pattern.

  Note that this is not the same as testing whether the
  model involves terms that depend on the marks (i.e. whether the
  fitted model ignores the marks in the data).
  See the Examples for a trick to do this.

  If this function returns \code{TRUE}, the implications are
  (for example) that
  any simulation of this model will require simulation of random marks
  as well as random point locations.
}
\seealso{
  \code{\link[spatstat.geom]{is.marked}}.
}
\examples{
  fit <- lppm(chicago ~ x)
  is.marked(fit)
  ## result is TRUE, i.e. the data are marked

  ## To check whether the model involves marks:
  "marks" \%in\% spatstat.utils::variablesinformula(formula(fit))
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{manip}
\keyword{models}
