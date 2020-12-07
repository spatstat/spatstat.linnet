\name{is.multitype.lppm}
\alias{is.multitype.lppm}
\title{Test Whether A Point Process Model is Multitype}
\description{
  Tests whether a fitted point process model on a network involves ``marks''
  attached to the points that classify the points into several types.
}
\usage{
  \method{is.multitype}{lppm}(X, \dots) 
}
\arguments{
  \item{X}{
    Fitted point process model on a linear network
    (object of class \code{"lppm"})
    usually obtained from \code{\link{lppm}}.
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  Logical value, equal to \code{TRUE} if
  \code{X} is a model that was fitted to a multitype point pattern dataset.
}
\details{
  ``Marks'' are observations attached to each point of a point pattern.
  For example the \code{\link[spatstat.data]{chicago}} dataset contains
  the locations of crimes, each crime location
  being marked by the type of crime.

  The argument \code{X} is a fitted point process model
  on a network (an object of class \code{"lppm"}) typically obtained
  by fitting a model to point pattern data on a network
  using \code{\link{lppm}}.

  This function returns \code{TRUE} if the \emph{original data}
  (to which the model \code{X} was fitted) were a multitype point pattern.

  Note that this is not the same as testing whether the
  model involves terms that depend on the marks (i.e. whether the
  fitted model ignores the marks in the data).
  See the Examples for a trick for doing this.

  If this function returns \code{TRUE}, the implications are
  (for example) that
  any simulation of this model will require simulation of random marks
  as well as random point locations.
}
\seealso{
  \code{\link[spatstat.geom]{is.multitype}},
  \code{\link{is.multitype.lpp}}
}
\examples{
  fit <- lppm(chicago ~ x)
  is.multitype(fit)
  # TRUE because chicago data are multitype

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
