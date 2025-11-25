\name{qqplot.lppm}
\alias{qqplot.lppm}
\title{
  Q-Q Plot of Residuals from Fitted Point Process Model on a Linear Network
}
\description{
  Given a point process model fitted to a point pattern
  on a linear network, produce a Q-Q plot based on residuals from the model.
}
\usage{
  qqplot.lppm(fit, nsim=100, expr=NULL, \dots, type="raw",
             style="mean", fast=TRUE, verbose=TRUE, plot.it=TRUE,
             probs=NULL, 
             saveall=FALSE,
             monochrome=FALSE,
             limcol=if(monochrome) "black" else "red",
             maxerr=max(100, ceiling(nsim/10)),
             envir.expr)
}
\arguments{
  \item{fit}{
    The fitted point process model on a network, which is to be assessed
    using the Q-Q plot. An object of class \code{"lppm"}.
    Smoothed residuals obtained from this fitted model will provide the
    \dQuote{data} quantiles for the Q-Q plot.
  }
  \item{nsim}{
    The number of simulations from the \dQuote{reference} point process model.
  }
  \item{expr}{
    Determines the simulation mechanism
    which provides the \dQuote{theoretical} quantiles for the
    Q-Q plot. See Details.
  }
  \item{\dots}{
    Arguments passed to \code{\link{diagnose.lppm}} influencing the
    computation of residuals.
  }
  \item{type}{
    String indicating the type of residuals or weights to be used.
    Current options are \code{"eem"}
    for the Stoyan-Grabarnik exponential energy weights,
    \code{"raw"} for the raw residuals,
    \code{"inverse"} for the inverse-lambda residuals,
    and \code{"pearson"} for the Pearson residuals.
    A partial match is adequate.
  }
  \item{style}{
    Character string controlling the type of Q-Q plot.
    Options are \code{"classical"} and \code{"mean"}.
    See Details.
  }
  \item{fast}{
    Logical flag controlling the speed and accuracy of computation.
    Use \code{fast=TRUE} for interactive use and \code{fast=FALSE}
    for publication standard plots. See Details.
  }
  \item{verbose}{
    Logical value controlling whether the algorithm prints progress
    reports during long computations.
  }
  \item{plot.it}{
    Logical value controlling whether the function produces a plot
    or simply returns a value (silently).
  }
  \item{probs}{
    Probabilities for which the quantiles should be calculated and
    plotted. A numeric vector of values between 0 and 1.
  }
  \item{saveall}{
    Logical flag indicating whether to save all the intermediate
    calculations. 
  }
  \item{monochrome}{
    Logical flag indicating whether the plot should be
    in black and white (\code{monochrome=TRUE}), or in colour
    (\code{monochrome=FALSE}).
  }
  \item{limcol}{
    String. The colour to be used when plotting the 95\% limit
    curves.
  }
  \item{maxerr}{
    Maximum number of failures tolerated while generating
    simulated realisations. See Details.
  }
  \item{envir.expr}{
    Optional. An environment in which the expression
    \code{expr} should be evaluated.
  }
}
\value{
  An object of class \code{"qqlppm"} containing the information
  needed to reproduce the Q-Q plot.
  Entries \code{x} and \code{y} are numeric vectors containing 
  quantiles of the simulations and of the data, respectively.
}
\details{
  This function generates a Q-Q plot of the residuals from a
  fitted point process model on a linear network.
  It is an addendum to the suite of
  diagnostic plots produced by the function \code{\link{diagnose.lppm}},
  kept separate because it is computationally intensive. The
  quantiles of the theoretical distribution are estimated by simulation.

  In classical statistics, a Q-Q plot of residuals is a useful
  diagnostic for checking the distributional assumptions. Analogously,
  in spatial statistics, a Q-Q plot of the (smoothed) residuals from a
  fitted point process model is a useful way
  to check the interpoint interaction part of the model
  (Baddeley et al, 2005). The systematic part of the model
  (spatial trend, covariate effects, etc) is assessed using
  other plots made by \code{\link{diagnose.lppm}}.

  The argument \code{fit} represents the fitted point process
  model. It must be an object of class \code{"lppm"} (typically produced
  by the maximum pseudolikelihood fitting algorithm \code{\link{lppm}}).
  Residuals will be computed for this fitted model using
  \code{\link{residuals.lppm}},
  and the residuals will be kernel-smoothed to produce a
  \dQuote{residual field}. The values of this residual field will provide the
  \dQuote{data} quantiles for the Q-Q plot.

  The argument \code{expr} is not usually specified.
  It provides a way to modify the \dQuote{theoretical} or \dQuote{reference}
  quantiles for the Q-Q plot.

  In normal usage we set \code{expr=NULL}. The default
  is to generate \code{nsim} simulated realisations
  of the fitted model \code{fit} using
  \code{\link[spatstat.linnet]{simulate.lppm}}, 
  re-fit this model to each of the simulated patterns,
  evaluate the residuals from these fitted models,
  and use the kernel-smoothed residual field
  from these fitted models as a sample from the reference distribution
  for the Q-Q plot.
  
  In advanced use, \code{expr} may be an \code{expression}.
  It will be re-evaluated \code{nsim} times, and should include
  random computations so that the results are not identical
  each time. The result of evaluating \code{expr}
  should be either a point pattern on a network (object of class
  \code{"lpp"}) or a fitted point process model on a network (object of class
  \code{"lppm"}). If the value is a point pattern, then the
  original fitted model \code{fit} will be re-fitted to this new point
  pattern using \code{\link{update.lppm}}, to yield another fitted
  model. Smoothed residuals obtained from these \code{nsim} fitted
  models will yield the \dQuote{theoretical} quantiles for the Q-Q plot.

  Alternatively \code{expr} can be a list of point patterns,
  or an \code{envelope} object that contains a list of point patterns
  (typically generated by calling \code{\link[spatstat.linnet]{envelope.lpp}}
  or \code{\link[spatstat.linnet]{envelope.lppm}} with \code{savepatterns=TRUE}).
  These point patterns will be used as the simulated patterns. 

  The argument \code{type} selects the type of residual or weight
  that will be computed. For options, see \code{\link{diagnose.lppm}}.

  The argument \code{style} determines the type of Q-Q plot. It is
  highly recommended to use the default, \code{style="mean"}.
  \describe{
    \item{\code{style="classical"}}{
      The quantiles of the residual field for the data (on the \eqn{y}
      axis) are plotted against the
      quantiles of the \bold{pooled} simulations (on the \eqn{x} axis).
      This plot is biased, and therefore difficult to interpret,
      because of strong autocorrelations in the residual field
      and the large differences in sample size.
    }
    \item{\code{style="mean"}}{
      Quantiles of the residual field for the original data are plotted
      against the sample means, over the \code{nsim} simulations,
      of the corresponding quantiles of the residual field
      for the simulated datasets.
      Dotted lines show the 2.5 and 97.5 percentiles, over the
      \code{nsim} simulations, of each of these quantiles.
    }
  }

  The argument \code{fast} is a simple way to control
  the accuracy and speed of computation.
  If \code{fast=FALSE}, the residual field is computed on
  a fine grid of pixels (by default 100 by 100 pixels, see below)
  and the Q-Q plot is based on the complete set of order statistics
  (up to 10,000 quantiles).
  If \code{fast=TRUE}, the residual field is computed on a coarse
  grid (at most 40 by 40 pixels) and the Q-Q plot is based on the
  \emph{percentiles} only. This is about 7 times faster.
  It is recommended to use \code{fast=TRUE} for interactive data
  analysis and \code{fast=FALSE} for definitive plots for
  publication.

  Since the computation is so time-consuming, \code{qqplot.lppm} returns
  a list containing all the data necessary to re-display the Q-Q plot.
  It is advisable to assign the result of \code{qqplot.lppm} to something
  (or use \code{.Last.value} if you forgot to.)
  The return value is an object of class \code{"qqlppm"}. There are methods for
  \code{plot} and \code{print}. See the Examples.

  The argument \code{saveall} is usually set to \code{FALSE}.
  If \code{saveall=TRUE}, then the intermediate results of calculation for each
  simulated realisation are saved and returned. The return value
  includes a 3-dimensional array \code{sim} containing the
  smoothed residual field images for each of the \code{nsim} 
  realisations. When \code{saveall=TRUE}, the return value is an object of very
  large size, and should not be saved on disk.

  Errors may occur during the simulation process, because
  random data are generated. For example:
  \itemize{
    \item one of the simulated patterns may be empty.
    \item one of the simulated patterns may 
    cause an error in the code that fits the point process model.
    \item the user-supplied argument \code{expr} may have a bug.
  }
  Empty point patterns do not cause a problem for the code,
  but they are reported.
  Other problems that would lead to a crash are trapped; 
  the offending simulated data are discarded, and the simulation is
  retried. The argument \code{maxerr} determines the maximum number of
  times that such errors will be tolerated (mainly as a
  safeguard against an infinite loop).
}

\section{Side Effects}{
  Produces a Q-Q plot if \code{plot.it=TRUE} (the default).
}
\section{Warning messages}{
  A warning message will be issued if any of the simulations
  trapped an error (a potential crash).
  
  A warning message will be issued if all, or many, of the
  simulated point patterns are empty.
  This usually indicates a problem with the simulation procedure.
}

\references{
  Baddeley, A., Turner, R., \Moller, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Stoyan, D. and Grabarnik, P. (1991)
  Second-order characteristics for stochastic structures connected with
  Gibbs point processes.
  \emph{Mathematische Nachrichten}, 151:95--100.
}
\seealso{
 \code{\link{diagnose.lppm}},
 \code{\link{lurking}},
 \code{\link{residuals.lppm}},
 \code{\link{lppm}}
}
\examples{
    X <- unmark(chicago)
    if(interactive()) X <- X[c(FALSE,TRUE)]

    fit <- lppm(X ~ x + y)

    ns <- if(interactive()) 100 else 4
    qqplot.lppm(fit, nsim=ns)
    qqplot.lppm(fit, nsim=ns, type="pearson")

    \donttest{
      # capture the plot coordinates
      mypreciousdata <- qqplot.lppm(fit, ns, type="pearson")
      ## or use the idiom .Last.value if you forgot to assign them
      qqplot.lppm(fit, ns, type="pearson")
      mypreciousdata <- .Last.value
      plot(mypreciousdata)
    }
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{models}
\keyword{hplot}
