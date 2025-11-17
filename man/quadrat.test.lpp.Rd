\name{quadrat.test.lpp} 
\alias{quadrat.test.lpp} 
\alias{quadrat.test.lppm} 
\alias{quadrat.test.linearquadratcount} 
\title{Dispersion Test for Point Pattern on a Network Based on Quadrat Counts}
\description{
  Performs a test of Complete Spatial Randomness
  for a given point pattern on a linear network, based on quadrat counts.
  Alternatively performs a goodness-of-fit test of a fitted
  inhomogeneous Poisson model on a network.
  By default performs chi-squared tests; can also perform
  power-divergence tests and Monte Carlo tests.
}
\usage{
\method{quadrat.test}{lpp}(X, \dots,
                           tess=NULL,
                           nx=5, ny=nx,
                           xbreaks=NULL, ybreaks=NULL, 
                           alternative=c("two.sided", "regular", "clustered"),
                           method=c("Chisq", "MonteCarlo"),
                           conditional=TRUE, CR=1,
                           lambda=NULL, df.est=NULL,
                           nsim=1999)
\method{quadrat.test}{lppm}(X,\dots,
                           tess=NULL,
                           nx=5, ny=nx,
                           xbreaks=NULL, ybreaks=NULL, 
                           alternative=c("two.sided", "regular", "clustered"),
                           method=c("Chisq", "MonteCarlo"),
                           conditional=TRUE, CR=1, df.est=NULL,
                           nsim=1999)
\method{quadrat.test}{linearquadratcount}(X, \dots,
                           alternative=c("two.sided", "regular", "clustered"),
                           method=c("Chisq", "MonteCarlo"),
                           conditional=TRUE, CR=1,
                           lambda=NULL, df.est=NULL,
                           nsim=1999)
}
\arguments{
  \item{X}{
    A point pattern on a network (object of class \code{"lpp"})
    to be subjected to the goodness-of-fit test.
    Alternatively a fitted point process model on a network (object of class
    \code{"lppm"}) to be tested.
    Alternatively \code{X} can be the result of applying
    \code{\link[spatstat.geom]{quadratcount}} to a point pattern on a network.
  }
  \item{tess}{
    Tessellation (object of class \code{"tess"} or \code{"lintess"})
    determining the quadrats. Incompatible with \code{nx,ny,xbreaks,ybreaks}.
  }
  \item{nx,ny}{
    Numbers of quadrats in the \eqn{x} and \eqn{y} directions.
    Incompatible with \code{xbreaks} and \code{ybreaks}.
  }
  \item{xbreaks}{
    Optional. Numeric vector giving the \eqn{x} coordinates of the
    boundaries of the quadrats. Incompatible with \code{nx}.
  }
  \item{ybreaks}{
    Optional. Numeric vector giving the \eqn{y} coordinates of the
    boundaries of the quadrats. Incompatible with \code{ny}.
  }
  \item{alternative}{
    Character string (partially matched) specifying the alternative
    hypothesis.
  }
  \item{method}{
    Character string (partially matched) specifying the test to use:
    either \code{method="Chisq"} for the chi-squared test (the default),
    or \code{method="MonteCarlo"} for a Monte Carlo test.
  }
  \item{conditional}{
    Logical.  Should the Monte Carlo test be conducted
    conditionally upon the observed number of points of the pattern?
    Ignored if \code{method="Chisq"}.
  }
  \item{CR}{
    Optional. Numerical value. The exponent
    for the Cressie-Read test statistic. See Details.
  }
  \item{lambda}{
    Optional. Pixel image (object of class \code{"im"} or \code{"linim"})
    or function (class \code{"funxy"}) giving the predicted
    intensity of the point process.
  }
  \item{df.est}{
    Optional. Advanced use only.
    The number of fitted parameters, or the
    number of degrees of freedom lost by estimation of
    parameters. 
  }
  \item{\dots}{Ignored.}
  \item{nsim}{
    The number of simulated samples to generate when
    \code{method="MonteCarlo"}. 
  }
}
\details{
  These functions perform \eqn{\chi^2}{chi^2} tests or Monte Carlo tests
  of goodness-of-fit for a point process model on a linear network,
  based on quadrat counts.

  The function \code{\link[spatstat.explore]{quadrat.test}} is generic, with methods for
  many classes. This page documents the methods for data on a linear
  network.

  \itemize{
    \item
    if \code{X} is a point pattern on a network (object of class
    \code{"lpp"}),
    we test the null hypothesis
    that the data pattern is a realisation of Complete Spatial
    Randomness (the uniform Poisson point process) on the network.
    Marks in the point
    pattern are ignored. (If \code{lambda} is given then the null
    hypothesis is the Poisson process with intensity \code{lambda}.)
    \item
    If \code{X} is a fitted point process model on a network (object of
    class \code{"lppm"}), then it should be
    a Poisson point process model. The 
    data to which this model was fitted are extracted from the model
    object, and are treated as the data point pattern for the test.
    We test the null hypothesis 
    that the data pattern is a realisation of the (inhomogeneous) Poisson point
    process specified by \code{X}.
  }

  First the network is divided into pieces to form a
  tessellation (object of class \code{"lintess"}) as follows:
  \itemize{
    \item
    By default, if none of the arguments
    \code{nx, ny, xbreaks, ybreaks, tess} is given,
    every segment of the network is taken as a separate piece.
    The number of points in each segment of the network is counted.
    \item
    If \code{nx, ny} are given, the window containing
    the point pattern \code{X} is divided into
    an \code{nx * ny} grid of rectangular tiles or `quadrats'.
    These tiles are then intersected with the network on which \code{X}
    is defined. The number of points falling in each rectangle
    is counted.
    \item
    If \code{xbreaks} is given, 
    the window containing the point pattern \code{X} will be divided
    into rectangles, with \code{xbreaks} and \code{ybreaks}
    giving the \eqn{x} and \eqn{y} coordinates of the rectangle
    boundaries, respectively.
    The lengths of \code{xbreaks} and \code{ybreaks} may be different.
    \item 
    The argument \code{tess} can be a tessellation on the network
    (object of class \code{"lintess"}) whose tiles will serve as the
    quadrats.
    \item
    Alternatively \code{tess} can be a two-dimensional tessellation
    (object of class \code{"tess"}) which will be intersected with the
    network to determine the tessellation of the network.
  }

  Next the number of data points in each tile of the tessellation is 
  counted.

  The expected number of points in each quadrat is also calculated,
  as determined by CSR (in the first case) or by the fitted model
  (in the second case).

  Then the Pearson \eqn{X^2} statistic 
  \deqn{
    X^2 = sum((observed - expected)^2/expected)
  }
  is computed.
  
  If \code{method="Chisq"} then a \eqn{\chi^2}{chi^2} test of
  goodness-of-fit is performed by comparing the test statistic
  to the \eqn{\chi^2}{chi^2} distribution
  with \eqn{m-k} degrees of freedom, where \code{m} is the number of
  quadrats and \eqn{k} is the number of fitted parameters
  (equal to 1 for \code{quadrat.test.ppp}). The default is to
  compute the \emph{two-sided} \eqn{p}-value, so that the test will
  be declared significant if \eqn{X^2} is either very large or very
  small. One-sided \eqn{p}-values can be obtained by specifying the
  \code{alternative}. An important requirement of the
  \eqn{\chi^2}{chi^2} test is that the expected counts in each quadrat
  be greater than 5.

  If \code{method="MonteCarlo"} then a Monte Carlo test is performed,
  obviating the need for all expected counts to be at least 5.  In the
  Monte Carlo test, \code{nsim} random point patterns are generated
  from the null hypothesis (either CSR or the fitted point process
  model). The Pearson \eqn{X^2} statistic is computed as above.
  The \eqn{p}-value is determined by comparing the \eqn{X^2}
  statistic for the observed point pattern, with the values obtained
  from the simulations. Again the default is to
  compute the \emph{two-sided} \eqn{p}-value.

  If \code{conditional} is \code{TRUE} then the simulated samples are
  generated from the multinomial distribution with the number of \dQuote{trials}
  equal to the number of observed points and the vector of probabilities
  equal to the expected counts divided by the sum of the expected counts.
  Otherwise the simulated samples are independent Poisson counts, with
  means equal to the expected counts.

  If the argument \code{CR} is given, then instead of the
  Pearson \eqn{X^2} statistic, the Cressie-Read (1984) power divergence
  test statistic
  \deqn{
    2nI = \frac{2}{CR(CR+1)}
           \sum_i \left[ \left( \frac{X_i}{E_i} \right)^CR - 1 \right]
  }{
    2nI = (2/(CR * (CR+1))) * sum((X[i]/E[i])^CR - 1)
  }
  is computed, where \eqn{X_i}{X[i]} is the \eqn{i}th observed count
  and \eqn{E_i}{E[i]} is the corresponding expected count.
  The value \code{CR=1} gives the Pearson \eqn{X^2} statistic;
  \code{CR=0} gives the likelihood ratio test statistic \eqn{G^2};
  \code{CR=-1/2} gives the Freeman-Tukey statistic \eqn{T^2};
  \code{CR=-1} gives the modified likelihood ratio test statistic \eqn{GM^2};
  and \code{CR=-2} gives Neyman's modified statistic \eqn{NM^2}.
  In all cases the asymptotic distribution of this test statistic is
  the same \eqn{\chi^2}{chi^2} distribution as above.

  The return value is an object of class \code{"htest"}.
  Printing the object gives comprehensible output
  about the outcome of the test.

  The return value also belongs to
  the special class \code{"quadrat.test"}. Plotting the object
  will display the quadrats, annotated by their observed and expected
  counts and the Pearson residuals. See the examples.
}
\seealso{
  \code{\link[spatstat.explore]{quadrat.test}},
  \code{\link[spatstat.linnet]{quadratcount.lpp}},
  \code{\link[spatstat.linnet]{lintess}}.
}
\value{
  An object of class \code{"htest"}. See \code{\link[stats]{chisq.test}}
  for explanation. 

  The return value is also an object of the special class
  \code{"quadrattest"}, and there is a plot method for this class.
  See the examples.
}
\references{
  Cressie, N. and Read, T.R.C. (1984)
  Multinomial goodness-of-fit tests.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{46}, 440--464. 
}
\examples{
  X <- runiflpp(100, simplenet)

  quadrat.test(X)

  quadrat.test(X, nx=2)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{htest}
\concept{Test of randomness}
\concept{Test of clustering}
