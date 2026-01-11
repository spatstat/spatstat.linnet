# spatstat.linnet

## Spatial analysis on a linear network, for the spatstat family

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spatstat.linnet)](http://CRAN.R-project.org/package=spatstat.linnet) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.linnet)](https://github.com/spatstat/spatstat.linnet)

You are viewing the GitHub repository which holds
the latest **development version** of `spatstat.linnet`.
For the latest public release on CRAN, click the green badge above.

 - [Overview of `spatstat.linnet`](#overview)
 - [Where to find data](#data)
 - [Detailed contents of package](#detailed)
 - [Installing the package](#installing)
 - [Bug reports](#bugreports)
 - [Questions](#questions)
 - [Proposing changes to code](#proposing)
 - [Future development](#future)

____

### <a name="overview"></a> Overview of `spatstat.linnet`

The original `spatstat` package has been split into several sub-packages
(See [spatstat/spatstat](https://github.com/spatstat/spatstat)).

This package `spatstat.linnet` is one of the sub-packages. 
It contains the subset of the functionality of `spatstat`
that deals with **data on linear networks**. It supports

 - network geometry
 - point patterns on a network
 - spatial covariates on a network
 - simulation
 - exploratory data analysis
 - parametric modelling and formal inference
 - informal model diagnostics

There is also an extension package
[spatstat.Knet](https://github.com/spatstat/spatstat.Knet)
which contains additional algorithms for linear networks.

___

### <a name="data"></a> Where to find data

Examples of datasets on linear networks are
the point patterns `chicago`, `dendrite` and `spiders` provided in the
[spatstat.data](https://github.com/spatstat/spatstat.data)
package (available when `spatstat.linnet` is loaded)
and the point pattern `wacrashes` provided in the extension package
[spatstat.Knet](https://github.com/spatstat/spatstat.Knet)
(which must be loaded separately).

___

### <a name="detailed"></a> Detailed contents of `spatstat.linnet`

`spatstat.linnet` supports

#### Network geometry

- creation of linear networks from coordinate data
- extraction of networks from tessellations
- modification of networks 
- interactive editing of networks
- geometrical operations and measurement on networks
- construction of the disc in the shortest-path metric
- trees, tree branch labels, tree pruning

#### Point patterns on a network

- creation of point patterns on a network from coordinate data
- extraction of sub-patterns
- shortest-path distance measurement

#### Covariates on a network

- create pixel images and functions on a network
- arithmetic operators for pixel images on a network
- plot pixel images on a network (colour/thickness/perspective)
- tessellation on a network

#### Simulation

- completely random (uniform Poisson) point patterns on a network
- nonuniform random (Poisson) point patterns on a network
- Switzer-type point process
- log-Gaussian Cox process

#### Exploratory analysis of point patterns on a network

- kernel density estimation on a network
- bandwidth selection
- kernel smoothing on a network
- estimation of intensity as a function of a covariate
- ROC curves
- Berman-Waller-Lawson test
- CDF test
- variable selection by Sufficient Dimension Reduction
- K function on a network (shortest path or Euclidean distance)
- pair correlation function on a network (shortest path or Euclidean distance)
- inhomogeneous K function and pair correlation function
- inhomogeneous F, G and J functions
- simulation envelopes of summary functions

#### Parametric modelling and inference on a network

- fit point process model on a network
- fitted/predicted intensity
- analysis of deviance for point process model
- simulate fitted model

#### Informal model diagnostics

- lurking variable plot
- residuals
- leverage and influence
- four-panel diagnostic plot
- residual Q-Q plot

___

### <a name="installing"></a> Installing the package

This repository contains the _development version_ of
`spatstat.linnet`. The easiest way to install the development version
is to start R and type

```R
repo <- c('https://spatstat.r-universe.dev', 'https://cloud.r-project.org')
install.packages("spatstat.linnet", dependencies=TRUE, repos=repo)
```

To install the latest _public release_ of `spatstat.linnet`,
type

```R
install.packages("spatstat.linnet")
```


___

## <a name="bugreports"></a> Bug reports 

Users are encouraged to report bugs.
If you find a bug in a `spatstat` function,
please identify the sub-package containing that function.
Visit the GitHub repository for the sub-package, 
click the `Issues` tab at the top of the page, 
and press *new issue* to start a new bug report, documentation correction
or feature request.

**Please do not post questions** on the Issues pages,
because they are too clunky for correspondence.

## <a name="questions"></a> Questions about spatstat

For questions about the `spatstat` package family, first check 
the question-and-answer website
[stackoverflow](http://stackoverflow.com/questions/tagged/spatstat)
to see whether your question has already been asked and answered.
If not, you can either post your question at stackoverflow, or
email the authors.

## <a name="proposing"></a> Proposing changes to the code

Feel free to fork `spatstat.linnet`, make changes to the code,
and ask us to include them in the package by making a github *pull request*. 

## <a name="future"></a> Future development

The `spatstat` package family is the result of 30 years of software development
and contains over 200,000 lines of code.
It is still under development,
motivated by the needs of researchers in many fields,
and driven by innovations in statistical science.
We welcome contributions of code, and suggestions
for improvements.

