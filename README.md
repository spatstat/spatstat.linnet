# spatstat.linnet

## Spatial analysis on a linear network, for the spatstat family

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat.linnet)](http://CRAN.R-project.org/package=spatstat.linnet) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.linnet)](https://github.com/spatstat/spatstat.linnet)

The original `spatstat` package has been split into several sub-packages
(See [spatstat/spatstat](https://github.com/spatstat/spatstat)).

This package `spatstat.linnet` is one of the sub-packages. 
It contains the subset of the functionality of `spatstat`
that deals with data on linear networks.

There is also an extension package
[spatstat.Knet](https://github.com/spatstat/spatstat.Knet)
which contains additional algorithms for linear networks.

### Where to find data

Examples of datasets on linear networks are
the point patterns `chicago`, `dendrite` and `spiders` provided in the
[spatstat.data](https://github.com/spatstat/spatstat.data)
package (available when `spatstat.linnet` is loaded)
and the point pattern `wacrashes` provided in the extension package
[spatstat.Knet](https://github.com/spatstat/spatstat.Knet)
(which must be loaded separately).

### Overview of `spatstat.linnet`

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




