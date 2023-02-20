
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fake: Flexible Data Simulation Using The Multivariate Normal Distribution <img src="man/figures/logo.png" align="right" width="174" height="200"/>

<!-- badges: start -->

[![CRAN
version](https://www.r-pkg.org/badges/version/fake)](https://cran.r-project.org/package=fake)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-month/fake?color=blue)](https://r-pkg.org/pkg/fake)
![GitHub last
commit](https://img.shields.io/github/last-commit/barbarabodinier/fake?logo=GitHub&style=flat-square)
<!-- badges: end -->

## Description

> This R package can be used to generate artificial data conditionally
> on pre-specified (simulated or user-defined) relationships between the
> variables and/or observations. Each observation is drawn from a
> multivariate Normal distribution where the mean vector and covariance
> matrix reflect the desired relationships. Outputs can be used to
> evaluate the performances of variable selection, graphical modelling,
> or clustering approaches by comparing the true and estimated
> structures.

## Installation

The released version of the package can be installed from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("fake")
```

The development version can be installed from
[GitHub](https://github.com/):

``` r
remotes::install_github("barbarabodinier/fake")
```

## Main functions

### Linear model

``` r
library(fake)

set.seed(1)
simul <- SimulateRegression(n = 100, pk = 20)
head(simul$xdata)
head(simul$ydata)
```

### Logistic model

``` r
set.seed(1)
simul <- SimulateRegression(n = 100, pk = 20, family = "binomial")
head(simul$ydata)
```

### Structural causal model

``` r
set.seed(1)
simul <- SimulateStructural(n = 100, pk = c(3, 2, 3))
head(simul$data)
```

### Gaussian graphical model

``` r
set.seed(1)
simul <- SimulateGraphical(n = 100, pk = 20)
head(simul$data)
```

### Gaussian mixture model

``` r
set.seed(1)
simul <- SimulateClustering(n = c(10, 10, 10), pk = 20)
head(simul$data)
```

## Extraction and visualisation of the results

The true model structure is returned in the output of any of the main
functions in:

``` r
simul$theta
```

The functions `print()`, `summary()` and `plot()` can be used on the
outputs from the main functions.

## Reference

- Barbara Bodinier, Sarah Filippi, Therese Haugdahl Nost, Julien Chiquet
  and Marc Chadeau-Hyam. Automated calibration for stability selection
  in penalised regression and graphical models: a multi-OMICs network
  application exploring the molecular response to tobacco
  smoking. (2021) arXiv.
  [link](https://doi.org/10.48550/arXiv.2106.02521)

## Other resources

- R scripts to reproduce the simulation study (Bodinier et al. 2021)
  conducted using the functions in **fake**
  [link](https://github.com/barbarabodinier/stability_selection)

- R package **sharp** for stability selection and consensus clustering
  [link](https://github.com/barbarabodinier/sharp)
