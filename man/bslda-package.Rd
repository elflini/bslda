\name{bslda-package}
\alias{bslda-package}
\docType{package}

\title{
Sparse Linear Discriminant Analysis using Network-Guided Block Covariance Matrix
}

\description{
This package provides functions for fitting bslda, a sparse linear discriminant analysis using network-guided block covariance matrix, which allows researchers to utilize prior network information to guide structure of covariance matrix.
}

\details{
\tabular{ll}{
Package: \tab bslda\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2019-02-28\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

This package contains a main class, \code{bslda}, which represents bslda model fit.

The main methods for the bslda is \code{bslda}, which fits sparse linear discriminant analysis using network-guided covariance matrix.
}

\author{
Jin Hyun Nam

Maintainer: Jin Hyun Nam <elflini@gmail.com>
}

\references{
Nam JH, Kim D, and Chung D (2019) ''Sparse linear discrimintn analysis
using network-guided block covariance matrix''.
}

\keyword{ package }

\seealso{
\code{\link{bslda}}
}


\examples{
library(bslda)

# Loading data

data(sdata)

# Fit BSLDA
fit.bslda <- bslda(x,y,network=new.path,alpha=0.5)
fit.bslda

# Extract discriminant coefficients
coeff(fit.bslda)

# Prediction of new individuals
predict(fit.bslda,xtest)

}
