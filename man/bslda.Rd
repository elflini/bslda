\name{bslda}
\alias{bslda}

\title{Sparse Linear Discriminant Analysis using Network-Guided Block Covariance Matrix}

\description{
Fit sparse linear discriminant analysis using
network-guided block covariance matrix
}

\usage{
bslda(x,y,network=NA,alpha=0.5)
}

\arguments{
\item{x}{Matrix with n samples by p variables}

\item{y}{class vector}

\item{network}{List. Variable list for the network information (default: NA)}

\item{alpha}{Weight of L2 norm ranges in (0,1) (default: 0.5)}
}

\details{
\code{bslda} fits the sparse linear discriminant analysis model using network-guided covariance matrix. It requires list for network information and if network information is missing then bslda fits with whole covariance matrix.
}

\value{
Construct \code{bslda} class object.
}

\references{
Nam JH, Kim D, and Chung D (2019) ''Sparse linear discrimintn analysis
using network-guided block covariance matrix''.
}


\author{ Jin Hyun Nam }

\examples{
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
