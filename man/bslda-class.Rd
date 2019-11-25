\name{bslda-class}
\alias{bslda-class}
\docType{class}

\title{Class "bslda"}

\description{
This class represents sparse linear discriminant analysis using
network-guided block covariance matrix fit.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bslda", ...)}.
}

\section{Slots}{
\describe{

\item{\code{data}}{Object of class \code{"list"},
    	representing the raw data for bslda model. }

\item{\code{init}}{Object of class \code{"list"},
    	representing the setting for bslda model fitting.}

\item{\code{result}}{Object of class \code{"list"},
    	representing the fitted bslda model.}
}}


\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "bslda")}: provide brief summary of the object. }

    \item{coeff}{\code{signature(object = "bslda")}:
		provide estimated discriminant coefficients.}

		    \item{predict}{\code{signature(object = "bslda")}:
		provide prediced class of new individuals.}

}}

\details{

}

\references{
Nam JH, Kim D, and Chung D (2019) ''Sparse linear discrimintn analysis
using network-guided block covariance matrix''.
}

\author{ Jin Hyun Nam }

\seealso{
\code{\link{bslda}}.
}

\examples{
showClass("bslda")

data(sdata)

fit.bslda <- bslda(x,y,network=new.path,alpha=0.5)
fit.bslda
coeff(fit.bslda)
predict(fit.bslda,xtest)

}

\keyword{classes}
