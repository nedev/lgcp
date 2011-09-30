\name{lgcppars}
\alias{lgcppars}
\title{lgcppars function...}
\usage{lgcppars(sigma, phi, theta)}
\description{lgcppars function}
\details{A function for setting the parameters sigma, phi and theta for \code{lgcpPredict}. Note that the returned
set of parameters also features mu=-0.5*sigma^2, gives mean(exp(Y)) = 1.}
\seealso{\link{lgcpPredict}}
\arguments{\item{sigma}{sigma parameter}
\item{phi}{phi parameter}
\item{theta}{this is beta parameter in Brix and Diggle (2001)}}
