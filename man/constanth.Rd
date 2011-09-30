\name{constanth}
\alias{constanth}
\title{constanth function...}
\usage{constanth(h)}
\description{constanth function}
\details{This function is used to set up a constant acceptance scheme in the argument 
\code{mcmc.control} of the function \code{lgcpPredict}. The scheme is set via
the \code{mcmcpars} function.}
\value{object of class constanth}
\seealso{\link{mcmcpars}, \link{lgcpPredict}}
\arguments{\item{h}{an object}}
\examples{constanth(0.01)}
