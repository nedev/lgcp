\name{fftgrid}
\alias{fftgrid}
\title{fftgrid function...}
\usage{fftgrid(xyt, M, N, spatial, sigma, phi, model, covpars)}
\description{fftgrid function}
\details{\bold{Advanced use only.} Computes various quantities for use in \code{lgcpPredict},
\code{lgcpSim} .}
\value{fft objects for use in MALA}
\arguments{\item{xyt}{object of class stppp}
\item{M}{number of centroids in x-direction}
\item{N}{number of centroids in y-direction}
\item{spatial}{an object of class spatialAtRisk}
\item{sigma}{scaling paramter for spatial covariance function, see Brix and Diggle (2001)}
\item{phi}{scaling paramter for spatial covariance function, see Brix and Diggle (2001)}
\item{model}{correlation type see ?CovarianceFct}
\item{covpars}{vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct}}
\examples{xyt <- stppp(ppp(),t=c(),tlim=c(0,1))
sar <- spatialAtRisk(function(x,y){return(1)},warn=FALSE)
grid <- fftgrid(xyt,64,64,sar,1,0.01,"exponential",covpars=c())}
