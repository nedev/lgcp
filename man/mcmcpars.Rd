\name{mcmcpars}
\alias{mcmcpars}
\title{mcmcpars function...}
\usage{mcmcpars(mala.length, burnin, retain, inits, MCMCdiag=0,
    adaptivescheme)}
\description{mcmcpars function}
\details{A function for setting MCMC options in a run of \code{lgcpPredict} for example.}
\value{mcmc parameters}
\seealso{\link{lgcpPredict}}
\arguments{\item{mala.length}{default = 100,}
\item{burnin}{default = floor(mala.length/2),}
\item{retain}{thinning parameter eg operated on chain every 'retain' iteration (eg store output or compute some posterior functional)}
\item{inits}{optional initial values for MCMC}
\item{MCMCdiag}{non-negative integer, if greater than zero saves information from the MCMC chain}
\item{adaptivescheme}{the type of adaptive mcmc to use, see ?constanth (constant h) or ?andrieuthomsh (adaptive MCMC of Andrieu and Thoms (2008))}}
