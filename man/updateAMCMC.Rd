\name{updateAMCMC}
\alias{updateAMCMC}
\title{updateAMCMC function...}
\usage{updateAMCMC(obj, ...)}
\description{updateAMCMC function}
\details{A generic to be used for the purpose of user-defined adaptive MCMC schemes,
updateAMCMC tells the MALA algorithm how to update the value of h. See
lgcp vignette, code{vignette("lgcp")}, for further details on writing adaptive MCMC schemes.}
\value{method updateAMCMC}
\seealso{\link{updateAMCMC.constanth}, \link{updateAMCMC.andrieuthomsh}}
\arguments{\item{obj}{an object}
\item{...}{additional arguments}}
