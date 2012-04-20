\name{setoutput}
\alias{setoutput}
\title{setoutput function...}
\usage{setoutput(gridfunction, gridmeans)}
\description{setoutput function}
\details{Sets output functionality for \link{lgcpPredict} via the main functions \link{dump2dir} and \link{MonteCarloAverage}. Note that it is possible for
the user to create their own \code{gridfunction} and \code{gridmeans} schemes.}
\value{output parameters}
\seealso{\link{lgcpPredict}, \link{dump2dir}, \link{MonteCarloAverage}}
\arguments{\item{gridfunction}{what to do with the latent field, but default this set to nothing, but could save output to a directory, see ?dump2dir}
\item{gridmeans}{list of Monte Carlo averages to compute, see ?MonteCarloAverage}}
