\name{computeGradtrunc}
\alias{computeGradtrunc}
\title{computeGradtrunc function...}
\usage{computeGradtrunc(nsims=100, scale=1, nis, theta, timediffs, mu,
    temporal, grid)}
\description{computeGradtrunc function}
\details{\bold{Advanced use only.} A function to compute a gradient truncation parameter for MALA via simulation. The function
requires an FFT 'grid' to be pre-computed, see \link{fftgrid}.}
\value{gradient truncation parameter}
\seealso{\link{fftgrid}}
\arguments{\item{nsims}{The number of simulations to use in computation of gradient truncation.}
\item{scale}{multiplicative scaling constant, returned value is scale*max(gradient over simulations). Default scale is 1.}
\item{nis}{cell counts on the extended grid}
\item{theta}{temporal correlation parameter}
\item{timediffs}{vector of time-differences between time intervals under consideration}
\item{mu}{mean parameter, see Brix and Diggle}
\item{temporal}{temporal fitted values}
\item{grid}{output from function fftgrid}}
