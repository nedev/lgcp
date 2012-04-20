\name{fftinterpolate.fromFunction}
\alias{fftinterpolate.fromFunction}
\title{fftinterpolate...}
\usage{\method{fftinterpolate}{fromFunction} (spatial, mcens, ncens, ext, ...)}
\description{fftinterpolate.fromFunction function}
\details{This method performs interpolation within the function \code{fftgrid} for \code{fromFunction} objects.}
\value{matrix of interpolated values}
\seealso{\link{fftgrid}, \link{spatialAtRisk.function}}
\arguments{\item{spatial}{objects of class spatialAtRisk}
\item{mcens}{x-coordinates of interpolation grid in extended space}
\item{ncens}{y-coordinates of interpolation grid in extended space}
\item{ext}{integer multiple by which grid should be extended, default is 2. Generally this will not need to be altered, but if the spatial correlation decays slowly, increasing 'ext' may be necessary.}
\item{...}{additional arguments}}
