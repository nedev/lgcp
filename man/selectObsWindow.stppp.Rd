\name{selectObsWindow.stppp}
\alias{selectObsWindow.stppp}
\title{selectObsWindow...}
\usage{\method{selectObsWindow}{stppp} (xyt, gridsize, ...)}
\description{selectObsWindow.stppp function}
\details{This function computes an appropriate observation window on which to perform prediction. Since the FFT grid
must have dimension 2^M by 2^N for some M and N, the window \code{xyt$window}, is extended to allow this to be fit in for a given cell width.}
\value{a resized stppp object together with grid sizes M and N ready for FFT}
\arguments{\item{xyt}{an object of class stppp}
\item{gridsize}{size of the grid spacing in chosen units (equivalent to the cell width argument in \link{lgcpPredict})}
\item{...}{additional arguments}}
