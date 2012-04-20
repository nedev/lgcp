\name{selectObsWindow.default}
\alias{selectObsWindow.default}
\title{selectObsWindow...}
\usage{\method{selectObsWindow}{default} (xyt, cellwidth, ...)}
\description{selectObsWindow.default function}
\details{Default method, note at present, there is only an implementation for stppp objects.

!!NOTE!! that this function also returns the grid ($xvals and $yvals) on which the FFT (and hence MALA) will be performed. It is useful to
define spatialAtRiskobjects on this grid to prevent loss of information from the bilinear interpolation that takes place as part of the fitting 
algorithm.}
\value{this is the same as selectObsWindow.stppp}
\seealso{\link{spatialAtRisk} \link{selectObsWindow.stppp}}
\arguments{\item{xyt}{an object}
\item{cellwidth}{size of the grid spacing in chosen units (equivalent to the cell width argument in \link{lgcpPredict})}
\item{...}{additional arguments}}
