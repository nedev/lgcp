\name{showGrid.stppp}
\alias{showGrid.stppp}
\title{showGrid...}
\usage{\method{showGrid}{stppp} (x, ...)}
\description{showGrid.stppp function}
\details{If an stppp object has been created via simulation, ie using the function \code{lgcpSim}, then
this function will display the grid centroids that were used in the simulation}
\value{plots grid centroids on the current graphics device. FOR SIMULATED DATA ONLY.}
\seealso{\link{lgcpSim}, \link{showGrid.default}, \link{showGrid.lgcpPredict}}
\arguments{\item{x}{an object of class stppp. Note this function oly applies to SIMULATED data.}
\item{...}{additional arguments  passed to points}}
\examples{\dontrun{xyt <- lgcpSim()}
\dontrun{plot(xyt)}
\dontrun{showGrid(xyt,pch="+",cex=0.5)}}
