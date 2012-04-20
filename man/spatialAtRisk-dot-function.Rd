\name{spatialAtRisk.function}
\alias{spatialAtRisk.function}
\title{spatialAtRisk...}
\usage{\method{spatialAtRisk}{function} (X, warn=TRUE, ...)}
\description{spatialAtRisk.function function}
\details{Creates a spatialAtRisk object from a function mapping R^2 onto the non negative reals. Note that for spatialAtRisk
objects defined in this manner, the user is responsible for ensurng that the integral of the function is 1 over the observation
window of interest.}
\value{object of class spatialAtRisk
NOTE The function provided is assumed to integrate to 1 over the observation window, the user is responsible for ensuring this is the case.
\enumerate{
\item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
\item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
}}
\seealso{\link{lgcpPredict}, link{lgcpSim}, \link{spatialAtRisk.default}, \link{spatialAtRisk.fromXYZ}, \link{spatialAtRisk.im}, \link{spatialAtRisk.SpatialGridDataFrame}, \link{spatialAtRisk.SpatialPolygonsDataFrame}, \link{spatialAtRisk.bivden}}
\arguments{\item{X}{a function with accepts arguments x and y that returns the at risk population at coordinate (x,y), which should be a numeric of length 1}
\item{warn}{whether to issue a warning or not}
\item{...}{additional arguments}}
