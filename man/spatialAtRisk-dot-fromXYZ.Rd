\name{spatialAtRisk.fromXYZ}
\alias{spatialAtRisk.fromXYZ}
\title{spatialAtRisk...}
\usage{\method{spatialAtRisk}{fromXYZ} (X, Y, Zm, ...)}
\description{spatialAtRisk.fromXYZ function}
\details{Creates a spatialAtRisk object from a list of X, Y, Zm giving respectively the x and y coordinates of the grid and the 'z' values ie
so that Zm[i,j] is proportional to the at-risk population at X[i], Y[j].}
\value{object of class spatialAtRisk
\enumerate{
\item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
\item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
}}
\seealso{\link{lgcpPredict}, link{lgcpSim}, \link{spatialAtRisk.default}, \link{spatialAtRisk.im}, \link{spatialAtRisk.function}, \link{spatialAtRisk.SpatialGridDataFrame}, \link{spatialAtRisk.SpatialPolygonsDataFrame}, \link{spatialAtRisk.bivden}}
\arguments{\item{X}{vector of x-coordinates}
\item{Y}{vector of y-coordinates}
\item{Zm}{matrix such that Zm[i,j] = f(x[i],y[j]) for some function f}
\item{...}{additional arguments}}
