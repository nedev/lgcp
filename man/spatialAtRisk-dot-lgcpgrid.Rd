\name{spatialAtRisk.lgcpgrid}
\alias{spatialAtRisk.lgcpgrid}
\title{spatialAtRisk...}
\usage{\method{spatialAtRisk}{lgcpgrid} (X, idx=length(X$grid), ...)}
\description{spatialAtRisk.lgcpgrid function}
\details{Creates a spatialAtRisk object from an lgcpgrid object}
\value{object of class spatialAtRisk
\enumerate{
\item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
\item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
}}
\seealso{\link{lgcpPredict}, link{lgcpSim}, \link{spatialAtRisk.default}, \link{spatialAtRisk.fromXYZ}, \link{spatialAtRisk.im}, \link{spatialAtRisk.function}, \link{spatialAtRisk.SpatialGridDataFrame}, \link{spatialAtRisk.SpatialPolygonsDataFrame}}
\arguments{\item{X}{an lgcpgrid object}
\item{idx}{in the case that X$grid is a list of length > 1, this argument specifies which element of the list to convert. By default, it is the last.}
\item{...}{additional arguments}}
