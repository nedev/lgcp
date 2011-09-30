\name{spatialparsEst}
\alias{spatialparsEst}
\title{spatialparsEst function...}
\usage{spatialparsEst(gk, sigma.range, phi.range)}
\description{spatialparsEst function}
\details{Having estimated either the pair correlation or K functions using respectively \link{ginhomAverage} or \link{KinhomAverage}, the spatial 
parameters sigma and phi can be estimated. This function provides a visual tool for this estimation procedure.}
\value{rpanel function to help choose sigma nad phi by eye}
\references{\enumerate{
\item Baddeley AJ, Moller J, Waagepetersen R (2000). Non-and semi-parametric estimation of interaction in inhomogeneous point patterns. Statistica Neerlandica, 54, 329-350.
\item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
\item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
}}
\seealso{\link{ginhomAverage}, \link{KinhomAverage}, \link{thetaEst}, \link{lambdaEst}, \link{muEst}}
\arguments{\item{gk}{an R object; output from the function KinhomAverage or ginhomAverage}
\item{sigma.range}{range of sigma values to consider}
\item{phi.range}{range of phi values to consider}}
