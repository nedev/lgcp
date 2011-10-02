\name{lambdaEst}
\alias{lambdaEst}
\title{lambdaEst function...}
\usage{lambdaEst(xyt, weights=c(), edge=TRUE, ...)}
\description{lambdaEst function}
\details{A tool for the visual estimation of lambda(s) via a 2 dimensional smoothing of the case locations. For parameter estimation, the alternative is
to estimate lambda(s) by some other means, convert it into a spatialAtRisk object and then into a pixel image object using the build in coercion 
methods, this \code{im} object can then be fed to \link{ginhomAverage}, \link{KinhomAverage} or \link{thetaEst} for instance.}
\value{This is an rpanel function for visual choice of lambda(s), the output is a variable, varname, with the density *per unit time* 
the variable varname can be fed to the function ginhomAverage or KinhomAverage as the argument density (see for example ?ginhomAverage), or into the 
function thetaEst as the argument spatial.intensity.}
\references{\enumerate{
\item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
\item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
}}
\seealso{\link{spatialAtRisk}, \link{ginhomAverage}, \link{KinhomAverage}, \link{spatialparsEst}, \link{thetaEst}, \link{muEst}}
\arguments{\item{xyt}{object of class stppp}
\item{weights}{Optional vector of weights to be attached to the points.  May include negative values. See ?density.ppp.}
\item{edge}{Logical flag: if TRUE, apply edge correction. See ?density.ppp.}
\item{...}{arguments to be passed to plot}}
