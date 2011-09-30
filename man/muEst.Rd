\name{muEst}
\alias{muEst}
\title{muEst function...}
\usage{muEst(xyt, ...)}
\description{muEst function}
\details{Computes a non-parametric estimate of mu(t). For the purposes of performing prediction, the alternatives are: (1) use a parameteric model as in Diggle P, Rowlingson B, Su T (2005),
or (2) a \link{constantInTime} model.}
\value{object of class temporalAtRisk giving the smoothed mut using the lowess function}
\references{\enumerate{
\item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
\item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
}}
\seealso{\link{temporalAtRisk}, \link{constantInTime}, \link{ginhomAverage}, \link{KinhomAverage}, \link{spatialparsEst}, \link{thetaEst}, \link{lambdaEst}}
\arguments{\item{xyt}{an stppp object}
\item{...}{additional arguments to be passed to lowess}}
