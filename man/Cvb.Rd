\name{Cvb}
\alias{Cvb}
\title{Cvb function...}
\usage{Cvb(xyt, spatial.intensity, N=100, spatial.covmodel, covpars)}
\description{Cvb function}
\details{This function is used in \code{thetaEst} to estimate the temporal correlation parameter, theta.}
\value{a function, see below.
Computes Monte carlo estimate of function C(v;beta) in Brix and Diggle 2001 pp 829 (... note later corrigendum to paper (2003) corrects the expression given in this paper)}
\references{\enumerate{
\item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
}}
\seealso{\link{thetaEst}}
\arguments{\item{xyt}{object of class stppp}
\item{spatial.intensity}{bivariate density estimate of lambda, an object of class im (produced from density.ppp for example)}
\item{N}{number of integration points}
\item{spatial.covmodel}{spatial covariance model}
\item{covpars}{additional covariance parameters}}
