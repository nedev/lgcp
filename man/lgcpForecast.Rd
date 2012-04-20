\name{lgcpForecast}
\alias{lgcpForecast}
\title{lgcpForecast function...}
\usage{lgcpForecast(lg, ptimes, spatial.intensity, temporal.intensity)}
\description{lgcpForecast function}
\details{Function to produce forecasts for the mean field \eqn{Y}{Y} at times beyond the last time point in the 
analysis (given by the argument \code{T} in the function \code{lgcpPredict}).}
\value{forcasted relative risk, Poisson intensities and Y values over grid, together with approximate variance.}
\references{Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.}
\seealso{\link{lgcpPredict}}
\arguments{\item{lg}{an object of class lgcpPredict}
\item{ptimes}{vector of time points for prediction. Must start strictly after last inferred time point.}
\item{spatial.intensity}{the fixed spatial component: an object of that can be coerced to one of class spatialAtRisk}
\item{temporal.intensity}{the fixed temporal component: either a numeric vector, or a function that can be coerced into an object of class temporalAtRisk}}
