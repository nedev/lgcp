\name{spatialAtRisk}
\alias{spatialAtRisk}
\title{spatialAtRisk function...}
\usage{spatialAtRisk(X, ...)}
\description{spatialAtRisk function}
\details{Generic function used in the construction of spatialAtRisk objects. The class of spatialAtRisk objects provide a framework for describing
the spatial inhomogeneity of the at-risk population, lambda(s). This is in contrast to the class of 
\link{temporalAtRisk} objects, which describe the global levels of the population at risk, mu(t).

Unless the user has specified lambda(s) directly by an R function (a mapping the from the real plane onto the non-negative real numbers, 
see ?spatialAtRisk.function), then it is only necessary to describe the population at risk up to a constant of proportionality, as the routines 
automatically normalise the lambda provided to integrate to 1.

For reference purposes, the following is a mathematical description of a log-Gaussian Cox Process, it is best viewed in the pdf version of the manual.

Let \eqn{\mathcal Y(s,t)}{\mathcal Y(s,t)} be a spatiotemporal Gaussian process, \eqn{W\subset R^2}{W\subset R^2} be an 
observation window in space and \eqn{T\subset R_{\geq 0}}{T\subset R_{\geq 0}} be an interval of time of interest. 
Cases occur at spatio-temporal positions \eqn{(x,t) \in W \times T}{(x,t) \in W \times T} 
according to an inhomogeneous spatio-temporal Cox process,
i.e. a Poisson process with a stochastic intensity \eqn{R(x,t)}{R(x,t)},
The number of cases, \eqn{X_{S,[t_1,t_2]}}{X_{S,[t_1,t_2]}}, arising in 
any \eqn{S \subseteq W}{S \subseteq W} during the interval \eqn{[t_1,t_2]\subseteq T}{[t_1,t_2]\subseteq T} is 
then Poisson distributed conditional on \eqn{R(\cdot)}{R(\cdot)},
\deqn{X_{S,[t_1,t_2]} \sim \mbox{Poisson}\left\{\int_S\int_{t_1}^{t_2} R(s,t)d sd t\right\}}{%
X_{S,[t_1,t_2]} \sim \mbox{Poisson}\left\{\int_S\int_{t_1}^{t_2} R(s,t)d sd t\right\}.}
Following Brix and Diggle (2001) and Diggle et al (2005), the intensity is decomposed multiplicatively as
\deqn{R(s,t) = \lambda(s)\mu(t)\exp\{\mathcal Y(s,t)\}.}{%
R(s,t) = \lambda(s)\mu(t)Exp\{\mathcal Y(s,t)\}.}
In the above, the fixed spatial component, \eqn{\lambda:R^2\mapsto R_{\geq 0}}{\lambda:R^2\mapsto R_{\geq 0}}, 
is a known function, proportional to the population at risk at each point in space and scaled so that
\deqn{\int_W\lambda(s)d s=1,}{%
\int_W\lambda(s)d s=1,}
whilst the fixed temporal component, 
\eqn{\mu:R_{\geq 0}\mapsto R_{\geq 0}}{\mu:R_{\geq 0}\mapsto R_{\geq 0}}, is also a known function with
\deqn{\mu(t) \delta t = E[X_{W,\delta t}],}{%
\mu(t) \delta t = E[X_{W,\delta t}],}
for \eqn{t}{t} in a small interval of time, \eqn{\delta t}{\delta t}, over which the rate of the process over \eqn{W}{W} can be considered constant.}
\value{method spatialAtRisk
\enumerate{
\item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
\item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
}}
\seealso{\link{lgcpPredict}, link{lgcpSim}, \link{spatialAtRisk.default}, \link{spatialAtRisk.fromXYZ}, \link{spatialAtRisk.im}, \link{spatialAtRisk.function}, \link{spatialAtRisk.SpatialGridDataFrame}, \link{spatialAtRisk.SpatialPolygonsDataFrame}, \link{spatialAtRisk.bivden}}
\arguments{\item{X}{an object}
\item{...}{additional arguments}}
