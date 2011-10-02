\name{lgcpSim}
\alias{lgcpSim}
\title{lgcpsim function...}
\usage{lgcpSim(owin, tlim=as.integer(c(0, 10)), spatial.intensity,
    temporal.intensity, cellwidth=0.05, model.parameters=lgcppars(sigma
    = 2, phi = 0.2, theta = 3), spatial.covmodel="exponential",
    covpars=c(), progressbar=TRUE, returnintensities=FALSE, plot=FALSE)}
\description{lgcpsim function}
\details{Approximate simulation from a spatiotemoporal log-Gaussian Cox Process. Returns an stppp object.

The following is a mathematical description of a log-Gaussian Cox Process, it is best viewed in the pdf version of the manual.

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
\value{returns an object of class stppp containing the simulated data}
\references{\enumerate{
\item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
\item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
\item Wood ATA, Chan G (1994). Simulation of Stationary Gaussian Processes in [0,1]d. Journal of Computational and Graphical Statistics, 3(4), 409-432.
\item Moller J, Syversveen AR, Waagepetersen RP (1998). Log Gaussian Cox Processes. Scandinavian Journal of Statistics, 25(3), 451-482.
}}
\seealso{\link{lgcpPredict}, \link{showGrid.stppp}, \link{stppp}}
\arguments{\item{owin}{a spatstat owin object, possibly with polygonal boundary}
\item{tlim}{time interval on which to simulate data}
\item{spatial.intensity}{object that can be coerced into a spatialAtRisk object. if NULL then uniform spatial is chosen}
\item{temporal.intensity}{the fixed temporal component: either a numeric vector, or a function that can be coerced into an object of class temporalAtRisk}
\item{cellwidth}{size of gris discretisation as in lgcpPredict}
\item{model.parameters}{parameters of hte Gaussian process, specify using function lgcppars,}
\item{spatial.covmodel}{spatial covariance model, default is "exponential"}
\item{covpars}{additional covariance parameters}
\item{progressbar}{logical, whether to print a progress bar. Default TRUE.}
\item{returnintensities}{logigal, whether to return the vector of simulation times and corresponding array of spatial intensities at each time. Default FALSE.}
\item{plot}{plot each time step? Only really useful for debugging purposes.}}
\examples{xyt <- lgcpSim()}
