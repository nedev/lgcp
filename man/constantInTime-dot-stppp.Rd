\name{constantInTime.stppp}
\alias{constantInTime.stppp}
\title{constantInTime...}
\usage{\method{constantInTime}{stppp} (obj, ...)}
\description{constantInTime.stppp function}
\details{Create a constant-in-time temporalAtRisk object from an stppp object. The returned temporalAtRisk object is scaled to return
mu(t) = E(number of cases in a unit time interval).

For further details of temporalAtRisk objects, see ?temporalAtRisk>}
\value{a function f(t) giving the (constant) temporal intensity at time t for integer t in the interval [tlim[1],tlim[2]] of class temporalAtRisk}
\seealso{\link{temporalAtRisk}, \link{spatialAtRisk}, \link{temporalAtRisk.numeric}, \link{temporalAtRisk.function}, \link{constantInTime}, \link{constantInTime.numeric}, \link{print.temporalAtRisk}, \link{plot.temporalAtRisk},}
\arguments{\item{obj}{an object of class stppp.}
\item{...}{additional arguments}}
