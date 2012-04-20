\name{constantInTime.numeric}
\alias{constantInTime.numeric}
\title{constantInTime...}
\usage{\method{constantInTime}{numeric} (obj, tlim, warn=TRUE, ...)}
\description{constantInTime.numeric function}
\details{Create a constant-in-time temporalAtRisk object from a numeric object of length 1. The returned temporalAtRisk object is assumed to have
been scaled correctly by the user so that mu(t) = E(number of cases in a unit time interval).

For further details of temporalAtRisk objects, see ?temporalAtRisk>}
\value{a function f(t) giving the (constant) temporal intensity at time t for integer t in the interval [tlim[1],tlim[2]] of class temporalAtRisk}
\seealso{\link{temporalAtRisk}, \link{spatialAtRisk}, \link{temporalAtRisk.numeric}, \link{temporalAtRisk.function}, \link{constantInTime}, \link{constantInTime.stppp}, \link{print.temporalAtRisk}, \link{plot.temporalAtRisk},}
\arguments{\item{obj}{numeric constant}
\item{tlim}{vector of length 2 giving time limits}
\item{warn}{Issue a warning if the given temporal intensity treated is treated as 'known'?}
\item{...}{additional arguments}}
