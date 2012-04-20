\name{plotExceed.lgcpPredict}
\alias{plotExceed.lgcpPredict}
\title{plotExceed...}
\usage{\method{plotExceed}{lgcpPredict} (obj, fun, nlevel=64, ask=TRUE,
    plotcases=FALSE, ...)}
\description{plotExceed.lgcpPredict function}
\details{Function for plotting exceedance probabilities stored in \code{lgcpPredict} ojects.}
\value{plot of exceedances}
\seealso{\link{lgcpPredict}, \link{MonteCarloAverage}, \link{setoutput}}
\arguments{\item{obj}{an object}
\item{fun}{the name of the function used to compute exceedances (character vector of length 1). Note that the named function must be in memory.}
\item{nlevel}{number of colour levels to use in plot, default is 64}
\item{ask}{whether or not to ask for a new plot between plotting exceedances at different thresholds.}
\item{plotcases}{whether or not to plot the cases on the map}
\item{...}{additional arguments passed to image.plot}}
\examples{\dontrun{exceedfun <- exceedProbs(c(1.5,2,4))}
\dontrun{
plot(lg,"exceedfun") # lg is an object of class lgcpPredict
# in which the Monte Carlo mean of 
# "exceedfun" was computed  
# see ?MonteCarloAverage and ?setoutput
}}
