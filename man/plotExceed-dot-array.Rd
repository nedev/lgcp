\name{plotExceed.array}
\alias{plotExceed.array}
\title{plotExceed...}
\usage{\method{plotExceed}{array} (obj, fun, lgcppredict, xvals, yvals,
    window, cases, nlevel=64, ask=TRUE, mapunderlay, alpha=1, ...)}
\description{plotExceed.array function}
\details{Function for plotting exceedance probabilities stored in array objects. Used in \code{plotExceed.lgcpPredict}.}
\value{generic function returning method plotExceed}
\seealso{\link{plotExceed.lgcpPredict}}
\arguments{\item{obj}{an object}
\item{fun}{the name of the function used to compute exceedances (character vector of length 1). Note that the named function must be in memory.}
\item{lgcppredict}{an object of class lgcpPredict that can be used to supply an observation window and x and y coordinates}
\item{xvals}{optional vector giving x coords of centroids of cells}
\item{yvals}{optional vector giving y coords of centroids of cells}
\item{window}{optional obervation window}
\item{cases}{optional xy (n x 2) matrix of locations of cases to plot}
\item{nlevel}{number of colour levels to use in plot, default is 64}
\item{ask}{whether or not to ask for a new plot between plotting exceedances at different thresholds.}
\item{mapunderlay}{optional underlay to plot underneath maps of exceedance probabilities. Use in conjunction with rainbow parameter 'alpha' (eg alpha=0.3) to set transparency of exceedance layer.}
\item{alpha}{graphical parameter takign values in [0,1] controlling transparency of exceedance layer. Default is 1.}
\item{...}{additional arguments passed to image.plot}}
