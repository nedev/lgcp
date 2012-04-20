\name{plot.lgcpPredict}
\alias{plot.lgcpPredict}
\title{plot...}
\usage{\method{plot}{lgcpPredict} (x, type="relrisk", sel=1:x$EY.mean$len,
    plotdata=TRUE, ask=TRUE, clipWindow=TRUE, ...)}
\description{plot.lgcpPredict function}
\details{Simple plotting function for objects of class \code{lgcpPredict}.}
\value{plots the Monte Carlo mean of quantities obtained via simulation. By default the mean relative risk is plotted.}
\seealso{\link{lgcpPredict}}
\arguments{\item{x}{an object of class lgcpPredict}
\item{type}{Character string: what type of plot to produce. Choices are "relrisk" (=exp(Y)); "serr" (standard error of relative risk); or "intensity" (=lambda*mu*exp(Y)).}
\item{sel}{vector of integers between 1 and grid$len: which grids to plot. Default NULL, in which case all grids are plotted.}
\item{plotdata}{whether or not to overlay the data}
\item{ask}{logical; if TRUE the user is asked before each plot}
\item{clipWindow}{whether to plot grid cells outside the observation window}
\item{...}{additional arguments passed to image.plot}}
