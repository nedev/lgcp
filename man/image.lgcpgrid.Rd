\name{image.lgcpgrid}
\alias{image.lgcpgrid}
\title{image...}
\usage{\method{image}{lgcpgrid} (x, sel=1:x$len, ask=TRUE, ...)}
\description{image.lgcpgrid function}
\details{Produce an image plot of an lgcpgrid object.}
\value{grid plotting}
\seealso{\link{lgcpgrid.list}, \link{lgcpgrid.array}, \link{as.list.lgcpgrid}, \link{print.lgcpgrid}, 
\link{summary.lgcpgrid}, \link{quantile.lgcpgrid}, \link{plot.lgcpgrid}}
\arguments{\item{x}{an object of class lgcpgrid}
\item{sel}{vector of integers between 1 and grid$len: which grids to plot. Default NULL, in which case all grids are plotted.}
\item{ask}{logical; if TRUE the user is asked before each plot}
\item{...}{other arguments}}
