\name{plot.lgcpgrid}
\alias{plot.lgcpgrid}
\title{plot...}
\usage{\method{plot}{lgcpgrid} (x, sel=1:x$len, ask=TRUE, ...)}
\description{plot.lgcpgrid function}
\details{This is a wrapper function for image.lgcpgrid}
\value{an image-type plot}
\seealso{\link{lgcpgrid.list}, \link{lgcpgrid.array}, \link{as.list.lgcpgrid}, \link{print.lgcpgrid}, 
\link{summary.lgcpgrid},\link{quantile.lgcpgrid}, \link{image.lgcpgrid}}
\arguments{\item{x}{an object of class lgcpgrid}
\item{sel}{vector of integers between 1 and grid$len: which grids to plot. Default NULL, in which case all grids are plotted.}
\item{ask}{logical; if TRUE the user is asked before each plot}
\item{...}{other arguments}}
