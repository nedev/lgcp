\name{lgcpgrid.matrix}
\alias{lgcpgrid.matrix}
\title{lgcpgrid...}
\usage{\method{lgcpgrid}{matrix} (grid, xvals=1:nrow(grid),
    yvals=1:ncol(grid), ...)}
\description{lgcpgrid.matrix function}
\details{Creates an lgcp grid object from an 2-dimensional matrix.}
\value{an object of class lgcpgrid}
\seealso{\link{lgcpgrid.list}, \link{as.list.lgcpgrid}, \link{print.lgcpgrid}, 
\link{summary.lgcpgrid}, \link{quantile.lgcpgrid}, \link{image.lgcpgrid}, \link{plot.lgcpgrid}}
\arguments{\item{grid}{a three dimensional array object}
\item{xvals}{optional vector of x-coordinates associated to grid. By default, this is the cell index in the x direction.}
\item{yvals}{optional vector of y-coordinates associated to grid. By default, this is the cell index in the y direction.}
\item{...}{other arguments}}
