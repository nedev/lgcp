\name{lgcpgrid.list}
\alias{lgcpgrid.list}
\title{lgcpgrid...}
\usage{\method{lgcpgrid}{list} (grid, ...)}
\description{lgcpgrid.list function}
\details{Creates an lgcpgrid object from a list object. Note that each element of the list should be a matrix,
and that each matrix should have the same dimension.}
\value{an object of class lgcpgrid}
\seealso{\link{lgcpgrid.array}, \link{as.list.lgcpgrid}, \link{print.lgcpgrid}, 
\link{summary.lgcpgrid}, \link{quantile.lgcpgrid}, \link{image.lgcpgrid}, \link{plot.lgcpgrid}}
\arguments{\item{grid}{a list object with each member of the list being a numeric matrix, each matrix having the same dimension}
\item{...}{other arguments}}
