\name{showGrid.lgcpPredict}
\alias{showGrid.lgcpPredict}
\title{showGrid...}
\usage{\method{showGrid}{lgcpPredict} (x, ...)}
\description{showGrid.lgcpPredict function}
\details{This function displays the FFT grid used on a plot of an \code{lgcpPredict} object.
First plot the object using for example \code{plot(lg)}, where \code{lg} is an object
of class \code{lgcpPredict}, then for any of the plots produced, a call to
\code{showGrid(lg,pch=="+",cex=0.5)} will display the centroids of the FFT grid.}
\value{plots grid centroids on the current graphics device}
\seealso{\link{lgcpPredict}, \link{showGrid.default}, \link{showGrid.stppp}}
\arguments{\item{x}{an object of class lgcpPredict}
\item{...}{additional arguments  passed to points}}
