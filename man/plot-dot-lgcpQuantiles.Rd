\name{plot.lgcpQuantiles}
\alias{plot.lgcpQuantiles}
\title{plot...}
\usage{\method{plot}{lgcpQuantiles} (x, sel=1:dim(x)[3], ask=TRUE, crop=TRUE,
    plotwin=FALSE, ...)}
\description{plot.lgcpQuantiles function}
\details{Plots \code{lgcpQuantiles} objects: output from \code{quantiles.lgcpPredict}}
\value{grid plotting
This is a wrapper function for image.lgcpgrid}
\seealso{\link{quantile.lgcpPredict}}
\arguments{\item{x}{an object of class lgcpQuantiles}
\item{sel}{vector of integers between 1 and grid$len: which grids to plot. Default NULL, in which case all grids are plotted.}
\item{ask}{logical; if TRUE the user is asked before each plot}
\item{crop}{whether or not to crop to bounding box of observation window}
\item{plotwin}{logical whether to plot the window attr(x,"window"), default is FALSE}
\item{...}{other arguments  passed to image.plot}}
\examples{\dontrun{qtiles <- quantile(lg,qt=c(0.5,0.75,0.9),fun=exp)} 
# assumed that lg has class lgcpPredict
\dontrun{plot(qtiles)}}
