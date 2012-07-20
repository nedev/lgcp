\name{autocorr}
\alias{autocorr}
\title{autocorr function...}
\usage{autocorr(x, lags, tidx, inWindow=x$xyt$window, crop2parentwindow=TRUE,
    ...)}
\description{autocorr function}
\details{\bold{This function requires data to have been dumped to disk}: see \code{?dump2dir} and \code{?setoutput}. The routine \code{autocorr.lgcpPredict} 
computes cellwise selected autocorrelations of Y.
Since computing the quantiles is an expensive operation, the option to output the quantiles on a subregion of interest is also provided (by
setting the argument \code{inWindow}, which has a sensible default).}
\value{an array, the [,,i]th slice being the grid of cell-wise autocorrelations.}
\seealso{\link{lgcpPredict}, \link{dump2dir}, \link{setoutput}, \link{plot.lgcpAutocorr}}
\arguments{\item{x}{an object of class lgcpPredict}
\item{lags}{a vector of the required lags}
\item{tidx}{the index number of the the time interval of interest, default is the last time point.}
\item{inWindow}{an observation owin window on which to compute the autocorrelations, can speed up calculation. Default is x$xyt$window, set to NULL for full grid.}
\item{crop2parentwindow}{logical: whether to only compute autocorrelations for cells inside x$xyt$window (the 'parent window')}
\item{...}{additional arguments}}
