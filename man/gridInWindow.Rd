\name{gridInWindow}
\alias{gridInWindow}
\title{gridInWindow function...}
\usage{gridInWindow(xvals, yvals, win)}
\description{gridInWindow function}
\details{For the grid defined by x-coordinates, xvals, and y-coordinates, yvals, and an owin object W, this function just returns
a logical matrix M, whose [i,j] entry is TRUE if the point(xvals[i], yvals[j]) is inside the observation window.}
\value{matrix of TRUE/FALSE, which elements of the grid are inside the observation window win}
\arguments{\item{xvals}{x coordinates}
\item{yvals}{y coordinates}
\item{win}{owin object}}
