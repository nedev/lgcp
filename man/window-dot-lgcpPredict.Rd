\name{window.lgcpPredict}
\alias{window.lgcpPredict}
\title{window...}
\usage{\method{window}{lgcpPredict} (x, ...)}
\description{window.lgcpPredict function}
\details{Accessor function returning the observation window from objects of class \code{lgcpPredict}. Note that for
computational purposes, the window of an \code{stppp} object will be extended to accommodate the requirement that 
the dimensions must be powers of 2. The function \code{window.lgcpPredict} returns the extended window.}
\value{returns the observation window used durign computation}
\seealso{\link{lgcpPredict}}
\arguments{\item{x}{an object of class lgcpPredict}
\item{...}{additional arguments}}
