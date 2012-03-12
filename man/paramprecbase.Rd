\name{paramprecbase}
\alias{paramprecbase}
\title{paramprecbase function...}
\usage{paramprecbase(ns, M, N, inverse=FALSE)}
\description{paramprecbase function}
\details{A function to compute the parametrised base matrix of a precision matrix of a GMRF on an M x N toral grid with neighbourhood size ns. Note that 
the precision matrix is block circulant. The returned function operates on a parameter vector as in Rue and Held (2005) pp 187.}
\value{a functioin that returns the base matrix of the precision matrix}
\arguments{\item{ns}{neighbourhood size}
\item{M}{number of x cells}
\item{N}{number of y cells}
\item{inverse}{whether or not to compute the base matrix of the inverse precision matrix (ie the covariance matrix). default is FALSE}}
