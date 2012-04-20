\name{as.fromXYZ.fromFunction}
\alias{as.fromXYZ.fromFunction}
\title{as...}
\usage{\method{as.fromXYZ}{fromFunction} (X, xyt, M=100, N=100, ...)}
\description{as.fromXYZ.fromFunction function}
\details{Method for converting from the fromFunction class of objects to the fromXYZ class of objects. Clearly this requires the
user to specify a grid onto which to compute the discretised verion.}
\value{object of class im containing normalised intensities}
\seealso{\link{as.im.fromXYZ}, \link{as.im.fromSPDF}, \link{as.im.fromFunction}, \link{as.fromXYZ}}
\arguments{\item{X}{an object of class fromFunction}
\item{xyt}{and objects of class stppp}
\item{M}{number of cells in x direction}
\item{N}{number of cells in y direction}
\item{...}{additional arguments}}
