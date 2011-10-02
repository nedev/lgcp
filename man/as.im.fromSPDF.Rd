\name{as.im.fromSPDF}
\alias{as.im.fromSPDF}
\title{as...}
\usage{\method{as.im}{fromSPDF} (X, ncells=100, ...)}
\description{as.im.fromSPDF function}
\details{Convert an object of class fromSPDF (created by spatialAtRisk for example) into a spatstat im object.}
\value{object of class im containing normalised intensities}
\seealso{\link{as.im.fromXYZ}, \link{as.im.fromSPDF}, \link{as.im.fromFunction}, \link{as.fromXYZ}}
\arguments{\item{X}{an object of class fromSPDF}
\item{ncells}{number of cells to divide range into; default 100}
\item{...}{additional arguments}}
