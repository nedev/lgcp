\name{getRotation.stppp}
\alias{getRotation.stppp}
\title{getRotation...}
\usage{\method{getRotation}{stppp} (xyt, ...)}
\description{getRotation.stppp function}
\details{Compute  rotation matrix if observation window is a polygonal boundary}
\value{the optimal rotation matrix and rotated data and observation window. Note it may or may not be advantageous to rotate the window, this information is displayed prior to the MALA routine when using lgcpPredict}
\arguments{\item{xyt}{an object of class stppp}
\item{...}{additional arguments}}
